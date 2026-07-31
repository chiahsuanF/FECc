#!/usr/bin/env python3

import os
import sys
import logging

from fec_login import *


# -------------------------------------------------------
# 6. Simple-cycle detection
# -------------------------------------------------------

def normalize_cycle(cycle):
    nodes = list(cycle[:-1])

    if not nodes:
        return tuple()

    # Forward orientation, rotated to its smallest node.
    min_node = min(nodes)
    min_index = nodes.index(min_node)
    forward = nodes[min_index:] + nodes[:min_index]

    # Reverse orientation, also rotated to its smallest node.
    reverse_nodes = list(reversed(nodes))
    min_index = reverse_nodes.index(min_node)
    backward = (
        reverse_nodes[min_index:]
        + reverse_nodes[:min_index]
    )

    return tuple(min(forward, backward))


def dfs_visit(start, node, graph, path, cycles):
    path.append(node)

    for neighbor in graph.get(node, {}):
        if neighbor == start and len(path) > 1:
            cycles.append(path.copy() + [start])

        elif neighbor not in path:
            dfs_visit(
                start,
                neighbor,
                graph,
                path,
                cycles
            )

    path.pop()


def find_cycles(graph, outdir, sample):
    t0 = tic("Detecting cycles")

    detected_cycles = []

    for start_node in graph:
        dfs_visit(
            start_node,
            start_node,
            graph,
            [],
            detected_cycles
        )

    unique_cycles = {}

    for cycle in detected_cycles:
        key = normalize_cycle(cycle)

        if key and key not in unique_cycles:
            normalized_nodes = list(key)
            unique_cycles[key] = (
                normalized_nodes + [normalized_nodes[0]]
            )

    cycles = list(unique_cycles.values())

    output_file = os.path.join(
        outdir,
        f"06_{sample}_raw_cycles.txt"
    )

    with open(output_file, "w") as handle:
        handle.write("cycle_nodes\n")

        for cycle in cycles:
            handle.write(
                ",".join(map(str, cycle)) + "\n"
            )

    logging.info(
        "Wrote %d unique raw cycles -> %s",
        len(cycles),
        output_file
    )

    toc(t0, "Detecting cycles")

    return cycles


# -------------------------------------------------------
# Cycle and event sizes
# -------------------------------------------------------

def cycle_size(cycle, regions):
    total_size = 0

    for node in cycle[:-1]:
        chrom, start, end = regions[node]
        total_size += end - start + 1

    return total_size


def event_size(nodes, regions):
    return sum(
        regions[node][2] - regions[node][1] + 1
        for node in set(nodes)
    )


# -------------------------------------------------------
# Filter subset cycles
# -------------------------------------------------------

def filter_maximal_cycles(cycles):
    if not cycles:
        return []

    node_sets = [
        set(cycle[:-1])
        for cycle in cycles
    ]

    keep = []

    for index, current_nodes in enumerate(node_sets):
        is_strict_subset = any(
            index != other_index
            and current_nodes < other_nodes
            for other_index, other_nodes in enumerate(node_sets)
        )

        if not is_strict_subset:
            keep.append(cycles[index])

    return keep


# -------------------------------------------------------
# Merge overlapping cycles into eccDNA events
# -------------------------------------------------------

def merge_overlapping_cycles(cycles):
    cycle_sets = [
        set(cycle[:-1])
        for cycle in cycles
        if len(cycle) > 1
    ]

    merged_events = []

    while cycle_sets:
        event_nodes = cycle_sets.pop()

        changed = True

        while changed:
            changed = False
            remaining_sets = []

            for candidate_nodes in cycle_sets:
                if event_nodes.intersection(candidate_nodes):
                    event_nodes.update(candidate_nodes)
                    changed = True
                else:
                    remaining_sets.append(candidate_nodes)

            cycle_sets = remaining_sets

        merged_events.append(event_nodes)

    return merged_events

def expand_events_with_attached_nodes( events, graph, max_hops=1, min_edge_support=1):

    expanded_events = []

    for event_nodes in events:
        event_nodes = set(event_nodes)
        frontier = set(event_nodes)

        for _ in range(max_hops):
            new_nodes = set()

            for node in frontier:
                for neighbor, edge_data in graph.get(node, {}).items():

                    # Support dictionary, numeric support, or other
                    # existing graph edge representations.
                    if isinstance(edge_data, dict):
                        support = edge_data.get(
                            "support",
                            edge_data.get("count", 1)
                        )
                    elif isinstance(edge_data, (int, float)):
                        support = edge_data
                    else:
                        support = 1

                    if (
                        support >= min_edge_support
                        and neighbor not in event_nodes
                    ):
                        new_nodes.add(neighbor)

            if not new_nodes:
                break

            event_nodes.update(new_nodes)
            frontier = new_nodes

        expanded_events.append(event_nodes)

    # Expansion can make two previously separate events overlap.
    # Merge them again transitively.
    pending = list(expanded_events)
    final_events = []

    while pending:
        current_event = pending.pop()
        changed = True

        while changed:
            changed = False
            remaining = []

            for candidate_event in pending:
                if current_event.intersection(candidate_event):
                    current_event.update(candidate_event)
                    changed = True
                else:
                    remaining.append(candidate_event)

            pending = remaining

        final_events.append(current_event)

    return final_events

# -------------------------------------------------------
# Coordinate ordering
# -------------------------------------------------------

def order_event_nodes(event_nodes, cycles):
    event_nodes = set(event_nodes)
    ordered_nodes = []
    seen = set()

    for cycle in cycles:
        cycle_nodes = set(cycle[:-1])

        if not cycle_nodes.intersection(event_nodes):
            continue

        for node in cycle[:-1]:
            if node in event_nodes and node not in seen:
                ordered_nodes.append(node)
                seen.add(node)

    for node in sorted(event_nodes):
        if node not in seen:
            ordered_nodes.append(node)

    return ordered_nodes


# -------------------------------------------------------
# Write individual simple cycles
# -------------------------------------------------------

def write_cycles_with_coords( cycles, regions, outdir, sample, label):

    maximal_cycles = filter_maximal_cycles(cycles)

    if label == "filtered":
        output_file = os.path.join(
            outdir,
            "summary.txt"
        )
    else:
        output_file = os.path.join(
            outdir,
            f"07_{sample}_{label}_cycles_coords.txt"
        )

    with open(output_file, "w") as handle:
        handle.write(
            "cycle_id\tsize\tfragment_count\tpath\n"
        )

        for cycle_id, cycle in enumerate(
            maximal_cycles,
            start=1
        ):
            size = cycle_size(cycle, regions)

            path = ";".join(
                (
                    f"{regions[node][0]}:"
                    f"{regions[node][1]}-"
                    f"{regions[node][2]}"
                )
                for node in cycle
            )

            handle.write(
                f"{cycle_id}\t"
                f"{size}\t"
                f"{len(cycle) - 1}\t"
                f"{path}\n"
            )

    logging.info(
        "Wrote %d %s simple cycles -> %s",
        len(maximal_cycles),
        label,
        output_file
    )


# -------------------------------------------------------
# Write merged eccDNA events
# -------------------------------------------------------
def write_merged_eccdna_events( cycles, regions, graph, outdir, sample,
    label,ninclude_attached=True, max_hops=1, min_edge_support=1):

    events = merge_overlapping_cycles(cycles)

    if include_attached:
        events = expand_events_with_attached_nodes(
            events=events,
            graph=graph,
            max_hops=max_hops,
            min_edge_support=min_edge_support
        )

    # Stable event ordering: largest event first.
    events.sort(
        key=lambda nodes: (
            -event_size(nodes, regions),
            min(nodes)
        )
    )

    if label == "filtered":
        output_file = os.path.join(
            outdir,
            "summary.txt"
        )
    else:
        output_file = os.path.join(
            outdir,
            f"08_{sample}_{label}_merged_events.txt"
        )

    with open(output_file, "w") as handle:
        handle.write(
            "event_id\t"
            "total_size\t"
            "fragment_count\t"
            "coordinates\n"
        )

        for event_id, event_nodes in enumerate(
            events,
            start=1
        ):
            ordered_nodes = order_event_nodes(
                event_nodes,
                cycles
            )

            total_size = event_size(
                event_nodes,
                regions
            )

            coordinates = ";".join(
                f"{regions[node][0]}:"
                f"{regions[node][1]}-"
                f"{regions[node][2]}"
                for node in ordered_nodes
            )

            handle.write(
                f"eccDNA_{sample}_{event_id}\t"
                f"{total_size}\t"
                f"{len(event_nodes)}\t"
                f"{coordinates}\n"
            )

    logging.info(
        "Wrote %d merged %s eccDNA events -> %s",
        len(events),
        label,
        output_file
    )

    return events
# -------------------------------------------------------
# Diagnostics
# -------------------------------------------------------

def log_cycle_diagnostics(cycles, regions, graph):
    cyclic_nodes = set()

    for cycle_id, cycle in enumerate(cycles, start=1):
        nodes = cycle[:-1]
        cyclic_nodes.update(nodes)

        coordinates = [
            (
                f"{regions[node][0]}:"
                f"{regions[node][1]}-"
                f"{regions[node][2]}"
            )
            for node in nodes
        ]

        logging.info(
            "Cycle %d: nodes=%d size=%d path=%s",
            cycle_id,
            len(nodes),
            cycle_size(cycle, regions),
            " -> ".join(coordinates)
        )

    graph_nodes = set(graph)
    noncyclic_nodes = graph_nodes - cyclic_nodes

    logging.info(
        "Detected %d unique nodes in cycles",
        len(cyclic_nodes)
    )

    logging.info(
        "Found %d graph nodes not included in any cycle",
        len(noncyclic_nodes)
    )

    for node in sorted(noncyclic_nodes):
        if node < 0 or node >= len(regions):
            logging.warning(
                "Graph node %s is outside the regions array",
                node
            )
            continue

        chrom, start, end = regions[node]

        logging.info(
            "Non-cyclic node %s: %s:%s-%s; neighbors=%s",
            node,
            chrom,
            start,
            end,
            list(graph.get(node, {}).keys())
        )
