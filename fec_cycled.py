#!/usr/bin/env python3

import os
import logging

from fec_login import *


# -------------------------------------------------------
# 6. Simple-cycle detection
# -------------------------------------------------------


def normalize_cycle(cycle):
    """
    Canonicalize a cycle.

    Removes:
      1. rotation differences
      2. reverse-direction differences
      3. strand differences

    Keeps only region IDs.
    """

    nodes = [node[0] for node in cycle[:-1]]

    if not nodes:
        return tuple()

    n = len(nodes)

    representations = []

    #
    # forward rotations
    #
    for i in range(n):
        representations.append(
            tuple(nodes[i:] + nodes[:i])
        )

    #
    # reverse rotations
    #
    rev = list(reversed(nodes))

    for i in range(n):
        representations.append(
            tuple(rev[i:] + rev[:i])
        )

    return min(representations)


def cycle_support(cycle, graph):
    """
    Total split-read support of the junctions a cycle traverses.
    """
    return sum(
        graph.get(left, {}).get(right, 0)
        for left, right in zip(cycle, cycle[1:])
    )


def dfs_visit(
    start,
    node,
    graph,
    path,
    path_segment_ids,
    cycles,
    max_cycle_len
):

    path.append(node)
    path_segment_ids.add(node[0])

    for neighbor in graph.get(node, {}):

        if neighbor == start and len(path) > 1: #1

            full_cycle = path + [start]

            key = normalize_cycle(
                full_cycle
            )

            # The same regions can be traversed in several orientation
            # patterns, each using different junctions. normalize_cycle
            # gives them all one key, so keep the best-supported one
            # rather than whichever the search reached first.
            candidate = (
                cycle_support(full_cycle, graph),
                tuple(full_cycle)
            )

            best = cycles.get(key)

            if best is None or candidate > best[0]:
                cycles[key] = (candidate, full_cycle)

        elif (
            len(path) < max_cycle_len
            and neighbor[0] not in path_segment_ids
            #and neighbor >= start
        ):

            dfs_visit(
                start,
                neighbor,
                graph,
                path,
                path_segment_ids,
                cycles, max_cycle_len
            )

    path.pop()
    path_segment_ids.remove(node[0])


def find_cycles(graph, outdir, sample, max_cycle_len=20):

    t0 = tic("Detecting cycles")

    logging.info(
        f"Cycle graph: "
        f"{len(graph):,} nodes, "
        f"{sum(len(v) for v in graph.values()):,} edges"
    )

    detected_cycles = {}
    starts = [ n for n, nbrs in graph.items() if nbrs]

    for i, start_node in enumerate(starts):

        if i % 100 == 0:
            logging.info(
                f"Cycle search progress "
                f"{i:,}/{len(starts):,}"
            )

        dfs_visit(
            start_node,
            start_node,
            graph,
            [],
            set(),
            detected_cycles, max_cycle_len
        )

    cycles = [
        cycle
        for _, cycle in detected_cycles.values()
    ]


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
        f"Wrote {len(cycles):,} unique raw cycles -> "
        f"{output_file}"
    )

    toc(t0, "Detecting cycles")

    return cycles

# -------------------------------------------------------
# Cycle and event sizes
# -------------------------------------------------------

def cycle_size(cycle, regions):
    """
    Return the total base-pair size of a simple cycle.

    The final node is excluded because it repeats the first node.
    """
    total_size = 0

    for node in cycle[:-1]:
        node_id = node[0]
        chrom, start, end = regions[node_id]
        total_size += end - start + 1

    return total_size


def event_size(nodes, regions):
    """
    Return the total base-pair size of an eccDNA event.

    Each region node is counted exactly once.
    """
    return sum(
        regions[node[0]][2] - regions[node[0]][1] + 1
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

def expand_events_with_attached_nodes(
    events,
    graph,
    max_hops=1,
    min_edge_support=1
):

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
def write_cycles_with_coords(
    cycles,
    regions,
    outdir,
    sample,
    label
):
    """
    Write simple cycles with genomic coordinates.

    Raw output can be restricted to maximal cycles.
    Filtered output preserves every quality-filtered cycle, including
    valid two-region cycles that are subsets of larger cycles.
    """

    if label == "filtered":
        output_cycles = cycles

        output_file = os.path.join(
            outdir,
            "summary.txt"
        )
    else:
        output_cycles = filter_maximal_cycles(cycles)

        output_file = os.path.join(
            outdir,
            f"07_{sample}_{label}_cycles_coords.txt"
        )

    with open(output_file, "w") as handle:

        handle.write(
            "cycle_id\t"
            "size\t"
            "fragment_count\t"
            "path\n"
        )

        for cycle_id, cycle in enumerate(
            output_cycles,
            start=1
        ):
            size = cycle_size(
                cycle,
                regions
            )

            path = ";".join(
                (
                    f"{regions[node[0]][0]}:"
                    f"{regions[node[0]][1]}-"
                    f"{regions[node[0]][2]}"
                    f"{node[1]}"
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
        len(output_cycles),
        label,
        output_file
    )

    return output_cycles
# -------------------------------------------------------
# Write merged eccDNA events
# -------------------------------------------------------

def write_merged_eccdna_events(
    cycles,
    regions,
    graph,
    outdir,
    sample,
    label,
    include_attached=True,
    max_hops=1,
    min_edge_support=1
):
    """
    Merge overlapping cycles into eccDNA events.

    If include_attached=True, graph nodes directly connected to a
    cyclic event are also included, even if those nodes are not part
    of a simple cycle.
    """
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
                f"{regions[node[0]][0]}:"
                f"{regions[node[0]][1]}-"
                f"{regions[node[0]][2]}"
                f"{node[1]}"
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
    """
    Log the coordinates in every cycle and report graph nodes that are
    not included in any detected cycle.
    """
    cyclic_nodes = set()

    for cycle_id, cycle in enumerate(cycles, start=1):
        nodes = cycle[:-1]
        cyclic_nodes.update(nodes)

        coordinates = [
            (
                f"{regions[node[0]][0]}:"
                f"{regions[node[0]][1]}-"
                f"{regions[node[0]][2]}"
                f"{node[1]}"
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
        node_id = node[0]
        if node_id < 0 or node_id >= len(regions):
            logging.warning(
                "Graph node %s is outside the regions array",
                node
            )
            continue

        chrom, start, end = regions[node_id]

        logging.info(
            "Non-cyclic node %s: %s:%s-%s; neighbors=%s",
            node,
            chrom,
            start,
            end,
            list(graph.get(node, {}).keys())
        )