#!/usr/bin/env python3

import logging
import pysam
from collections import Counter, defaultdict

def parse_nodes(path):
    """
    Parse a semicolon-separated cycle path.

    Example:
        A+;B-;C+;A+

    Returns:
        ["A+", "B-", "C+", "A+"]
    """
    return [
        node.strip()
        for node in path.strip().split(";")
        if node.strip()
    ]


def strandless_node(node):
    """
    Return a region without its terminal strand symbol.

    Example:
        chr2:100-200+ -> chr2:100-200
        chr2:100-200- -> chr2:100-200
    """
    if node.endswith("+") or node.endswith("-"):
        return node[:-1]

    return node


def remove_one(values, target):
    """
    Remove one occurrence of target from a list.

    Returns True if an occurrence was removed.
    """
    try:
        values.remove(target)
        return True
    except ValueError:
        return False


def remove_leaf_nodes(graph):
    """
    Iteratively remove nodes with degree zero or one.

    The graph is an undirected multigraph represented by adjacency
    lists. Duplicate neighbors represent repeated edges.
    """
    while True:
        leaves = [
            node
            for node, neighbors in graph.items()
            if len(neighbors) <= 1
        ]

        if not leaves:
            break

        for leaf in leaves:
            if leaf not in graph:
                continue

            neighbors = list(graph[leaf])

            for neighbor in neighbors:
                if neighbor in graph:
                    remove_one(graph[neighbor], leaf)

            del graph[leaf]

    return graph


def graph_has_edges(graph):
    """
    Return True if the graph contains at least one edge.
    """
    return any(neighbors for neighbors in graph.values())


def euler_walk(graph):
    """
    Run Hierholzer's algorithm on an undirected multigraph.

    Duplicate edges are retained, so the resulting walk may revisit
    the same region multiple times.

    Returns:
        Ordered node list. For an Eulerian cycle, the first and final
        nodes are identical.
    """
    work_graph = {
        node: list(neighbors)
        for node, neighbors in graph.items()
    }

    if not graph_has_edges(work_graph):
        return []

    # Prefer a node with odd degree for an Eulerian path.
    # If all degrees are even, this is an Eulerian cycle and any node
    # with an edge can be used.
    odd_nodes = [
        node
        for node, neighbors in work_graph.items()
        if len(neighbors) % 2 == 1
    ]

    if odd_nodes:
        start = sorted(odd_nodes)[0]
    else:
        start = sorted(
            node
            for node, neighbors in work_graph.items()
            if neighbors
        )[0]

    stack = [start]
    walk = []

    while stack:
        current = stack[-1]

        if work_graph.get(current):
            neighbor = work_graph[current].pop()

            # Remove the matching reverse copy of this undirected edge.
            if neighbor not in work_graph:
                logging.warning(
                    "Missing graph node while traversing edge %s -> %s",
                    current,
                    neighbor
                )
                continue

            if not remove_one(work_graph[neighbor], current):
                logging.warning(
                    "Missing reverse edge while traversing %s -> %s",
                    current,
                    neighbor
                )
                continue

            stack.append(neighbor)

        else:
            walk.append(stack.pop())

    walk.reverse()

    # Confirm that every edge in this connected graph was consumed.
    remaining_edges = sum(
        len(neighbors)
        for neighbors in work_graph.values()
    )

    if remaining_edges:
        logging.warning(
            "Euler traversal left %d adjacency entries unused",
            remaining_edges
        )

    return walk


def connected_edge_components(graph):
    """
    Split a graph into connected components containing edges.

    This protects against a strandless cycle group containing multiple
    disconnected strand-aware subgraphs.
    """
    components = []
    visited = set()

    for start in sorted(graph):

        if start in visited or not graph[start]: continue

        component_nodes = set()
        stack = [start]
        visited.add(start)


        while stack:
            node = stack.pop()
            component_nodes.add(node)

            for neighbor in graph.get(node, []):
                if neighbor not in visited:
                    visited.add(neighbor)
                    stack.append(neighbor)

        component_graph = defaultdict(list)

        for node in component_nodes:
            component_graph[node] = [
                neighbor
                for neighbor in graph.get(node, [])
                if neighbor in component_nodes
            ]

        components.append(component_graph)

    return components


def path_size(path):
    """
    Calculate the summed fragment length of an Euler walk.

    The final node is excluded when it repeats the first node.
    Repeated internal nodes are counted repeatedly because they
    represent repeated traversal/copy usage.
    """
    if not path:
        return 0

    nodes = path[:-1] if path[0] == path[-1] else path

    total = 0

    for node in nodes:
        region = strandless_node(node)

        try:
            coordinate_text = region.rsplit(":", 1)[1]
            start_text, end_text = coordinate_text.split("-", 1)

            start = int(start_text)
            end = int(end_text)

            total += end - start + 1

        except (IndexError, ValueError):
            logging.warning(
                "Could not parse region size from node: %s",
                node
            )

    return total

def parse_region(node):
    """
    Parse a strand-aware node.

    Example:
        chr2:11551952-12459661+

    Returns:
        ("chr2", 11551952, 12459661)
    """
    region = strandless_node(node)

    chrom, coordinates = region.rsplit(":", 1)
    start_text, end_text = coordinates.split("-", 1)

    return (
        chrom,
        int(start_text),
        int(end_text)
    )


def calculate_region_depth(
    bam,
    node,
    min_mapping_quality=0
):
    """
    Calculate average read depth across a genomic region.

    Zero-coverage positions are included because total aligned depth
    is divided by the full region length.
    """
    chrom, start, end = parse_region(node)

    region_length = end - start + 1

    if region_length <= 0:
        return 0.0

    total_depth = 0

    try:
        for column in bam.pileup(
            chrom,
            start,
            end,
            truncate=True,
            stepper="samtools",
            min_mapping_quality=min_mapping_quality,
            max_depth=1000000
        ):
            total_depth += column.nsegments

    except ValueError:
        logging.warning(
            "Could not calculate depth for %s",
            node
        )
        return 0.0

    return total_depth / float(region_length)


def calculate_copy_number(
    region_depth,
    baseline_depth
):
    """
    Estimate integer copy number from region depth.

    Examples:
        region depth 18, baseline 18 -> CN 1
        region depth 37, baseline 18 -> CN 2
        region depth 74, baseline 18 -> CN 4
    """
    if baseline_depth <= 0:
        raise ValueError(
            "baseline_depth must be greater than zero"
        )

    bam = pysam.AlignmentFile(
        bam_path,
        "rb"
    )

    region_depths = {}
    copy_numbers = {}

    all_regions = {
        strandless_node(node)
        for cycle in cycles
        for node in cycle["nodes"][:-1]
    }

    logging.info(
        "Calculating depth for %d unique cycle regions",
        len(all_regions)
    )

    for region_number, region in enumerate(
        sorted(all_regions),
        start=1
    ):
        depth = calculate_region_depth(
            bam=bam,
            node=region,
            min_mapping_quality=min_mapping_quality
        )

        copy_number = calculate_copy_number(
            region_depth=depth,
            baseline_depth=baseline_depth
        )

        region_depths[region] = depth
        copy_numbers[region] = copy_number

        logging.info(
            "Region depth %s: depth=%.2f baseline=%.2f CN=%d",
            region,
            depth,
            baseline_depth,
            copy_number
        )

        if region_number % 100 == 0:
            logging.info(
                "Calculated depth for %d/%d regions",
                region_number,
                len(all_regions)
            )

    bam.close()

def calculate_copy_number(
    region_depth,
    baseline_depth
):
    """
    Estimate integer copy number from region depth.

    Examples:
        region depth 18, baseline 18 -> CN 1
        region depth 37, baseline 18 -> CN 2
        region depth 74, baseline 18 -> CN 4
    """
    if baseline_depth <= 0:
        raise ValueError(
            "baseline_depth must be greater than zero"
        )

    if region_depth <= 0:
        return 0

    return max(
        1,
        int(round(region_depth / baseline_depth))
    )


def calculate_region_copy_numbers(
    cycles,
    bam_path,
    baseline_depth,
    min_mapping_quality=0
):
    """
    Calculate depth and estimated copy number once for every unique
    genomic region found in the small cycles.
    """
    all_regions = {
        strandless_node(node)
        for cycle in cycles
        for node in cycle["nodes"][:-1]
    }

    region_depths = {}
    copy_numbers = {}

    logging.info(
        "Calculating depth for %d unique cycle regions",
        len(all_regions)
    )

    bam = pysam.AlignmentFile(
        bam_path,
        "rb"
    )

    try:
        for region_number, region in enumerate(
            sorted(all_regions),
            start=1
        ):
            depth = calculate_region_depth(
                bam=bam,
                node=region,
                min_mapping_quality=min_mapping_quality
            )

            copy_number = calculate_copy_number(
                region_depth=depth,
                baseline_depth=baseline_depth
            )

            region_depths[region] = depth
            copy_numbers[region] = copy_number

            logging.info(
                "Region %s: depth=%.2f baseline=%.2f CN=%d",
                region,
                depth,
                baseline_depth,
                copy_number
            )

            if region_number % 100 == 0:
                logging.info(
                    "Calculated region depth %d/%d",
                    region_number,
                    len(all_regions)
                )

    finally:
        bam.close()

    return region_depths, copy_numbers


def count_walk_regions(walk):
    """
    Count how many times each strandless genomic region occurs in an
    Euler walk.

    The repeated closing node is excluded.
    """
    if not walk:
        return Counter()

    if walk[0] == walk[-1]:
        effective_walk = walk[:-1]
    else:
        effective_walk = walk

    return Counter(
        strandless_node(node)
        for node in effective_walk
    )


def validate_walk_depth(
    walk,
    copy_numbers,
    copy_tolerance=0
):
    """
    Compare Euler-walk region occurrences with depth-estimated copy
    numbers.

    Returns:
        valid
        occurrences
        violations
    """
    occurrences = count_walk_regions(walk)
    violations = {}

    for region, observed_count in occurrences.items():
        allowed_count = copy_numbers.get(region, 0)

        if observed_count > allowed_count + copy_tolerance:
            violations[region] = {
                "observed": observed_count,
                "allowed": allowed_count
            }

    return (
        not violations,
        occurrences,
        violations
    )


def load_cycles(summary_file):
    """
    Read small cycles from summary.txt.
    """
    cycles = []

    with open(summary_file, "r") as handle:
        header = handle.readline()

        if not header:
            return cycles

        for line_number, line in enumerate(
            handle,
            start=2
        ):
            line = line.rstrip("\n")

            if not line:
                continue

            fields = line.split("\t")

            if len(fields) < 4:
                logging.warning(
                    "Skipping malformed line %d in %s",
                    line_number,
                    summary_file
                )
                continue

            try:
                cycle_size_value = int(fields[1])
                fragment_count = int(fields[2])
            except ValueError:
                logging.warning(
                    "Skipping line %d with invalid numeric fields",
                    line_number
                )
                continue

            nodes = parse_nodes(
                fields[3]
            )

            if len(nodes) < 2:
                continue

            if nodes[0] != nodes[-1]:
                nodes.append(nodes[0])

            cycles.append({
                "id": fields[0],
                "size": cycle_size_value,
                "frag_count": fragment_count,
                "path": fields[3],
                "nodes": nodes
            })

    return cycles


def build_cycle_groups(cycles):
    """
    Group small cycles when they share at least one genomic region.

    Strand is ignored only when deciding which cycles belong to the
    same event.
    """
    node_sets = [
        {
            strandless_node(node)
            for node in cycle["nodes"][:-1]
        }
        for cycle in cycles
    ]

    groups = []
    assigned = set()

    for index in range(len(cycles)):
        if index in assigned:
            continue

        component = {index}
        stack = [index]

        while stack:
            current = stack.pop()

            for other in range(len(cycles)):
                if other in component:
                    continue

                if node_sets[current].intersection(
                    node_sets[other]
                ):
                    component.add(other)
                    stack.append(other)

        assigned.update(component)
        groups.append(
            sorted(component)
        )

    return groups


def format_region_values(
    regions,
    values,
    numeric_format=None
):
    """
    Format a region-value dictionary for output.
    """
    result = []

    for region in sorted(regions):
        value = values.get(region, 0)

        if numeric_format is not None:
            value_text = numeric_format.format(value)
        else:
            value_text = str(value)

        result.append(
            "{}={}".format(
                region,
                value_text
            )
        )

    return ",".join(result)


def format_violations(violations):
    """
    Format depth-copy-number violations for output.
    """
    if not violations:
        return "none"

    return ",".join(
        "{}:observed={},allowed={}".format(
            region,
            values["observed"],
            values["allowed"]
        )
        for region, values in sorted(
            violations.items()
        )
    )


def merge_small_cycles(
    summary_file,
    outfile,
    bam_path,
    baseline_depth,
    min_mapping_quality=0,
    copy_tolerance=0,
    require_depth_support=False
):
    """
    Merge overlapping small cycles and construct strand-aware Euler
    walks.

    Region depth is calculated from the BAM file. The number of times
    each genomic region occurs in the Euler walk is compared with its
    depth-estimated copy number.

    If require_depth_support=False, all reconstructed cycles are
    written and marked as depth-consistent or depth-inconsistent.

    If require_depth_support=True, depth-inconsistent cycles are
    omitted.
    """
    if baseline_depth <= 0:
        raise ValueError(
            "baseline_depth must be greater than zero"
        )

    cycles = load_cycles(
        summary_file
    )

    header = (
        "bigcycle_id\t"
        "source_cycles\t"
        "total_size\t"
        "fragment_count\t"
        "depth_consistent\t"
        "region_depths\t"
        "region_copy_numbers\t"
        "region_occurrences\t"
        "depth_violations\t"
        "path\n"
    )

    if not cycles:
        with open(outfile, "w") as out:
            out.write(header)

        logging.warning(
            "No cycles found in %s",
            summary_file
        )

        return []

    region_depths, copy_numbers = (
        calculate_region_copy_numbers(
            cycles=cycles,
            bam_path=bam_path,
            baseline_depth=baseline_depth,
            min_mapping_quality=min_mapping_quality
        )
    )

    groups = build_cycle_groups(
        cycles
    )

    output_records = []
    bigcycle_id = 1

    for group in groups:
        source_ids = [
            cycles[index]["id"]
            for index in group
        ]

        graph = defaultdict(list)

        # Add every edge from every small cycle. Duplicate edges are
        # retained and can therefore cause repeated traversal.
        for index in group:
            nodes = cycles[index]["nodes"]

            for position in range(
                len(nodes) - 1
            ):
                left = nodes[position]
                right = nodes[position + 1]

                if left == right:
                    continue

                graph[left].append(right)
                graph[right].append(left)

        graph = remove_leaf_nodes(
            graph
        )

        components = connected_edge_components(
            graph
        )

        if not components:
            logging.warning(
                "Cycle group %s has no cyclic core",
                ",".join(source_ids)
            )
            continue

        for component_graph in components:
            odd_count = sum(
                1
                for neighbors in component_graph.values()
                if len(neighbors) % 2 == 1
            )

            if odd_count not in (0, 2):
                logging.warning(
                    "Skipping source cycles %s because the graph "
                    "has %d odd-degree nodes",
                    ",".join(source_ids),
                    odd_count
                )
                continue

            walk = euler_walk(
                component_graph
            )

            if len(walk) < 2:
                continue

            occurrences = count_walk_regions(
                walk
            )

            valid, occurrences, violations = (
                validate_walk_depth(
                    walk=walk,
                    copy_numbers=copy_numbers,
                    copy_tolerance=copy_tolerance
                )
            )

            if require_depth_support and not valid:
                logging.warning(
                    "Skipping depth-inconsistent cycle from %s: %s",
                    ",".join(source_ids),
                    violations
                )
                continue

            used_regions = sorted(
                occurrences
            )

            depth_text = format_region_values(
                used_regions,
                region_depths,
                numeric_format="{:.2f}"
            )

            copy_number_text = (
                format_region_values(
                    used_regions,
                    copy_numbers
                )
            )

            occurrence_text = (
                format_region_values(
                    used_regions,
                    occurrences
                )
            )

            violation_text = format_violations(
                violations
            )

            if walk[0] == walk[-1]:
                fragment_count = len(walk) - 1
            else:
                fragment_count = len(walk)

            total_size = path_size(
                walk
            )

            output_records.append({
                "bigcycle_id": bigcycle_id,
                "source_cycles": ",".join(source_ids),
                "total_size": total_size,
                "fragment_count": fragment_count,
                "depth_consistent": (
                    "yes" if valid else "no"
                ),
                "region_depths": depth_text,
                "region_copy_numbers": copy_number_text,
                "region_occurrences": occurrence_text,
                "depth_violations": violation_text,
                "path": " -> ".join(walk)
            })

            bigcycle_id += 1

    with open(outfile, "w") as out:
        out.write(header)

        for record in output_records:
            out.write(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                    record["bigcycle_id"],
                    record["source_cycles"],
                    record["total_size"],
                    record["fragment_count"],
                    record["depth_consistent"],
                    record["region_depths"],
                    record["region_copy_numbers"],
                    record["region_occurrences"],
                    record["depth_violations"],
                    record["path"]
                )
            )

    logging.info(
        "Wrote %d depth-checked big cycles -> %s",
        len(output_records),
        outfile
    )

    return output_records