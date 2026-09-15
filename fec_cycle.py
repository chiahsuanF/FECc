
#!/usr/bin/env python3

import logging
import math, pysam
from collections import Counter, defaultdict

import pysam


def parse_nodes(path):
    """
    Parse a semicolon-separated cycle path.
    """
    return [
        node.strip()
        for node in path.strip().split(";")
        if node.strip()
    ]


def strandless_node(node):
    """
    Remove the terminal strand symbol from a graph node.
    """
    if node.endswith(("+", "-")):
        return node[:-1]

    return node


def parse_region(node):
    """
    Parse a strand-aware genomic node.

    Returns:
        chromosome, start, end
    """
    region = strandless_node(node)

    chromosome, coordinates = region.rsplit(":", 1)
    start_text, end_text = coordinates.split("-", 1)

    return (
        chromosome,
        int(start_text),
        int(end_text)
    )


def path_size(path):
    """
    Calculate the total size of the fragments in a closed path.

    This version assumes zero-based, half-open intervals. If the
    summary file uses one-based inclusive intervals, change
    end - start to end - start + 1.
    """
    if not path:
        return 0

    effective_path = (
        path[:-1]
        if path[0] == path[-1]
        else path
    )

    total_size = 0

    for node in effective_path:
        try:
            _, start, end = parse_region(node)

            if end <= start:
                logging.warning(
                    "Invalid region coordinates: %s",
                    node
                )
                continue

            total_size += end - start

        except (IndexError, ValueError):
            logging.warning(
                "Could not parse region size from node: %s",
                node
            )

    return total_size


def calculate_region_depth(
    bam,
    node,
    min_mapping_quality=0
):
    """
    Calculate average read depth across a half-open genomic interval.

    Positions with zero coverage are included because the summed
    pileup depth is divided by the full interval length.
    """
    chromosome, start, end = parse_region(node)

    region_length = end - start

    if region_length <= 0:
        logging.warning(
            "Invalid region length for %s",
            node
        )
        return 0.0

    total_depth = 0

    try:
        for column in bam.pileup(
            chromosome,
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
    Estimate an integer copy number from depth.
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
    min_mapping_quality=30
):
    """
    Calculate depth and copy number once for each unique region.
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

    with pysam.AlignmentFile(bam_path, "rb") as bam:
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

    return region_depths, copy_numbers


def count_walk_regions(walk):
    """
    Count strandless fragment occurrences in a closed path.

    The final closing node is excluded.
    """
    if not walk:
        return Counter()

    effective_walk = (
        walk[:-1]
        if walk[0] == walk[-1]
        else walk
    )

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
    Compare path occurrences with depth-estimated copy numbers.
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

    Expected columns:
        cycle_id
        total_size
        fragment_count
        path

    Additional columns are ignored.
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

            nodes = parse_nodes(fields[3])

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

def flip_strand(node):
    """
    Flip the terminal strand orientation of a node.
    """
    if node.endswith("+"):
        return node[:-1] + "-"

    if node.endswith("-"):
        return node[:-1] + "+"

    return node


def canonical_cycle_key(nodes):
    """
    Canonicalize a closed strand-aware cycle.

    Cycles are considered equivalent when they differ only by:
      - circular rotation;
      - reverse-complement traversal.

    The repeated closing node is excluded from the key.
    """
    if not nodes:
        return tuple()

    body = list(nodes)

    if len(body) > 1 and body[0] == body[-1]:
        body = body[:-1]

    if not body:
        return tuple()

    forward_rotations = [
        tuple(
            body[index:]
            + body[:index]
        )
        for index in range(len(body))
    ]

    reverse_complement_body = [
        flip_strand(node)
        for node in reversed(body)
    ]

    reverse_complement_rotations = [
        tuple(
            reverse_complement_body[index:]
            + reverse_complement_body[:index]
        )
        for index in range(
            len(reverse_complement_body)
        )
    ]

    return min(
        forward_rotations
        + reverse_complement_rotations
    )




def find_duplicate_cycle_regions(cycle):
    """
    Return strandless genomic regions repeated internally.

    The final closing node is excluded because it intentionally
    repeats the first node.
    """
    if not cycle:
        return {}

    effective_cycle = (
        cycle[:-1]
        if cycle[0] == cycle[-1]
        else cycle
    )

    occurrences = Counter(
        strandless_node(node)
        for node in effective_cycle
    )

    return {
        region: count
        for region, count in occurrences.items()
        if count > 1
    }

def deduplicate_loaded_cycles(cycles):
    """
    Collapse complete cycles that differ only by circular rotation or
    reverse-complement traversal.

    The first cycle is retained as the representative. IDs of all
    equivalent cycles are preserved as duplicate evidence.
    """
    retained_cycles = {}
    duplicates_removed = 0

    for cycle in cycles:
        key = canonical_cycle_key(
            cycle["nodes"]
        )

        if key not in retained_cycles:
            retained_cycle = dict(cycle)

            retained_cycle["duplicate_cycle_ids"] = [
                cycle["id"]
            ]

            retained_cycle["duplicate_count"] = 1

            retained_cycles[key] = retained_cycle

        else:
            retained_cycles[key][
                "duplicate_cycle_ids"
            ].append(
                cycle["id"]
            )

            retained_cycles[key][
                "duplicate_count"
            ] += 1

            duplicates_removed += 1

    unique_cycles = list(
        retained_cycles.values()
    )

    logging.info(
        "Cycle deduplication: raw=%d unique=%d "
        "duplicates_removed=%d",
        len(cycles),
        len(unique_cycles),
        duplicates_removed
    )

    return unique_cycles


def canonical_undirected_edge(left, right):
    """
    Return a consistent strand-aware undirected edge key.
    """
    return tuple(sorted((left, right)))


def cycle_edge_set(cycle):
    """
    Return unique strand-aware undirected edges in one small cycle.
    """
    nodes = cycle["nodes"]

    return {
        canonical_undirected_edge(
            nodes[index],
            nodes[index + 1]
        )
        for index in range(len(nodes) - 1)
        if nodes[index] != nodes[index + 1]
    }


def build_cycle_groups(
    cycles,
    minimum_shared_regions=2
):
    """
    Group cycles transitively when they share enough strandless
    genomic regions.

    This preserves compatibility with the existing pipeline. Large
    groups are later resolved by simple-cycle extraction.
    """
    if minimum_shared_regions < 1:
        raise ValueError(
            "minimum_shared_regions must be at least 1"
        )

    region_sets = [
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

                shared_regions = (
                    region_sets[current]
                    & region_sets[other]
                )

                if (
                    len(shared_regions)
                    >= minimum_shared_regions
                ):
                    component.add(other)
                    stack.append(other)

        assigned.update(component)
        groups.append(sorted(component))

    return groups


def build_deduplicated_cycle_graph(
    group,
    cycles
):
    """
    Build a simple undirected strand-aware graph.

    Each structural connection is stored once. IDs of all small-cycle
    records supporting each connection are preserved as evidence.
    """
    unique_edges = set()
    edge_source_cycles = defaultdict(set)

    for cycle_index in group:
        cycle = cycles[cycle_index]
        nodes = cycle["nodes"]

        source_cycle_ids = cycle.get(
            "duplicate_cycle_ids",
            [cycle["id"]]
        )

        for position in range(len(nodes) - 1):
            left = nodes[position]
            right = nodes[position + 1]

            if left == right:
                continue

            edge = canonical_undirected_edge(
                left,
                right
            )

            unique_edges.add(
                edge
            )

            edge_source_cycles[edge].update(
                source_cycle_ids
            )

    graph = defaultdict(set)

    for left, right in unique_edges:
        graph[left].add(
            right
        )

        graph[right].add(
            left
        )

    logging.info(
        "Built deduplicated graph: nodes=%d edges=%d "
        "raw_edge_evidence=%d",
        len(graph),
        len(unique_edges),
        sum(
            len(source_ids)
            for source_ids
            in edge_source_cycles.values()
        )
    )

    return (
        graph,
        edge_source_cycles
    )

def graph_has_edges(graph):
    """
    Return True when a graph contains at least one edge.
    """
    return any(
        neighbors
        for neighbors in graph.values()
    )

def connected_edge_components(graph):
    """
    Split an undirected graph into connected components containing
    at least one edge.
    """
    components = []
    visited = set()

    for start in sorted(graph):
        if start in visited:
            continue

        if not graph.get(start, set()):
            visited.add(start)
            continue

        component_nodes = set()
        stack = [start]
        visited.add(start)

        while stack:
            node = stack.pop()
            component_nodes.add(node)

            for neighbor in graph.get(node, set()):
                if neighbor not in visited:
                    visited.add(neighbor)
                    stack.append(neighbor)

        component_graph = defaultdict(set)

        for node in component_nodes:
            component_graph[node] = {
                neighbor
                for neighbor in graph.get(node, set())
                if neighbor in component_nodes
            }

        if graph_has_edges(component_graph):
            components.append(component_graph)

    return components

def find_duplicate_cycle_edges(cycle):
    """
    Return undirected structural connections repeated within a cycle.
    """
    edge_counts = Counter()

    for index in range(len(cycle) - 1):
        edge = canonical_undirected_edge(
            cycle[index],
            cycle[index + 1]
        )

        edge_counts[edge] += 1

    return {
        edge: count
        for edge, count in edge_counts.items()
        if count > 1
    }

def canonical_simple_cycle(cycle):
    """
    Canonicalize a closed strand-aware candidate cycle.

    Circular rotations and reverse-complement traversal are treated
    as equivalent.
    """
    return canonical_cycle_key(
        cycle
    )

def is_valid_simple_cycle(cycle):
    """
    Validate a simple closed structural cycle.

    Requirements:
      - at least two distinct genomic regions;
      - first and final nodes are identical;
      - no internal strand-aware node repeats;
      - no internal strandless region repeats;
      - no structural connection repeats.

    Accepted:
        A-B-A
        A-B-C-A
        A-B-C-D-A

    Rejected:
        A-B-C-B-A
    """
    if len(cycle) < 3:
        return False

    if cycle[0] != cycle[-1]:
        return False

    body = cycle[:-1]

    if len(body) < 2:
        return False

    if len(body) != len(set(body)):
        return False

    strandless_regions = [
        strandless_node(node)
        for node in body
    ]

    if (
        len(strandless_regions)
        != len(set(strandless_regions))
    ):
        return False

    # A-B-A is a valid special two-region circle. The same undirected
    # edge appears twice because the path must return to A.
    if len(body) == 2:
        return True

    if find_duplicate_cycle_edges(cycle):
        return False

    return True


def find_two_region_cycles(
    group,
    cycles
):
    """
    Retain unique original A-B-A two-region cycles.

    The canonical representation is used only as a duplicate key.
    The original path and orientation are preserved.
    """
    detected = {}

    for cycle_index in group:
        cycle = cycles[cycle_index]
        nodes = cycle["nodes"]

        if len(nodes) != 3:
            continue

        if nodes[0] != nodes[-1]:
            continue

        strandless_regions = {
            strandless_node(node)
            for node in nodes[:-1]
        }

        if len(strandless_regions) != 2:
            continue

        if not is_valid_simple_cycle(nodes):
            continue

        key = canonical_simple_cycle(nodes)

        if key not in detected:
            detected[key] = list(nodes)

    return list(detected.values())


def find_simple_cycles(
    graph,
    maximum_cycle_length=30,
    maximum_candidate_cycles=100000
):
    """
    Enumerate unique simple cycles containing at least three nodes.

    No internal strand-aware or strandless region may repeat.
    """
    if maximum_cycle_length < 3:
        raise ValueError(
            "maximum_cycle_length must be at least 3"
        )

    if maximum_candidate_cycles < 1:
        raise ValueError(
            "maximum_candidate_cycles must be at least 1"
        )

    detected_cycles = {}
    search_stopped = False

    def search(
        start,
        current,
        path,
        used_nodes,
        used_regions
    ):
        nonlocal search_stopped

        if search_stopped:
            return

        if (
            len(detected_cycles)
            >= maximum_candidate_cycles
        ):
            search_stopped = True
            return

        if len(path) > maximum_cycle_length:
            return

        for neighbor in sorted(
            graph.get(current, set())
        ):
            if neighbor == start:
                if len(path) >= 3:
                    candidate = path + [start]

                    if not is_valid_simple_cycle(candidate):
                        continue

                    key = canonical_simple_cycle(candidate)

                    if key not in detected_cycles:
                        detected_cycles[key] = list(candidate)

                continue

            if neighbor in used_nodes:
                continue

            neighbor_region = strandless_node(
                neighbor
            )

            if neighbor_region in used_regions:
                continue

            search(
                start=start,
                current=neighbor,
                path=path + [neighbor],
                used_nodes=used_nodes | {neighbor},
                used_regions=(
                    used_regions
                    | {neighbor_region}
                )
            )

            if search_stopped:
                return

    for start in sorted(graph):
        search(
            start=start,
            current=start,
            path=[start],
            used_nodes={start},
            used_regions={
                strandless_node(start)
            }
        )

        if search_stopped:
            logging.warning(
                "Stopped simple-cycle search after %d candidates",
                maximum_candidate_cycles
            )
            break

    result = list(detected_cycles.values())

    logging.info(
        "Detected %d unique simple cycles",
        len(result)
    )

    return result


def supporting_cycles_for_candidate(
    candidate_cycle,
    edge_source_cycles
):
    """
    Return source-cycle IDs supporting at least one selected edge.
    """
    source_cycle_ids = set()

    for index in range(len(candidate_cycle) - 1):
        edge = canonical_undirected_edge(
            candidate_cycle[index],
            candidate_cycle[index + 1]
        )

        source_cycle_ids.update(
            edge_source_cycles.get(
                edge,
                set()
            )
        )

    return sorted(
        source_cycle_ids,
        key=lambda value: (
            0,
            int(value)
        )
        if str(value).isdigit()
        else (
            1,
            str(value)
        )
    )


def estimate_automatic_depth_range(
    region_depths,
    minimum_fold_gap=1.4,
    upper_fold_gap=2.0,
    maximum_upper_tail_fraction=0.10
):
    """
    Estimate one diagnostic depth range for a preliminary group.

    The first qualifying lower gap is used as the lower cutoff.
    An upper gap is used only when it isolates a small upper tail.
    """
    depths = sorted(
        float(depth)
        for depth in region_depths.values()
        if depth > 0
    )

    if not depths:
        return (
            0.0,
            0.0,
            {
                "reason": "no_positive_depths",
                "lower_gap": None,
                "upper_gap": None,
                "minimum_depth": 0.0,
                "maximum_depth": 0.0,
                "depth_count": 0
            }
        )

    if len(depths) < 3:
        minimum_depth = min(depths)
        maximum_depth = max(depths)

        return (
            minimum_depth,
            maximum_depth,
            {
                "reason": (
                    "not_enough_depths_for_gap_detection"
                ),
                "lower_gap": None,
                "upper_gap": None,
                "minimum_depth": minimum_depth,
                "maximum_depth": maximum_depth,
                "depth_count": len(depths)
            }
        )

    gaps = []

    for index in range(len(depths) - 1):
        lower_value = depths[index]
        upper_value = depths[index + 1]

        gaps.append({
            "index": index,
            "lower_value": lower_value,
            "upper_value": upper_value,
            "fold_gap": (
                upper_value / lower_value
            )
        })

    lower_candidates = [
        gap
        for gap in gaps
        if gap["fold_gap"] >= minimum_fold_gap
    ]

    if lower_candidates:
        lower_gap = min(
            lower_candidates,
            key=lambda gap: gap["index"]
        )

        minimum_depth = math.sqrt(
            lower_gap["lower_value"]
            * lower_gap["upper_value"]
        )
    else:
        lower_gap = None
        minimum_depth = min(depths)

    upper_candidates = []

    for gap in gaps:
        upper_tail_count = (
            len(depths)
            - gap["index"]
            - 1
        )

        upper_tail_fraction = (
            upper_tail_count
            / float(len(depths))
        )

        if (
            gap["fold_gap"] >= upper_fold_gap
            and upper_tail_fraction
            <= maximum_upper_tail_fraction
            and (
                lower_gap is None
                or gap["index"] > lower_gap["index"]
            )
        ):
            candidate = dict(gap)
            candidate["upper_tail_count"] = upper_tail_count
            candidate["upper_tail_fraction"] = (
                upper_tail_fraction
            )
            upper_candidates.append(candidate)

    if upper_candidates:
        upper_gap = max(
            upper_candidates,
            key=lambda gap: gap["fold_gap"]
        )

        maximum_depth = math.sqrt(
            upper_gap["lower_value"]
            * upper_gap["upper_value"]
        )
    else:
        upper_gap = None
        maximum_depth = max(depths)

    if minimum_depth > maximum_depth:
        minimum_depth = min(depths)
        maximum_depth = max(depths)
        lower_gap = None
        upper_gap = None
        reason = "invalid_range_use_observed_values"

    elif lower_gap is not None and upper_gap is not None:
        reason = "lower_gap_and_upper_outlier_gap"

    elif lower_gap is not None:
        reason = "lower_gap_and_observed_maximum"

    elif upper_gap is not None:
        reason = "observed_minimum_and_upper_outlier_gap"

    else:
        reason = "observed_minimum_and_maximum"

    return (
        minimum_depth,
        maximum_depth,
        {
            "reason": reason,
            "lower_gap": lower_gap,
            "upper_gap": upper_gap,
            "minimum_depth": minimum_depth,
            "maximum_depth": maximum_depth,
            "depth_count": len(depths)
        }
    )


def evaluate_group_depth_range(
    graph,
    region_depths,
    minimum_fold_gap=1.4,
    upper_fold_gap=2.0,
    maximum_upper_tail_fraction=0.10
):
    """
    Calculate one depth range for all unique regions in a group graph.

    The graph is not modified.
    """
    group_regions = {
        strandless_node(node)
        for node in graph
    }

    positive_depths = {
        region: region_depths[region]
        for region in group_regions
        if (
            region in region_depths
            and region_depths[region] > 0
        )
    }

    missing_regions = {
        region
        for region in group_regions
        if region not in positive_depths
    }

    if not positive_depths:
        removed_regions = {
            region: {
                "depth": 0.0,
                "reason": "missing_region_depth"
            }
            for region in missing_regions
        }

        return (
            0.0,
            0.0,
            {},
            removed_regions,
            {
                "reason": "no_positive_depths",
                "lower_gap": None,
                "upper_gap": None,
                "depth_count": 0
            }
        )

    if len(positive_depths) <= 3:
        minimum_depth = min(positive_depths.values())
        maximum_depth = max(positive_depths.values())

        kept_regions = dict(positive_depths)
        removed_regions = {}

        details = {
            "reason": "short_group_keep_positive_depths",
            "lower_gap": None,
            "upper_gap": None,
            "minimum_depth": minimum_depth,
            "maximum_depth": maximum_depth,
            "depth_count": len(positive_depths)
        }

    else:
        (
            minimum_depth,
            maximum_depth,
            details
        ) = estimate_automatic_depth_range(
            region_depths=positive_depths,
            minimum_fold_gap=minimum_fold_gap,
            upper_fold_gap=upper_fold_gap,
            maximum_upper_tail_fraction=(
                maximum_upper_tail_fraction
            )
        )

        kept_regions = {}
        removed_regions = {}

        for region, depth in positive_depths.items():
            if minimum_depth <= depth <= maximum_depth:
                kept_regions[region] = depth
            else:
                removed_regions[region] = {
                    "depth": depth,
                    "reason": (
                        "below_minimum_depth"
                        if depth < minimum_depth
                        else "above_maximum_depth"
                    )
                }

    for region in missing_regions:
        removed_regions[region] = {
            "depth": 0.0,
            "reason": "missing_region_depth"
        }

    return (
        minimum_depth,
        maximum_depth,
        kept_regions,
        removed_regions,
        details
    )

def score_simple_cycle(
    candidate_cycle,
    edge_source_cycles
):
    """
    Score one structurally valid simple cycle.

    Priority:
      1. More unique genomic regions.
      2. More edges supported by multiple source cycles.
      3. Stronger weakest-edge evidence.
      4. Greater total edge evidence.

    Depth does not affect candidate selection.
    """
    if not is_valid_simple_cycle(candidate_cycle):
        return None

    candidate_regions = {
        strandless_node(node)
        for node in candidate_cycle[:-1]
    }

    edge_evidence_counts = []

    for index in range(len(candidate_cycle) - 1):
        edge = canonical_undirected_edge(
            candidate_cycle[index],
            candidate_cycle[index + 1]
        )

        evidence_count = len(
            edge_source_cycles.get(
                edge,
                set()
            )
        )

        edge_evidence_counts.append(
            evidence_count
        )

    if not edge_evidence_counts:
        return None

    # A candidate extracted from this graph should never have a
    # missing edge. Reject it if the original-path preservation failed.
    if min(edge_evidence_counts) == 0:
        logging.warning(
            "Rejecting candidate with missing edge evidence: %s",
            " -> ".join(candidate_cycle)
        )
        return None

    multiple_cycle_edge_count = sum(
        count >= 2
        for count in edge_evidence_counts
    )

    minimum_edge_evidence = min(
        edge_evidence_counts
    )

    total_edge_evidence = sum(
        edge_evidence_counts
    )

    return {
        "cycle": candidate_cycle,
        "candidate_regions": candidate_regions,
        "unique_region_count": len(candidate_regions),
        "multiple_cycle_edge_count": (
            multiple_cycle_edge_count
        ),
        "minimum_edge_evidence": (
            minimum_edge_evidence
        ),
        "total_edge_evidence": (
            total_edge_evidence
        ),
        "edge_evidence_counts": (
            edge_evidence_counts
        ),
        "score": (
            len(candidate_regions),
            multiple_cycle_edge_count,
            minimum_edge_evidence,
            total_edge_evidence
        )
    }

def format_region_values(
    regions,
    values,
    numeric_format=None
):
    """
    Format region-value mappings for output.
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
    Format copy-number violations.
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


def format_depth_range_violations(
    removed_regions
):
    """
    Format candidate-level depth-range violations.
    """
    if not removed_regions:
        return "none"

    return ",".join(
        "{}:depth={:.2f},reason={}".format(
            region,
            values["depth"],
            values["reason"]
        )
        for region, values in sorted(
            removed_regions.items()
        )
    )

def merge_small_cycles(
    summary_file,
    outfile,
    bam_path,
    baseline_depth,
    min_mapping_quality=0,
    copy_tolerance=0,
    require_depth_support=False,
    minimum_shared_regions=2,
    maximum_cycle_length=30,
    maximum_candidate_cycles=100000,
    maximum_reported_candidates=1
):
    """
    Reconstruct simple big-cycle candidates from small cycles.

    Candidate selection uses structural evidence only.

    Depth and copy number are calculated for reporting. They do not
    affect candidate ranking unless require_depth_support is True.

    Workflow:
      1. Load and deduplicate source cycles.
      2. Calculate region depth and copy number.
      3. Group cycles by shared strandless regions.
      4. Build a deduplicated structural graph.
      5. Enumerate two-region and larger simple cycles.
      6. Preserve original graph paths during deduplication.
      7. Rank candidates by fragment count and edge evidence.
      8. Report the requested number of candidates.
    """
    if baseline_depth <= 0:
        raise ValueError(
            "baseline_depth must be greater than zero"
        )

    if minimum_shared_regions < 1:
        raise ValueError(
            "minimum_shared_regions must be at least 1"
        )

    if maximum_cycle_length < 3:
        raise ValueError(
            "maximum_cycle_length must be at least 3"
        )

    if maximum_candidate_cycles < 1:
        raise ValueError(
            "maximum_candidate_cycles must be at least 1"
        )

    if maximum_reported_candidates < 1:
        raise ValueError(
            "maximum_reported_candidates must be at least 1"
        )

    raw_cycles = load_cycles(
        summary_file
    )

    cycles = deduplicate_loaded_cycles(
        raw_cycles
    )

    output_fields = (
        "bigcycle_id",
        "group_cycles",
        "source_cycles",
        "candidate_rank",
        "candidate_count",
        "total_size",
        "fragment_count",
        "depth_range_consistent",
        "depth_consistent",
        "automatic_min_depth",
        "automatic_max_depth",
        "minimum_edge_evidence",
        "total_edge_evidence",
        "region_depths",
        "region_copy_numbers",
        "region_occurrences",
        "depth_range_violations",
        "depth_violations",
        "path"
    )

    header = (
        "\t".join(output_fields)
        + "\n"
    )

    if not cycles:
        with open(outfile, "w") as output_handle:
            output_handle.write(
                header
            )

        logging.warning(
            "No cycles found in %s",
            summary_file
        )

        return []

    logging.info(
        "Loaded summary cycles: raw=%d unique=%d",
        len(raw_cycles),
        len(cycles)
    )

    region_depths, copy_numbers = (
        calculate_region_copy_numbers(
            cycles=cycles,
            bam_path=bam_path,
            baseline_depth=baseline_depth,
            min_mapping_quality=min_mapping_quality
        )
    )

    groups = build_cycle_groups(
        cycles=cycles,
        minimum_shared_regions=minimum_shared_regions
    )

    logging.info(
        "Created %d preliminary cycle groups",
        len(groups)
    )

    output_records = []
    bigcycle_id = 1

    for group_number, group in enumerate(
        groups,
        start=1
    ):
        group_source_ids = sorted(
            {
                source_id
                for cycle_index in group
                for source_id in cycles[
                    cycle_index
                ].get(
                    "duplicate_cycle_ids",
                    [cycles[cycle_index]["id"]]
                )
            },
            key=lambda value: (
                0,
                int(value)
            )
            if str(value).isdigit()
            else (
                1,
                str(value)
            )
        )

        graph, edge_source_cycles = (
            build_deduplicated_cycle_graph(
                group=group,
                cycles=cycles
            )
        )

        if not graph_has_edges(graph):
            logging.warning(
                "Group %d has no graph edges",
                group_number
            )
            continue

        logging.info(
            "Group %d: cycles=%s nodes=%d edges=%d",
            group_number,
            ",".join(
                str(cycle_id)
                for cycle_id in group_source_ids
            ),
            len(graph),
            sum(
                len(neighbors)
                for neighbors in graph.values()
            ) // 2
        )

        group_candidates = []

        two_region_candidates = (
            find_two_region_cycles(
                group=group,
                cycles=cycles
            )
        )

        group_candidates.extend(
            two_region_candidates
        )

        components = connected_edge_components(
            graph
        )

        logging.info(
            "Group %d: two_region_candidates=%d "
            "connected_components=%d",
            group_number,
            len(two_region_candidates),
            len(components)
        )

        for component_number, component_graph in enumerate(
            components,
            start=1
        ):
            component_candidates = (
                find_simple_cycles(
                    graph=component_graph,
                    maximum_cycle_length=(
                        maximum_cycle_length
                    ),
                    maximum_candidate_cycles=(
                        maximum_candidate_cycles
                    )
                )
            )

            logging.info(
                "Group %d component %d: candidates=%d",
                group_number,
                component_number,
                len(component_candidates)
            )

            group_candidates.extend(
                component_candidates
            )

        logging.info(
            "Group %d: raw candidates=%d",
            group_number,
            len(group_candidates)
        )

        # Use canonical representations only as duplicate keys.
        # Preserve the original paths found in the graph.
        unique_candidates = {}

        for candidate in group_candidates:
            if not is_valid_simple_cycle(
                candidate
            ):
                continue

            key = canonical_simple_cycle(
                candidate
            )

            if key not in unique_candidates:
                unique_candidates[key] = list(
                    candidate
                )

        group_candidates = list(
            unique_candidates.values()
        )

        logging.info(
            "Group %d: valid unique candidates=%d",
            group_number,
            len(group_candidates)
        )

        if not group_candidates:
            logging.warning(
                "No valid simple cycles found for group %d",
                group_number
            )
            continue

        scored_candidates = []

        for candidate in group_candidates:
            scored = score_simple_cycle(
                candidate_cycle=candidate,
                edge_source_cycles=(
                    edge_source_cycles
                )
            )

            if scored is not None:
                scored_candidates.append(
                    scored
                )

        logging.info(
            "Group %d: scoreable candidates=%d",
            group_number,
            len(scored_candidates)
        )

        if not scored_candidates:
            logging.warning(
                "No scoreable candidates found for group %d",
                group_number
            )
            continue

        scored_candidates.sort(
            key=lambda candidate: (
                candidate["score"],
                canonical_simple_cycle(
                    candidate["cycle"]
                )
            ),
            reverse=True
        )

        candidate_count = len(
            scored_candidates
        )

        accepted_rank = 0

        for scored in scored_candidates:
            if accepted_rank >= maximum_reported_candidates:
                break

            selected_cycle = scored[
                "cycle"
            ]

            source_cycle_ids = (
                supporting_cycles_for_candidate(
                    candidate_cycle=selected_cycle,
                    edge_source_cycles=(
                        edge_source_cycles
                    )
                )
            )

            valid, occurrences, violations = (
                validate_walk_depth(
                    walk=selected_cycle,
                    copy_numbers=copy_numbers,
                    copy_tolerance=copy_tolerance
                )
            )

            if require_depth_support and not valid:
                logging.warning(
                    "Rejected group %d candidate because of "
                    "copy-number violations: %s",
                    group_number,
                    format_violations(
                        violations
                    )
                )
                continue

            accepted_rank += 1

            used_regions = sorted(
                occurrences
            )

            depth_values = [
                region_depths.get(
                    region,
                    0.0
                )
                for region in used_regions
                if region_depths.get(
                    region,
                    0.0
                ) > 0
            ]

            if depth_values:
                minimum_region_depth = min(
                    depth_values
                )

                maximum_region_depth = max(
                    depth_values
                )
            else:
                minimum_region_depth = 0.0
                maximum_region_depth = 0.0

            record = {
                "bigcycle_id": bigcycle_id,
                "group_cycles": ",".join(
                    str(cycle_id)
                    for cycle_id in group_source_ids
                ),
                "source_cycles": ",".join(
                    str(cycle_id)
                    for cycle_id in source_cycle_ids
                ),
                "candidate_rank": accepted_rank,
                "candidate_count": candidate_count,
                "total_size": path_size(
                    selected_cycle
                ),
                "fragment_count": (
                    len(selected_cycle) - 1
                ),
                "depth_range_consistent": (
                    "not_applied"
                ),
                "depth_consistent": (
                    "yes"
                    if valid
                    else "no"
                ),
                "automatic_min_depth": (
                    "{:.2f}".format(
                        minimum_region_depth
                    )
                ),
                "automatic_max_depth": (
                    "{:.2f}".format(
                        maximum_region_depth
                    )
                ),
                "minimum_edge_evidence": scored[
                    "minimum_edge_evidence"
                ],
                "total_edge_evidence": scored[
                    "total_edge_evidence"
                ],
                "region_depths": format_region_values(
                    used_regions,
                    region_depths,
                    numeric_format="{:.2f}"
                ),
                "region_copy_numbers": (
                    format_region_values(
                        used_regions,
                        copy_numbers
                    )
                ),
                "region_occurrences": (
                    format_region_values(
                        used_regions,
                        occurrences
                    )
                ),
                "depth_range_violations": (
                    "not_applied"
                ),
                "depth_violations": (
                    format_violations(
                        violations
                    )
                ),
                "path": " -> ".join(
                    selected_cycle
                )
            }

            output_records.append(
                record
            )

            logging.info(
                "Accepted big cycle %d: group=%d "
                "rank=%d/%d fragments=%d "
                "minimum_edge_evidence=%d "
                "total_edge_evidence=%d",
                bigcycle_id,
                group_number,
                accepted_rank,
                candidate_count,
                len(selected_cycle) - 1,
                scored["minimum_edge_evidence"],
                scored["total_edge_evidence"]
            )

            bigcycle_id += 1

    with open(outfile, "w") as output_handle:
        output_handle.write(
            header
        )

        for record in output_records:
            output_handle.write(
                "\t".join(
                    str(record[field])
                    for field in output_fields
                )
                + "\n"
            )

    logging.info(
        "Wrote %d simple big-cycle candidates -> %s",
        len(output_records),
        outfile
    )

    return output_records