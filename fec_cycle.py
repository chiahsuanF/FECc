#!/usr/bin/env python3

import glob
import logging
import math
import os
import pysam
import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from collections import Counter, defaultdict

def parse_nodes(path):
    return [
        node.strip()
        for node in path.strip().split(";")
        if node.strip()
    ]


def strandless_node(node):
    if node.endswith("+") or node.endswith("-"):
        return node[:-1]

    return node

def cycle_regions(cycle):
    return {
        strandless_node(node)
        for node in cycle["nodes"][:-1]
    }


def select_big_cycles(cycles):
    if not cycles:
        return []

    regions = {
        index: cycle_regions(cycle)
        for index, cycle in enumerate(cycles)
    }

    def rank(index):
        return (
            -len(regions[index]),
            -len(cycles[index]["nodes"]),
            str(cycles[index]["id"])
        )

    order = sorted(range(len(cycles)), key=rank)

    selected = []
    taken = []

    for index in order:
        covering = next(
            (
                chosen
                for chosen in taken
                if regions[index] <= regions[chosen]
            ),
            None
        )

        if covering is None:
            taken.append(index)
            selected.append((index, []))
        else:
            position = next(
                slot
                for slot, (chosen, _) in enumerate(selected)
                if chosen == covering
            )
            selected[position][1].append(index)

    return [
        (
            cycles[index],
            [cycles[other] for other in sorted(subsumed)]
        )
        for index, subsumed in selected
    ]

LEFT = "L"
RIGHT = "R"

def split_strand(node):
    if node.endswith("+") or node.endswith("-"):
        return node[:-1], node[-1]

    return node, "+"


def entry_end(strand):
    return LEFT if strand == "+" else RIGHT


def exit_end(strand):
    return RIGHT if strand == "+" else LEFT


def cycle_junctions(nodes):
    junctions = []

    for position in range(len(nodes) - 1):
        left_region, left_strand = split_strand(nodes[position])
        right_region, right_strand = split_strand(nodes[position + 1])

        junctions.append(
            tuple(sorted((
                (left_region, exit_end(left_strand)),
                (right_region, entry_end(right_strand))
            )))
        )

    return junctions


def build_region_groups(cycles):
    owners = defaultdict(set)

    for index, cycle in enumerate(cycles):
        for region in cycle_regions(cycle):
            owners[region].add(index)

    parent = list(range(len(cycles)))

    def find(index):
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    def union(left, right):
        left_root = find(left)
        right_root = find(right)

        if left_root != right_root:
            parent[max(left_root, right_root)] = min(left_root, right_root)

    for members in owners.values():
        ordered = sorted(members)

        for other in ordered[1:]:
            union(ordered[0], other)

    groups = defaultdict(list)

    for index in range(len(cycles)):
        groups[find(index)].append(index)

    return [members for _, members in sorted(groups.items())]

def estimate_depth_cutoff(depths, minimum_fold_gap=1.4):
    values = sorted(
        float(depth)
        for depth in depths
        if depth > 0
    )

    if len(values) < 3:
        return 0.0

    widest = None

    for position in range(len(values) - 1):
        if len(values) - position - 1 < 2:
            break

        below = values[position]
        above = values[position + 1]
        fold_gap = above / below

        if fold_gap < minimum_fold_gap:
            continue

        if widest is None or fold_gap > widest[0]:
            widest = (fold_gap, below, above)

    if widest is None:
        return 0.0

    return math.sqrt(widest[1] * widest[2])

def partition_group_by_depth(
    cycles,
    group,
    region_depths,
    minimum_fold_gap=1.4
):
    group_regions = set()

    for index in group:
        group_regions |= cycle_regions(cycles[index])

    cutoff = estimate_depth_cutoff(
        (
            region_depths.get(region, 0.0)
            for region in group_regions
        ),
        minimum_fold_gap
    )

    if cutoff <= 0:
        return list(group), [], 0.0, set()

    low_regions = {
        region
        for region in group_regions
        if region_depths.get(region, 0.0) < cutoff
    }

    core = []
    low = []

    for index in group:
        if cycle_regions(cycles[index]) & low_regions:
            low.append(index)
        else:
            core.append(index)

    if not core:
        # Every cycle reaches below the cutoff, so the gap separates
        # depth levels inside one structure rather than an amplicon from
        # its background.
        return list(group), [], 0.0, set()

    return core, low, cutoff, low_regions


def reuse_supported(multiplicities, copy_numbers, copy_tolerance=0):
    violations = {}

    for region, needed in sorted(multiplicities.items()):
        allowed = copy_numbers.get(region, 0)

        if needed > allowed + copy_tolerance:
            violations[region] = {
                "observed": needed,
                "allowed": allowed
            }

    return not violations, violations


def junction_multiplicities(cycles, group):
    counts = {}

    for index in group:
        for junction, count in Counter(cycles[index]["junctions"]).items():
            counts[junction] = max(counts.get(junction, 0), count)

    return counts


def segment_multiplicities(junction_counts):
    incidence = defaultdict(int)

    for (end_a, end_b), count in junction_counts.items():
        incidence[end_a] += count
        incidence[end_b] += count

    multiplicities = {}
    imbalances = []

    for region in sorted({region for region, _ in incidence}):
        left = incidence.get((region, LEFT), 0)
        right = incidence.get((region, RIGHT), 0)

        if left != right:
            imbalances.append((region, left, right))
        else:
            multiplicities[region] = left

    return multiplicities, imbalances


def opposite_end(end):
    return RIGHT if end == LEFT else LEFT


def build_junction_uses(junction_counts):
    uses = []
    slots = defaultdict(list)

    for (end_a, end_b), count in sorted(junction_counts.items()):
        for _ in range(count):
            use_id = len(uses)
            uses.append((end_a, end_b))
            slots[end_a].append(use_id)
            slots[end_b].append(use_id)

    return uses, slots


def other_endpoint(uses, use_id, end):
    end_a, end_b = uses[use_id]

    return end_b if end_a == end else end_a


def initial_transitions(slots, multiplicities):
    transitions = {}

    for region in sorted(multiplicities):
        left_slots = slots.get((region, LEFT), [])
        right_slots = slots.get((region, RIGHT), [])

        if len(left_slots) != len(right_slots):
            return None

        for left_use, right_use in zip(left_slots, right_slots):
            transitions[(region, left_use)] = right_use
            transitions[(region, right_use)] = left_use

    return transitions


def decompose_transition_walks(uses, transitions):
    walks = []
    consumed = set()

    for start_use in range(len(uses)):
        if start_use in consumed:
            continue

        walk = []
        use_id = start_use
        end = uses[start_use][0]

        while use_id not in consumed:
            consumed.add(use_id)

            region, arrival = end
            exit_use = transitions.get((region, use_id))

            if exit_use is None:
                return None

            walk.append((
                region,
                "+" if arrival == LEFT else "-",
                use_id if arrival == LEFT else exit_use,
                exit_use if arrival == LEFT else use_id
            ))

            end = other_endpoint(
                uses,
                exit_use,
                (region, opposite_end(arrival))
            )

            use_id = exit_use

        if walk:
            walks.append(walk)

    return walks


def merge_transition_walks(uses, transitions):
    while True:
        walks = decompose_transition_walks(uses, transitions)

        if walks is None:
            return None

        if len(walks) <= 1:
            return walks[0] if walks else None

        owners = defaultdict(list)

        for index, walk in enumerate(walks):
            for region, _, left_use, right_use in walk:
                owners[region].append((index, left_use, right_use))

        spliced = False

        for region in sorted(owners):
            entries = owners[region]
            first_index = entries[0][0]

            other = next(
                (entry for entry in entries if entry[0] != first_index),
                None
            )

            if other is None:
                continue

            _, left_a, right_a = entries[0]
            _, left_b, right_b = other

            transitions[(region, left_a)] = right_b
            transitions[(region, right_b)] = left_a
            transitions[(region, left_b)] = right_a
            transitions[(region, right_a)] = left_b

            spliced = True
            break

        if not spliced:
            # No segment is shared, so the walks are genuinely separate
            # molecules rather than one.
            return None


def reconstruct_big_cycle(junction_counts, multiplicities):
    if not junction_counts or not multiplicities:
        return []

    uses, slots = build_junction_uses(junction_counts)

    transitions = initial_transitions(slots, multiplicities)

    if transitions is None:
        return []

    walk = merge_transition_walks(uses, transitions)

    if not walk:
        return []

    if len(walk) != len(uses):
        return []

    path = [region + strand for region, strand, _, _ in walk]

    return path + [path[0]]

def solve_appearance_counts(
    junctions,
    end_support,
    region_depths,
    min_support=3,
    max_uses=6,
    target_cap_factor=2.0
):
    import numpy as np
    from scipy.optimize import Bounds, LinearConstraint, milp

    has_support = bool(end_support)

    candidates = sorted(
        junction
        for junction in set(junctions)
        if not has_support or end_support.get(junction, 0) >= min_support
    )

    if not candidates:
        return None, {}

    regions = sorted({
        region
        for junction in candidates
        for region, _ in junction
    })

    depths = {
        region: region_depths.get(region, 0.0)
        for region in regions
    }

    positive = [depth for depth in depths.values() if depth > 0]
    unit = min(positive) if positive else 1.0

    raw_targets = {
        region: (depth / unit if depth > 0 else 1.0)
        for region, depth in depths.items()
    }

    cap = target_cap_factor * float(np.median(list(raw_targets.values())))

    targets = {
        region: min(target, cap)
        for region, target in raw_targets.items()
    }

    junction_count = len(candidates)
    region_count = len(regions)
    row_of = {region: index for index, region in enumerate(regions)}

    # Variables: x_j, then the scale v, then e_j per junction, then d_r
    # per region.
    scale_column = junction_count
    junction_error = junction_count + 1
    region_error = junction_error + junction_count
    variable_count = region_error + region_count

    balance = np.zeros((region_count, variable_count))
    appearance = np.zeros((region_count, variable_count))

    for column, junction in enumerate(candidates):
        for region, end in junction:
            row = row_of[region]

            if end == LEFT:
                balance[row, column] += 1
                appearance[row, column] += 1
            else:
                balance[row, column] -= 1

    rows = []
    upper = []

    if has_support:
        # |x_j - v * reads_j| <= e_j
        for column, junction in enumerate(candidates):
            reads = end_support[junction]

            for sign in (1, -1):
                row = np.zeros(variable_count)
                row[column] = sign
                row[scale_column] = -sign * reads
                row[junction_error + column] = -1
                rows.append(row)
                upper.append(0.0)

    # |t_r - target_r| <= d_r
    for region, row_index in row_of.items():
        over = appearance[row_index].copy()
        over[region_error + row_index] = -1
        rows.append(over)
        upper.append(targets[region])

        under = -appearance[row_index]
        under[region_error + row_index] = -1
        rows.append(under)
        upper.append(-targets[region])

    used = np.zeros(variable_count)
    used[:junction_count] = 1

    cost = np.zeros(variable_count)
    cost[junction_error:region_error] = 1
    cost[region_error:] = 1

    lower_bounds = np.zeros(variable_count)
    upper_bounds = np.full(variable_count, np.inf)
    upper_bounds[:junction_count] = max_uses

    if not has_support:
        # Without reads to weigh them, every proposed junction is taken
        # as real.
        lower_bounds[:junction_count] = 1

    integrality = np.zeros(variable_count)
    integrality[:junction_count] = 1

    result = milp(
        cost,
        constraints=[
            LinearConstraint(balance, 0, 0),
            LinearConstraint(np.array(rows), -np.inf, upper),
            LinearConstraint(used, 1, np.inf),
        ],
        integrality=integrality,
        bounds=Bounds(lower_bounds, upper_bounds)
    )

    if not result.success:
        return None, targets

    counts = {
        junction: int(round(value))
        for junction, value in zip(candidates, result.x[:junction_count])
    }

    return (
        {
            junction: count
            for junction, count in counts.items()
            if count > 0
        },
        targets
    )


def path_size(path):
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
    min_mapping_quality=30
):
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
                "nodes": nodes,
                "junctions": cycle_junctions(nodes)
            })

    return cycles


def format_region_values(
    regions,
    values,
    numeric_format=None
):
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

def build_region_index(regions):
    per_chromosome = defaultdict(list)

    for region in regions:
        chrom, start, end = parse_region(region)
        per_chromosome[chrom].append((start, end, region))

    for intervals in per_chromosome.values():
        intervals.sort()

    return per_chromosome


def find_containing_region(region_index, chrom, position):
    intervals = region_index.get(chrom)

    if not intervals:
        return None

    low = 0
    high = len(intervals) - 1

    while low <= high:
        middle = (low + high) // 2

        start, end, region = intervals[middle]

        if start <= position <= end:
            return region

        if position < start:
            high = middle - 1
        else:
            low = middle + 1

    return None


def count_junction_support(
    breakpoints,
    regions,
    min_mapping_quality=30
):
    region_index = build_region_index(regions)

    reads = defaultdict(set)

    for breakpoint in breakpoints:
        (
            qname,
            chrom_1, position_1, quality_1, _strand_1,
            chrom_2, position_2, quality_2, _strand_2
        ) = breakpoint

        if quality_1 < min_mapping_quality:
            continue

        if quality_2 < min_mapping_quality:
            continue

        region_1 = find_containing_region(
            region_index,
            chrom_1,
            position_1
        )

        if region_1 is None:
            continue

        region_2 = find_containing_region(
            region_index,
            chrom_2,
            position_2
        )

        if region_2 is None:
            continue

        if region_1 == region_2:
            continue

        reads[tuple(sorted((region_1, region_2)))].add(qname)

    return {
        pair: len(qnames)
        for pair, qnames in reads.items()
    }


def nearest_end(region, position):
    _, start, end = parse_region(region)

    return LEFT if position - start <= end - position else RIGHT


def collect_end_junction_reads(
    breakpoints,
    regions,
    min_mapping_quality=30
):
    region_index = build_region_index(regions)

    reads = defaultdict(set)

    for breakpoint in breakpoints:
        (
            qname,
            chrom_1, position_1, quality_1, _strand_1,
            chrom_2, position_2, quality_2, _strand_2
        ) = breakpoint

        if quality_1 < min_mapping_quality:
            continue

        if quality_2 < min_mapping_quality:
            continue

        region_1 = find_containing_region(
            region_index,
            chrom_1,
            position_1
        )

        if region_1 is None:
            continue

        region_2 = find_containing_region(
            region_index,
            chrom_2,
            position_2
        )

        if region_2 is None:
            continue

        if region_1 == region_2:
            continue

        junction = tuple(sorted((
            (region_1, nearest_end(region_1, position_1)),
            (region_2, nearest_end(region_2, position_2))
        )))

        reads[junction].add(qname)

    return reads


def count_end_junction_support(
    breakpoints,
    regions,
    min_mapping_quality=30
):
    return {
        junction: len(qnames)
        for junction, qnames in collect_end_junction_reads(
            breakpoints,
            regions,
            min_mapping_quality
        ).items()
    }

SOFT_CLIP = 4
HARD_CLIP = 5


def alignment_breakpoint(alignment):
    cigar = alignment.cigartuples

    if not cigar:
        return None

    def leading_clip(operations):
        total = 0

        for operation, length in operations:
            if operation not in (SOFT_CLIP, HARD_CLIP):
                break
            total += length

        return total

    left = leading_clip(cigar)
    right = leading_clip(reversed(cigar))

    start = alignment.reference_start + 1

    if left > right:
        return start

    if right > left:
        return start + alignment.reference_length

    return None


def exact_junction_breakpoints(
    bam_path,
    junctions,
    junction_reads,
    min_mapping_quality=30
):
    wanted = defaultdict(set)

    for junction in junctions:
        for region, _ in junction:
            wanted[region] |= junction_reads.get(junction, set())

    positions = defaultdict(Counter)

    bam = pysam.AlignmentFile(bam_path, "rb")

    try:
        for region, qnames in sorted(wanted.items()):
            if not qnames:
                continue

            chrom, start, end = parse_region(region)

            per_read = defaultdict(list)

            try:
                alignments = bam.fetch(chrom, max(0, start - 1), end)

                for alignment in alignments:
                    if alignment.query_name not in qnames:
                        continue

                    if (
                        alignment.is_unmapped
                        or alignment.is_secondary
                        or alignment.is_qcfail
                        or alignment.mapping_quality < min_mapping_quality
                    ):
                        continue

                    position = alignment_breakpoint(alignment)

                    if position is not None and start <= position <= end + 1:
                        per_read[alignment.query_name].append(position)

            except ValueError:
                logging.warning(
                    "Could not fetch %s for exact breakpoints",
                    region
                )
                continue

            for junction in junctions:
                if region not in (junction[0][0], junction[1][0]):
                    continue

                for qname in junction_reads.get(junction, set()):
                    for position in per_read.get(qname, []):
                        positions[(junction, region)][position] += 1

    finally:
        bam.close()

    return {
        key: "{}:{}".format(
            parse_region(key[1])[0],
            min(
                counter,
                key=lambda position: (-counter[position], position)
            )
        )
        for key, counter in positions.items()
        if counter
    }


def format_path_breakpoints(walk, breakpoints):
    a_side = []
    b_side = []

    for position in range(len(walk) - 1):
        junction = cycle_junctions([walk[position], walk[position + 1]])[0]

        a_side.append(breakpoints.get(
            (junction, strandless_node(walk[position])), "NA"
        ))
        b_side.append(breakpoints.get(
            (junction, strandless_node(walk[position + 1])), "NA"
        ))

    return ",".join(a_side), ",".join(b_side)


def load_breakpoints(path):
    breakpoints = []

    with open(path, "r") as handle:
        header = handle.readline()

        if not header:
            return breakpoints

        for line_number, line in enumerate(handle, start=2):
            fields = line.rstrip("\n").split("\t")

            if len(fields) < 9:
                continue

            try:
                breakpoints.append((
                    fields[0],
                    fields[1],
                    int(fields[2]),
                    int(fields[3]),
                    fields[4],
                    fields[5],
                    int(fields[6]),
                    int(fields[7]),
                    fields[8]
                ))
            except ValueError:
                logging.warning(
                    "Skipping breakpoint line %d with invalid numeric fields",
                    line_number
                )

    return breakpoints


def find_breakpoint_table(outfile):
    directory = os.path.dirname(os.path.abspath(outfile))

    candidates = sorted(
        glob.glob(os.path.join(directory, "03_*_breakpoints.txt"))
    )

    if not candidates:
        return None

    if len(candidates) > 1:
        logging.warning(
            "Several breakpoint tables in %s; using %s",
            directory,
            candidates[0]
        )

    return candidates[0]


def calculate_path_junction_supports(
    walk,
    end_support
):
    supports = []

    for position in range(len(walk) - 1):
        left = walk[position]
        right = walk[position + 1]

        junction = cycle_junctions([left, right])[0]

        supports.append(
            "{}->{}={}".format(
                left,
                right,
                end_support.get(junction, 0)
            )
        )

    return ",".join(supports)

def find_low_depth_variants(
    base_walk,
    low_regions,
    end_support,
    min_support=3,
    max_uses=6
):
    base_regions = {strandless_node(node) for node in base_walk}

    strong = {
        junction: reads
        for junction, reads in end_support.items()
        if reads >= min_support
    }

    neighbours = defaultdict(set)

    for (end_a, end_b) in strong:
        if end_a[0] != end_b[0]:
            neighbours[end_a[0]].add(end_b[0])
            neighbours[end_b[0]].add(end_a[0])

    remaining = {
        region
        for region in low_regions
        if region not in base_regions
    }

    low_sets = []

    while remaining:
        seed = min(remaining)
        component = {seed}
        frontier = [seed]

        while frontier:
            region = frontier.pop()

            for other in neighbours[region]:
                if other in remaining and other not in component:
                    component.add(other)
                    frontier.append(other)

        remaining -= component

        if any(neighbours[region] & base_regions for region in component):
            low_sets.append(component)

    base_counts = Counter(cycle_junctions(base_walk))
    variants = []

    for low_set in low_sets:
        allowed = base_regions | low_set

        candidates = sorted(
            junction
            for junction in strong
            if all(region in allowed for region, _ in junction)
        )

        regions = sorted({
            region
            for junction in candidates
            for region, _ in junction
        })

        if not candidates or not low_set <= set(regions):
            logging.info(
                "No supported junctions join low-depth region(s) %s",
                ",".join(sorted(low_set))
            )
            continue

        row_of = {region: index for index, region in enumerate(regions)}

        balance = np.zeros((len(regions), len(candidates)))
        appearance = np.zeros((len(regions), len(candidates)))

        for column, junction in enumerate(candidates):
            for region, end in junction:
                if end == LEFT:
                    balance[row_of[region], column] += 1
                    appearance[row_of[region], column] += 1
                else:
                    balance[row_of[region], column] -= 1

        low_rows = [row_of[region] for region in sorted(low_set)]

        most = max(strong[junction] for junction in candidates)

        cost = np.array([
            1.0 + 1e-3 * (1.0 - strong[junction] / most)
            for junction in candidates
        ])

        result = milp(
            cost,
            constraints=[
                LinearConstraint(balance, 0, 0),
                LinearConstraint(appearance[low_rows], 1, 1),
            ],
            integrality=np.ones(len(candidates)),
            bounds=Bounds(0, max_uses)
        )

        if not result.success:
            logging.info(
                "No closed add-on through low-depth region(s) %s",
                ",".join(sorted(low_set))
            )
            continue

        added = {
            junction: int(round(value))
            for junction, value in zip(candidates, result.x)
            if round(value) > 0
        }

        junction_counts = Counter(base_counts)
        junction_counts.update(added)

        multiplicities, imbalances = segment_multiplicities(junction_counts)

        walk = (
            []
            if imbalances
            else reconstruct_big_cycle(dict(junction_counts), multiplicities)
        )

        if len(walk) < 2:
            logging.info(
                "Add-on through %s does not join the big cycle into one "
                "molecule",
                ",".join(sorted(low_set))
            )
            continue

        variants.append((walk, low_set, added))

    return variants


def merge_small_cycles(
    summary_file,
    outfile,
    bam_path,
    baseline_depth,
    min_mapping_quality=0,
    copy_tolerance=0,
    require_depth_support=False,
    apply_depth_filter=True,
    depth_filter_fold_gap=1.4,
    breakpoints=None,
    junction_mapping_quality=30,
    min_junction_support=3,
    low_depth_variants=True
):
    if baseline_depth <= 0:
        raise ValueError(
            "baseline_depth must be greater than zero"
        )

    cycles = load_cycles(
        summary_file
    )

    header = (
        "bigcycle_id\t"
        "total_size\t"
        "fragment_count\t"
        "region_depths\t"
        "region_copy_numbers\t"
        "junction_supports\t"
        "path\t"
        "A_breakpoint\t"
        "B_breakpoint\n"
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

    cycle_region_set = {
        strandless_node(node)
        for cycle in cycles
        for node in cycle["nodes"][:-1]
    }

    if not breakpoints:
        table = find_breakpoint_table(outfile)

        if table is None:
            logging.warning(
                "No breakpoints passed in and no 03_*_breakpoints.txt "
                "beside %s; every junction support will report as 0",
                outfile
            )
            breakpoints = []
        else:
            logging.info(
                "No breakpoints passed in; reading them from %s",
                table
            )
            breakpoints = load_breakpoints(table)

    junction_reads = collect_end_junction_reads(
        breakpoints=breakpoints,
        regions=cycle_region_set,
        min_mapping_quality=junction_mapping_quality
    )

    end_support = {
        junction: len(qnames)
        for junction, qnames in junction_reads.items()
    }

    logging.info(
        "Counted split-read support for %d segment-end junction(s) "
        "among %d cycle region(s) from %d breakpoint(s) at MQ >= %d",
        len(end_support),
        len(cycle_region_set),
        len(breakpoints),
        junction_mapping_quality
    )

    if breakpoints and not end_support:
        logging.warning(
            "No breakpoint joins two cycle regions: the breakpoint table "
            "and the cycle coordinates are probably from different runs"
        )

    output_records = []
    held_back = []

    def add_record(walk, source_ids, cutoff=0.0):
        valid, occurrences, violations = validate_walk_depth(
            walk=walk,
            copy_numbers=copy_numbers,
            copy_tolerance=copy_tolerance
        )

        if require_depth_support and not valid:
            logging.warning(
                "Skipping depth-inconsistent cycle from %s: %s",
                ",".join(source_ids),
                violations
            )
            return

        used_regions = sorted(occurrences)

        junction_supports = calculate_path_junction_supports(
            walk=walk,
            end_support=end_support
        )

        output_records.append({
            "bigcycle_id": len(output_records) + 1,
            "total_size": path_size(walk),
            "fragment_count": (
                len(walk) - 1
                if walk[0] == walk[-1]
                else len(walk)
            ),
            "region_depths": format_region_values(
                used_regions,
                region_depths,
                numeric_format="{:.2f}"
            ),
            "region_copy_numbers": format_region_values(
                used_regions,
                copy_numbers
            ),
            "junction_supports": junction_supports,
            "path": " -> ".join(walk)
        })

    def ids_of(group):
        return [cycles[index]["id"] for index in group]

    def supported_cycles(group):
        kept = []

        for index in group:
            cycle = cycles[index]

            weak = [
                junction
                for junction in cycle["junctions"]
                if end_support.get(junction, 0) < min_junction_support
            ]

            if end_support and weak:
                logging.info(
                    "Dropping cycle %s: %d junction(s) below %d split "
                    "reads: %s",
                    cycle["id"],
                    len(weak),
                    min_junction_support,
                    "; ".join(
                        "{}{}~{}{}={}".format(
                            end_a[0], end_a[1], end_b[0], end_b[1],
                            end_support.get((end_a, end_b), 0)
                        )
                        for end_a, end_b in weak
                    )
                )
                continue

            kept.append(cycle)

        return kept

    def add_separately(group, cutoff, reason):
        chosen = select_big_cycles(supported_cycles(group))

        logging.info(
            "Reporting cycle(s) %s separately: %s",
            ",".join(cycle["id"] for cycle, _ in chosen),
            reason
        )

        for cycle, subsumed in chosen:
            add_record(
                cycle["nodes"],
                [cycle["id"]] + [other["id"] for other in subsumed],
                cutoff=cutoff
            )

    def merge_group(group, cutoff):
        if len(group) == 1:
            for cycle in supported_cycles(group):
                add_record(
                    cycle["nodes"],
                    ids_of(group),
                    cutoff=cutoff
                )
            return

        proposed = {
            junction
            for index in group
            for junction in cycles[index]["junctions"]
        }

        group_regions = set()
        for index in group:
            group_regions |= cycle_regions(cycles[index])

        proposed |= {
            junction
            for junction, reads in end_support.items()
            if reads >= min_junction_support
            and all(region in group_regions for region, _ in junction)
        }

        junction_counts, targets = solve_appearance_counts(
            junctions=proposed,
            end_support=end_support,
            region_depths=region_depths,
            min_support=min_junction_support
        )

        if not junction_counts:
            add_separately(
                group,
                cutoff,
                "no balanced use of the supported junctions exists"
            )
            return

        multiplicities, imbalances = segment_multiplicities(
            junction_counts
        )

        walk = (
            []
            if imbalances
            else reconstruct_big_cycle(junction_counts, multiplicities)
        )

        if len(walk) < 2:
            add_separately(
                group,
                cutoff,
                "the chosen junctions do not form a single closed molecule"
            )
            return

        unused = sorted(proposed - set(junction_counts))

        if unused:
            logging.info(
                "Not using %d proposed junction(s): %s",
                len(unused),
                "; ".join(
                    "{}{}~{}{}={}".format(
                        end_a[0], end_a[1], end_b[0], end_b[1],
                        end_support.get((end_a, end_b), 0)
                    )
                    for end_a, end_b in unused
                )
            )

        logging.info(
            "Connected cycle(s) %s into a %d-fragment big cycle; "
            "appearances (depth target): %s",
            ",".join(ids_of(group)),
            len(walk) - 1,
            ",".join(
                "{}={}({:.2f})".format(
                    region,
                    multiplicities[region],
                    targets.get(region, 0.0)
                )
                for region in sorted(multiplicities)
            )
        )

        add_record(walk, ids_of(group), cutoff=cutoff)

    for group in build_region_groups(cycles):

        if apply_depth_filter:
            core, low, cutoff, low_regions = partition_group_by_depth(
                cycles=cycles,
                group=group,
                region_depths=region_depths,
                minimum_fold_gap=depth_filter_fold_gap
            )
        else:
            core, low, cutoff, low_regions = list(group), [], 0.0, set()

        if low:
            logging.info(
                "Depth cutoff %.2fx splits a %d-cycle group into %d core "
                "cycle(s) and %d low-depth cycle(s); below the cutoff: %s",
                cutoff,
                len(group),
                len(core),
                len(low),
                format_region_values(
                    sorted(low_regions),
                    region_depths,
                    numeric_format="{:.2f}"
                )
            )

        for cycle, subsumed in select_big_cycles(
            supported_cycles(low)
        ):
            add_record(
                cycle["nodes"],
                [cycle["id"]] + [other["id"] for other in subsumed],
                cutoff=cutoff
            )

        core_cycles = [cycles[index] for index in core]

        core_start = len(output_records)

        for members in build_region_groups(core_cycles):
            merge_group(
                [core[position] for position in members],
                cutoff
            )

        if low_regions:
            held_back.append((
                low_regions,
                range(core_start, len(output_records)),
                cutoff
            ))

    original_count = len(output_records)
    variant_notes = {}

    if low_depth_variants and end_support:
        for low_regions, record_indices, cutoff in held_back:
            for index in record_indices:
                base = output_records[index]

                for walk, low_set, added in find_low_depth_variants(
                    base_walk=base["path"].split(" -> "),
                    low_regions=low_regions,
                    end_support=end_support,
                    min_support=min_junction_support
                ):
                    before = len(output_records)

                    add_record(
                        walk,
                        ["variant_of_{}".format(base["bigcycle_id"])],
                        cutoff=cutoff
                    )

                    if len(output_records) == before:
                        continue

                    added_text = ",".join(
                        "{}{}~{}{}={}{}".format(
                            end_a[0], end_a[1], end_b[0], end_b[1],
                            end_support.get((end_a, end_b), 0),
                            "" if count == 1 else "x{}".format(count)
                        )
                        for (end_a, end_b), count in sorted(added.items())
                    )

                    variant_notes[output_records[-1]["bigcycle_id"]] = (
                        base["bigcycle_id"],
                        ",".join(sorted(low_set)),
                        added_text
                    )

                    logging.info(
                        "Added big cycle %d: variant of big cycle %d "
                        "carrying low-depth region(s) %s via %s",
                        output_records[-1]["bigcycle_id"],
                        base["bigcycle_id"],
                        format_region_values(
                            sorted(low_set),
                            region_depths,
                            numeric_format="{:.2f}"
                        ),
                        added_text
                    )

    path_junctions = {
        junction
        for record in output_records
        for junction in cycle_junctions(record["path"].split(" -> "))
    }

    exact = (
        exact_junction_breakpoints(
            bam_path=bam_path,
            junctions=path_junctions,
            junction_reads=junction_reads,
            min_mapping_quality=junction_mapping_quality
        )
        if path_junctions and junction_reads
        else {}
    )

    for record in output_records:
        record["A_breakpoint"], record["B_breakpoint"] = (
            format_path_breakpoints(record["path"].split(" -> "), exact)
        )

    with open(outfile, "w") as out:
        out.write(header)

        for record in output_records:
            out.write(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                    record["bigcycle_id"],
                    record["total_size"],
                    record["fragment_count"],
                    record["region_depths"],
                    record["region_copy_numbers"],
                    record["junction_supports"],
                    record["path"],
                    record["A_breakpoint"],
                    record["B_breakpoint"]
                )
            )


    if low_depth_variants:
        variant_file = os.path.join(
            os.path.dirname(os.path.abspath(outfile)),
            "summary_bigcycle_variants.txt"
        )

        with open(variant_file, "w") as out:
            out.write(
                header.rstrip("\n")
                + "\tbase_bigcycle_id\tlow_regions\tadded_junctions\n"
            )

            for record in output_records[original_count:]:
                base_id, low_text, added_text = variant_notes[
                    record["bigcycle_id"]
                ]

                out.write(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                        record["bigcycle_id"],
                        record["total_size"],
                        record["fragment_count"],
                        record["region_depths"],
                        record["region_copy_numbers"],
                        record["junction_supports"],
                        record["path"],
                        record["A_breakpoint"],
                        record["B_breakpoint"],
                        base_id,
                        low_text,
                        added_text
                    )
                )

        logging.info(
            "Wrote %d low-depth variant big cycle(s) -> %s",
            len(output_records) - original_count,
            variant_file
        )

    logging.info(
        "Wrote %d depth-checked big cycles -> %s",
        len(output_records),
        outfile
    )

    return output_records
