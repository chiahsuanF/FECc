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


# ---------------------------------------------------------------
# Segment-end breakpoint graph
# ---------------------------------------------------------------
#
# A segment is a strandless genomic region with two ends, L and R, in
# reference orientation. Traversing it as "+" enters L and exits R;
# traversing as "-" enters R and exits L.
#
# A junction is the unordered pair of segment ends that a split read
# joins. Because the pair is unordered, the same physical junction
# canonicalizes identically regardless of which direction a cycle
# happens to traverse it: A+ -> B+ and B- -> A- both yield
# {(A, R), (B, L)}. That is what lets one molecule contain a segment in
# both orientations, which is what an inversion is.

LEFT = "L"
RIGHT = "R"


def split_strand(node):
    """
    Split a strand-aware node into its region and strand.

    Example:
        chr2:100-200+ -> ("chr2:100-200", "+")
    """
    if node.endswith("+") or node.endswith("-"):
        return node[:-1], node[-1]

    return node, "+"


def entry_end(strand):
    """Return the end a segment is entered through."""
    return LEFT if strand == "+" else RIGHT


def exit_end(strand):
    """Return the end a segment is left through."""
    return RIGHT if strand == "+" else LEFT


def cycle_junctions(nodes):
    """
    Convert a closed cycle path into its canonical junctions.

    Each junction is a sorted pair of (region, end) tuples.
    """
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


def build_junction_groups(cycles):
    """
    Group small cycles that share at least one junction.

    A shared junction is direct read evidence that two cycles describe
    the same molecule. A merely shared region is not: two independent
    ecDNA species can both contain the same amplified region without
    being the same molecule.
    """
    junction_owners = defaultdict(set)

    for index, cycle in enumerate(cycles):
        for junction in cycle["junctions"]:
            junction_owners[junction].add(index)

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

    for owners in junction_owners.values():
        ordered = sorted(owners)

        for other in ordered[1:]:
            union(ordered[0], other)

    groups = defaultdict(list)

    for index in range(len(cycles)):
        groups[find(index)].append(index)

    return [members for _, members in sorted(groups.items())]


def junction_multiplicities(cycles, group):
    """
    Count how many times each junction is used by a group of cycles.

    A junction listed by two different cycles describes one adjacency,
    not two, so counts are taken as the highest per-cycle count rather
    than summed. A cycle that genuinely uses a junction twice still
    contributes two.
    """
    counts = {}

    for index in group:
        for junction, count in Counter(cycles[index]["junctions"]).items():
            counts[junction] = max(counts.get(junction, 0), count)

    return counts


def segment_multiplicities(junction_counts):
    """
    Derive how many times each segment must be traversed.

    A closed molecule alternates segment and junction, so a segment's
    traversal count equals the junction count at each of its two ends.

    Returns:
        multiplicities  region -> traversal count
        imbalances      (region, left_count, right_count) for segments
                        whose ends disagree, which no single molecule
                        can explain
    """
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
    """Return the other end of the same segment."""
    return RIGHT if end == LEFT else LEFT


def build_junction_uses(junction_counts):
    """
    Expand junction multiplicities into individual uses.

    Returns:
        uses    use_id -> (end_a, end_b)
        slots   segment end -> list of use_ids arriving there
    """
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
    """Return the end at the far side of a junction use."""
    end_a, end_b = uses[use_id]

    return end_b if end_a == end else end_a


def initial_transitions(slots, multiplicities):
    """
    Pair every left-end slot of a segment with a right-end slot.

    Each pair is one traversal of that segment: it says which junction
    use the walk leaves by, having arrived on the other one. Any pairing
    is a valid starting point; merge_transition_walks then re-pairs to
    splice separate walks together.
    """
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
    """
    Split every junction use into closed alternating walks.

    The walk arrives at a segment end along one junction use, traverses
    that segment to its far end, and leaves along the use the transition
    system pairs it with. A junction use is crossed once no matter which
    of its two ends the walk happens to enter from, so consumption is
    tracked per use: entering from the other end would only retrace the
    same molecule in reverse.

    Returns a list of walks, each a list of
    (region, strand, left_use, right_use) traversals.
    """
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
    """
    Re-pair transitions until every junction use lies on one walk.

    Two walks that both traverse the same segment are spliced by
    swapping that segment's pairing, which is the constructive half of
    Kotzig's theorem. Each splice removes one walk, so this terminates.

    Returns the single remaining walk, or None when the segment ends
    cannot be joined into one molecule.
    """
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
    """
    Reconstruct one closed molecule from a group's junctions.

    Returns strand-aware region tokens with the first token repeated to
    mark the closure, or [] when the junctions cannot form a single
    molecule.
    """
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
                "nodes": nodes,
                "junctions": cycle_junctions(nodes)
            })

    return cycles


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
    Reconstruct big ecDNA cycles from the small cycles in summary.txt.

    Small cycles are grouped only when they share a junction, which is
    direct read evidence that they describe the same molecule. Each
    group is resolved on a segment-end breakpoint graph, so a single
    molecule may traverse a segment in both orientations.

    A group whose segment ends do not balance cannot be explained by one
    molecule - typically two species sharing a junction - so its member
    cycles are reported individually instead of being fused.

    Region depth is calculated from the BAM file and compared with how
    often each region is traversed. With require_depth_support=False
    every reconstruction is written and marked depth-consistent or not;
    with True, depth-inconsistent reconstructions are omitted.
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

    groups = build_junction_groups(
        cycles
    )

    logging.info(
        "Grouped %d small cycles into %d junction-linked group(s)",
        len(cycles),
        len(groups)
    )

    output_records = []

    def add_record(walk, source_ids):
        """Depth-check one reconstruction and queue it for output."""
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

        output_records.append({
            "bigcycle_id": len(output_records) + 1,
            "source_cycles": ",".join(source_ids),
            "total_size": path_size(walk),
            "fragment_count": (
                len(walk) - 1
                if walk[0] == walk[-1]
                else len(walk)
            ),
            "depth_consistent": "yes" if valid else "no",
            "region_depths": format_region_values(
                used_regions,
                region_depths,
                numeric_format="{:.2f}"
            ),
            "region_copy_numbers": format_region_values(
                used_regions,
                copy_numbers
            ),
            "region_occurrences": format_region_values(
                used_regions,
                occurrences
            ),
            "depth_violations": format_violations(violations),
            "path": " -> ".join(walk)
        })

    def add_members_individually(group, reason):
        """
        Report a group's cycles separately.

        Used when a group cannot be explained by one molecule, so the
        evidence is still reported rather than dropped.
        """
        logging.info(
            "Reporting cycles %s individually: %s",
            ",".join(cycles[index]["id"] for index in group),
            reason
        )

        for index in group:
            add_record(
                cycles[index]["nodes"],
                [cycles[index]["id"]]
            )

    for group in groups:
        source_ids = [
            cycles[index]["id"]
            for index in group
        ]

        if len(group) == 1:
            add_record(cycles[group[0]]["nodes"], source_ids)
            continue

        junction_counts = junction_multiplicities(cycles, group)

        multiplicities, imbalances = segment_multiplicities(
            junction_counts
        )

        if imbalances:
            add_members_individually(
                group,
                "segment ends do not balance ({})".format(
                    ", ".join(
                        "{} has {} junction(s) on its left end "
                        "and {} on its right".format(region, left, right)
                        for region, left, right in imbalances
                    )
                )
            )
            continue

        walk = reconstruct_big_cycle(junction_counts, multiplicities)

        if len(walk) < 2:
            add_members_individually(
                group,
                "the junctions do not form a single closed molecule"
            )
            continue

        logging.info(
            "Merged cycles %s into a %d-fragment big cycle",
            ",".join(source_ids),
            len(walk) - 1
        )

        add_record(walk, source_ids)

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