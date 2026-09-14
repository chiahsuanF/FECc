#!/usr/bin/env python3

import logging
import pysam
from collections import Counter

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
# Selecting big cycles
# ---------------------------------------------------------------
#
# fec_cycled's search never revisits a region inside one path, so a
# small cycle always lists distinct regions. A larger cycle therefore
# explains a smaller one whenever it covers all of the smaller one's
# regions: the small cycle is a shortcut across the same structure, not
# separate evidence. The big cycle for an event is the small cycle whose
# region set no other cycle covers.


def cycle_regions(cycle):
    """Return the set of strandless regions a cycle visits."""
    return {
        strandless_node(node)
        for node in cycle["nodes"][:-1]
    }


def select_big_cycles(cycles):
    """
    Keep the cycles whose region set is not contained in another's.

    Returns a list of (cycle, subsumed_cycles) pairs, ordered by
    descending region count so the largest reconstruction comes first.
    Ties on an identical region set keep the cycle with more fragments,
    then the lower id, so the choice does not depend on input order.
    """
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
                "nodes": nodes
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
    Build big ecDNA cycles from the small cycles in summary.txt.

    A small cycle whose regions are all covered by a larger cycle is a
    shortcut across the same structure rather than separate evidence, so
    the big cycle for an event is the small cycle no other cycle covers.
    The cycles it explains are listed alongside it.

    Region depth is read from the BAM and compared with how often each
    region is traversed. With require_depth_support=False every cycle is
    written and marked depth-consistent or not; with True, the
    depth-inconsistent ones are omitted.
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

    chosen = select_big_cycles(cycles)

    logging.info(
        "Selected %d big cycle(s) from %d small cycle(s)",
        len(chosen),
        len(cycles)
    )

    output_records = []

    for cycle, subsumed in chosen:
        walk = cycle["nodes"]

        source_ids = [cycle["id"]] + [
            other["id"]
            for other in subsumed
        ]

        if subsumed:
            logging.info(
                "Big cycle from cycle %s covers cycle(s) %s",
                cycle["id"],
                ",".join(other["id"] for other in subsumed)
            )

        valid, occurrences, violations = validate_walk_depth(
            walk=walk,
            copy_numbers=copy_numbers,
            copy_tolerance=copy_tolerance
        )

        if require_depth_support and not valid:
            logging.warning(
                "Skipping depth-inconsistent cycle %s: %s",
                cycle["id"],
                violations
            )
            continue

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