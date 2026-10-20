#!/usr/bin/env python3

import os, logging, pysam, random, statistics
from concurrent.futures import ThreadPoolExecutor

# -----------------------
# 1. Depth & MinSupport
# -----------------------

APPROX_CENTROMERE_MB = {
    "chr1": 123, "chr2": 94, "chr3": 91, "chr4": 50, "chr5": 49,
    "chr6": 60, "chr7": 60, "chr8": 45, "chr9": 43, "chr10": 40,
    "chr11": 53, "chr12": 36, "chr13": 18, "chr14": 17, "chr15": 19,
    "chr16": 37, "chr17": 25, "chr18": 18, "chr19": 26, "chr20": 28,
    "chr21": 12, "chr22": 15, "chrX": 61, "chrY": 10,
}

AVOID_MARGIN = 5_000_000


def _usable_read(read):
    return not (
        read.is_unmapped
        or read.is_secondary
        or read.is_qcfail
        or read.is_duplicate
    )


def _window_depth(args):
    bam_path, chrom, start, end = args

    bam = pysam.AlignmentFile(bam_path, "rb")

    coverage = bam.count_coverage(
        chrom,
        start,
        end,
        quality_threshold=0,
        read_callback=_usable_read
    )

    bam.close()

    return sum(sum(base) for base in coverage) / float(end - start)


def _near_centromere_or_telomere(chrom, length, start, end):
    if start < AVOID_MARGIN or end > length - AVOID_MARGIN:
        return True

    centre = APPROX_CENTROMERE_MB.get(chrom)

    if centre is None:
        return False

    centre *= 1_000_000

    return (
        start < centre + AVOID_MARGIN
        and end > centre - AVOID_MARGIN
    )


# -----------------------------
# sampled mean depth
# -----------------------------
def sample_region_depths(bam_path, num_regions=100, window=10000,
                         seed=1, threads=4, max_attempts=10000):

    bam = pysam.AlignmentFile(bam_path, "rb")

    allowed = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}

    contigs = [
        (chrom, length)
        for chrom, length in zip(bam.references, bam.lengths)
        if chrom in allowed and length > 2 * AVOID_MARGIN + window
    ]

    bam.close()

    if not contigs:
        raise ValueError("No recognized chromosomes found in BAM")

    rng = random.Random(seed)
    weights = [length for _, length in contigs]

    depths = []
    attempts = 0
    avoided = 0
    empty = 0

    with ThreadPoolExecutor(max_workers=threads) as pool:
        while len(depths) < num_regions and attempts < max_attempts:
            tasks = []

            while (
                len(tasks) < num_regions - len(depths)
                and attempts < max_attempts
            ):
                attempts += 1

                chrom, length = rng.choices(contigs, weights=weights)[0]
                start = rng.randrange(length - window)
                end = start + window

                if _near_centromere_or_telomere(chrom, length, start, end):
                    avoided += 1
                    continue

                tasks.append((bam_path, chrom, start, end))

            for depth in pool.map(_window_depth, tasks):
                if depth > 0:
                    depths.append(depth)
                else:
                    empty += 1

    logging.info(
        "Depth sampling: attempts=%d regions=%d/%d window=%d "
        "near_centromere_or_telomere=%d zero_depth=%d",
        attempts,
        len(depths),
        num_regions,
        window,
        avoided,
        empty
    )

    return depths


# -----------------------------
# main function
# -----------------------------
def estimate_depth_and_minsupp(bam_path, min_support, depth_override,
                               outdir, sample, threads=4,
                               num_regions=100, window=10000, seed=1):

    logging.info("Estimating depth & minsupp from %d random regions",
                 num_regions)

    stats = None

    if depth_override is not None:
        depth = depth_override
    else:
        depths = sorted(sample_region_depths(
            bam_path,
            num_regions=num_regions,
            window=window,
            seed=seed,
            threads=threads
        ))

        trimmed = depths[1:-1] if len(depths) > 2 else depths

        if not trimmed:
            raise ValueError("No covered regions found for depth estimate")

        depth = statistics.mean(trimmed)

        stats = {
            "sd": statistics.stdev(trimmed) if len(trimmed) > 1 else 0.0,
            "median": statistics.median(trimmed),
            "min": trimmed[0],
            "max": trimmed[-1],
            "used": len(trimmed),
            "sampled": len(depths),
        }

    if min_support is None:
        #minsupp = round(depth / 10) + 2
        minsupp = round(depth / 2)
    else:
        minsupp = min_support

    os.makedirs(outdir, exist_ok=True)
    fn = os.path.join(outdir, f"{sample}_depth_minsupp.txt")

    with open(fn, "w") as f:
        f.write(f"Estimated Depth: {depth:.2f}\n")
        f.write(f"MinSupport: {minsupp}\n")

        if stats is not None:
            f.write(f"Depth SD: {stats['sd']:.2f}\n")
            f.write(f"Depth Median: {stats['median']:.2f}\n")
            f.write(f"Depth Min (trimmed): {stats['min']:.2f}\n")
            f.write(f"Depth Max (trimmed): {stats['max']:.2f}\n")
            f.write(
                f"Regions used: {stats['used']} of {stats['sampled']} "
                f"(window {window} bp, max/min removed)\n"
            )

    logging.info(f"Wrote depth & minsupp → {fn}")

    return depth, minsupp
