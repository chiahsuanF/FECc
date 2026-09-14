#!/usr/bin/env python3

import  os,  logging,  pysam, csv
from bisect import bisect_left, bisect_right
from concurrent.futures import ThreadPoolExecutor
from fec_login import *

# -----------------------
# 1. Depth & MinSupport
# -----------------------
def _depth_worker(args):
    bam_path, chrom, start, end = args

    bam = pysam.AlignmentFile(bam_path, "rb")

    total = 0
    count = 0

    for col in bam.pileup(
        chrom,
        start,
        end,
        truncate=True,
        stepper="all"
    ):
        total += col.n
        count += 1

    bam.close()

    return total, count


# -----------------------------
# parallel mean depth
# -----------------------------
def compute_mean_depth_mt(bam_path, threads=4):

    bam = pysam.AlignmentFile(bam_path, "rb")

    contigs = list(zip(bam.references, bam.lengths))
    bam.close()

    tasks = [(bam_path, c, 0, l) for c, l in contigs]

    total_depth = 0
    total_pos = 0

    with ThreadPoolExecutor(max_workers=threads) as pool:
        results = pool.map(_depth_worker, tasks)

        for t, c in results:
            total_depth += t
            total_pos += c

    return total_depth / total_pos if total_pos > 0 else 0


# -----------------------------
# main function 
# -----------------------------
def estimate_depth_and_minsupp(bam_path, min_support, depth_override,
                               outdir, sample, threads=4):

    logging.info("Estimating depth & minsupp (multi-threaded)")

    # ---- depth ----
    if depth_override is not None:
        depth = depth_override
    else:
        depth = compute_mean_depth_mt(bam_path, threads=threads)

    # ---- minsupp ----
    if min_support is None:
        #minsupp = round(depth / 10) + 2
        minsupp = round(depth / 2)
    else:
        minsupp = min_support

    # ---- output ----
    os.makedirs(outdir, exist_ok=True)
    fn = os.path.join(outdir, f"{sample}_depth_minsupp.txt")

    with open(fn, "w") as f:
        f.write(f"Estimated Depth: {depth:.2f}\n")
        f.write(f"MinSupport: {minsupp}\n")

    logging.info(f"Wrote depth & minsupp → {fn}")

    return depth, minsupp