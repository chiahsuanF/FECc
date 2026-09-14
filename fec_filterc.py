#!/usr/bin/env python3

import os, logging, pysam, csv, time
from bisect import bisect_left, bisect_right
from fec_login import *

# -----------------------
# 7. Filter & Write Cycles
# -----------------------
def check_region_quality(region, bam_path, qthr=5, fthr=0.01):
    chrom,s,e = region
    bam = pysam.AlignmentFile(bam_path,'rb')
    tot=high=0
    for r in bam.fetch(chrom, s, e):
        tot+=1
        if r.mapping_quality>=qthr:
            high+=1
    bam.close()
    return tot>0 and (high/tot)>=fthr

def filter_and_write_cycles(
    cycles,
    regs,
    bam_path,
    outdir,
    sample
):

    t0 = tic("Filtering cycles")

    quality_cache = {}

    def good(node):

        idx = node[0]

        if idx not in quality_cache:

            quality_cache[idx] = check_region_quality(
                regs[idx],
                bam_path
            )

        return quality_cache[idx]

    filtered = [
        c
        for c in cycles
        if all(
            good(node)
            for node in c[:-1]
        )
    ]

    fn = os.path.join(
        outdir,
        f"08_{sample}_filtered_cycles.txt"
    )

    with open(fn, "w") as f:

        f.write("cycle_nodes\n")

        for c in filtered:

            f.write(
                ",".join(map(str, c))
                + "\n"
            )

    logging.info(
        f"Wrote {len(filtered):,} filtered cycles -> {fn}"
    )

    toc(t0, "Filtering cycles")

    return filtered

