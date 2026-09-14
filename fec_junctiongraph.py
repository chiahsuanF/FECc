#!/usr/bin/env python3

import os, logging
from collections import defaultdict
from functools import lru_cache
from fec_login import *


# -----------------------
# 5. Build Junction Graph
# -----------------------

MIN_SA_SUPPORT = 2
MIN_SA_RATIO = 0.2

def build_graph( regs, bps, outdir, sample, global_depth, min_mq=30):
    t0 = tic("Building junction graph")
    
    if global_depth <= 0:
        logging.error(
            "Global depth must be > 0"
        )
        return defaultdict(dict)

    # -----------------------------------
    # region index
    # -----------------------------------
    per_chr = defaultdict(list)

    for idx, (chrom, start, end) in enumerate(regs):
        per_chr[chrom].append((start, end, idx))

    for lst in per_chr.values():
        lst.sort()
    # -----------------------------------
    # cached binary search
    # -----------------------------------
    @lru_cache(maxsize=1000000)
    def find_idx(chrom, pos):

        lst = per_chr.get(chrom)

        if not lst:
            return None

        lo = 0
        hi = len(lst) - 1

        while lo <= hi:

            mid = (lo + hi) // 2

            start, end, idx = lst[mid]

            if start <= pos <= end:
                return idx

            if pos < start:
                hi = mid - 1
            else:
                lo = mid + 1

        return None


    # -----------------------------------
    # collect edge support
    # -----------------------------------

    edge_reads = defaultdict(set)

    processed = 0
    kept = 0

    for q, c1, b1, mq1, strand1, c2, b2, mq2, strand2 in bps:

        processed += 1

        if mq1 < min_mq or mq2 < min_mq:
            continue

        i1 = find_idx(c1, b1)
        if i1 is None:
            continue

        i2 = find_idx(c2, b2)
        if i2 is None:
            continue

        if i1 == i2:
            continue

        kept += 1

        edge_key = (i1, strand1, i2, strand2)

        edge_reads[edge_key].add(q)

    logging.info(
        f"Processed {processed:,} breakpoints "
        f"(kept {kept:,})"
    )

    edge_support = {
        edge: len(reads)
        for edge, reads in edge_reads.items()
    }

    # -----------------------------------
    # write edges
    # -----------------------------------
    fn = os.path.join(
        outdir,
        f"05_{sample}_graph_edges.txt"
    )

    graph = defaultdict(dict)

    with open(fn, "w") as f:

        f.write(
            "from\tto\tcount\tdepth\tratio\n"
        )

        for (i1, s1, i2, s2), sa_support in edge_support.items():

            if sa_support <= MIN_SA_SUPPORT:
                continue

            ratio = sa_support / global_depth

            if ratio <= MIN_SA_RATIO:
                continue

            graph[(i1, s1)][(i2, s2)] = sa_support

            f.write(
                f"{i1}{s1}\t"
                f"{i2}{s2}\t"
                f"{sa_support}\t"
                f"{global_depth:.2f}\t"
                f"{ratio:.3f}\n"
            )

    logging.info(
        f"Wrote graph edges -> {fn}"
    )

    toc(t0, "Building junction graph")

    return graph
