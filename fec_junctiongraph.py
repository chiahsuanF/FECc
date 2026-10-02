#!/usr/bin/env python3

import os, logging, pysam
from collections import defaultdict
from functools import lru_cache
from fec_login import *


# -----------------------
# 5. Build Junction Graph
# -----------------------

MIN_SA_SUPPORT = 2
MIN_SA_RATIO = 0.05#

# -----------------------------------
# local depth: mean coverage of the given fragments
# -----------------------------------
SKIP_FLAGS = 0x4 | 0x100 | 0x200 | 0x400

def aligned_bases(bam, chrom, s, e, min_mq):

    total = 0

    for read in bam.fetch(chrom, s, e):

        if read.flag & SKIP_FLAGS or read.mapping_quality < min_mq:
            continue

        if read.reference_start >= s and read.reference_end <= e:
            # M, =, X operations
            stats = read.get_cigar_stats()[0]
            total += int(stats[0] + stats[7] + stats[8])
            continue

        for bs, be in read.get_blocks():
            if be <= s or bs >= e:
                continue
            total += min(be, e) - max(bs, s)

    return total


def fragment_depths(regs, indices, bam_path, min_mq=30):

    depths = {}

    with pysam.AlignmentFile(bam_path, "rb") as bam:

        chrom_len = dict(zip(bam.references, bam.lengths))

        for idx in sorted(indices):

            chrom, start, end = regs[idx]

            s = max(0, int(start))
            e = min(int(end) + 1, chrom_len.get(chrom, 0))

            if e <= s:
                depths[idx] = 0.0
            else:
                depths[idx] = (
                    aligned_bases(bam, chrom, s, e, min_mq) / (e - s)
                )

            logging.info(
                f"Local depth fragment {idx} "
                f"{chrom}:{start}-{end} "
                f"(len={end - start + 1:,}) "
                f"depth={depths[idx]:.2f}"
            )

    return depths


def build_graph( regs, bps, bam_path, outdir, sample, min_mq=30):
    t0 = tic("Building junction graph")

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
    # depth only for fragments on an edge that can pass
    # -----------------------------------
    used = {
        i
        for (i1, _, i2, _), sa_support in edge_support.items()
        if sa_support > MIN_SA_SUPPORT
        for i in (i1, i2)
    }

    logging.info(
        f"Computing local depth for {len(used):,} "
        f"of {len(regs):,} fragments"
    )

    local_depth = fragment_depths(regs, used, bam_path, min_mq)

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

            depth = min(local_depth[i1], local_depth[i2])
            #depth = local_depth[i1]

            if depth <= 0:
                logging.warning(
                    f"Skip edge {i1}{s1}->{i2}{s2}: "
                    f"fragment depth is 0"
                )
                continue

            ratio = sa_support / depth

            if ratio <= MIN_SA_RATIO:
                continue

            graph[(i1, s1)][(i2, s2)] = sa_support

            f.write(
                f"{i1}{s1}\t"
                f"{i2}{s2}\t"
                f"{sa_support}\t"
                f"{depth:.2f}\t"
                f"{ratio:.3f}\n"
            )

    logging.info(
        f"Wrote graph edges -> {fn}"
    )

    toc(t0, "Building junction graph")

    return graph
