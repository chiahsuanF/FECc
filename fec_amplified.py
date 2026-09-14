#!/usr/bin/env python3

import argparse, os, sys, random, logging,  pysam, csv, time
from bisect import bisect_left, bisect_right
import numpy as np
from fec_login import *

# -----------------------
# 2. Identify Amplified Regions // find hight coverage region(correct)
# -----------------------

def merge_regions(regions):
    merged = []
    for chrom, start, end in sorted(regions):
        if not merged or chrom != merged[-1][0] or start > merged[-1][2]:
            merged.append([chrom, start, end])
        else:
            merged[-1][2] = max(merged[-1][2], end)
    return [(c, s, e) for c, s, e in merged]
def identify_amplified_regions(bam_path, outdir, sample,
                               min_len=100, num_samples=100): #100
    t0 = tic("Identifying amplified regions")
    bam = pysam.AlignmentFile(bam_path,'rb')
    
    allowed = {f"chr{i}" for i in range(1, 23)}
    allowed.update({"chrX", "chrY"})

    contigs = [
        (chrom, length)
        for chrom, length in zip(bam.references, bam.lengths)
        if chrom in allowed
    ]


    covs = []
    for _ in range(num_samples):
        chrom, L = random.choice(contigs)
        pos = random.randrange(L)
        for col in bam.pileup(chrom, pos, pos+1, truncate=True):
            if col.n > 0:
                mq60 = sum(p.alignment.mapping_quality >= 60 for p in col.pileups)
                if mq60 / col.n >= 0.95:
                    covs.append(col.n)
            break
    bam.close()
    if not covs:
        logging.error("No high‐quality coverage samples; aborting")
        sys.exit(1)

    arr = np.array(covs)
    amp_thresh = arr.mean() + 3 * arr.std()
    logging.info(f"Amplification threshold = {amp_thresh:.1f}")

    # scan genome
    regs = []
    bam = pysam.AlignmentFile(bam_path,'rb')
    for chrom, L in contigs:
        start = None
        for col in bam.pileup(chrom, 0, L, truncate=True):
            if col.n >= amp_thresh:
                if start is None:
                    start = col.pos
            else:
                if start is not None and (col.pos - start) >= min_len:
                    regs.append((chrom, start, col.pos - 1))
                start = None
        if start is not None:
            regs.append((chrom, start, L))
    bam.close()

    fn = os.path.join(outdir, f"{sample}_amplified_regions.txt")
    with open(fn,'w') as f:
        f.write("chrom\tstart\tend\n")
        for c,s,e in regs:
            f.write(f"{c}\t{s}\t{e}\n")
    logging.info(f"Wrote {len(regs)} amplified regions → {fn}")
    toc(t0, "Identifying amplified regions")
    return regs, amp_thresh