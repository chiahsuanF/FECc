#!/usr/bin/env python3

import os, logging, csv, time
from collections import defaultdict
from bisect import bisect_left, bisect_right
from fec_login import *

# -----------------------
# 4. Refine Regions by High-Q Clusters // 
# -----------------------
"""
def refine_regions_by_highQ_breakpoints(regions, bps, outdir, sample, window=120000):
    t0 = tic("Refining regions by high-Q breakpoints")
    new_regs = []
    per_chr = defaultdict(list)
    for c,s,e in regions:
        per_chr[c].append((s,e))
    for cl in per_chr.values():
        cl.sort()

    def find_block(c,pos):
        lst = per_chr.get(c,[])
        lo,hi=0,len(lst)-1
        while lo<=hi:
            mid=(lo+hi)//2
            s,e = lst[mid]
            if s<=pos< e: return (s,e) # if s<=pos<=e: return (s,e)
            if pos < s: hi=mid-1
            else: lo=mid+1
        return None

    for c,s,e in regions:
        pts = []
        for q,c1,b1,m1,c2,b2,m2 in bps:
            if c1==c and s<=b1<=e: pts.append((b1,m1))
            if c2==c and s<=b2<=e: pts.append((b2,m2))
        if not pts:
            new_regs.append((c,s,e))
            continue
        pts.sort()
        clusters, curr = [], [pts[0]]
        for pos,mq in pts[1:]:
            if pos - curr[-1][0] <= window:
                curr.append((pos,mq))
            else:
                clusters.append(curr)
                curr=[(pos,mq)]
        clusters.append(curr)
        reps = sorted({ max(cl, key=lambda x:x[1])[0] for cl in clusters })
        cuts = [s] + reps + [e]
        for a,b in zip(cuts, cuts[1:]):
            if a<b:
                new_regs.append((c,a,b))

    fn = os.path.join(outdir, f"04_{sample}_refined_regions.txt")
    with open(fn,'w') as f:
        f.write("chrom\tstart\tend\n")
        for c,s,e in new_regs:
            f.write(f"{c}\t{s}\t{e}\n")
    logging.info(f"Wrote {len(new_regs)} refined regions → {fn}")
    toc(t0, "Refining regions by high-Q breakpoints")
    return new_regs
"""

def refine_regions_by_highQ_breakpoints( regions, bps, outdir, sample, window=120000):

    t0 = tic("Refining regions by high-Q breakpoints")

    chr_pts = defaultdict(list)

    for q,c1,b1,m1,strand1,c2,b2,m2,strand2 in bps:
        chr_pts[c1].append((b1,m1))
        chr_pts[c2].append((b2,m2))

    for chrom in chr_pts:
        chr_pts[chrom].sort()

    for chrom in chr_pts:
        chr_pts[chrom].sort()

    chr_pos = {}

    for chrom, pts in chr_pts.items():
        chr_pos[chrom] = [x[0] for x in pts]

    new_regs = []

    for chrom,start,end in regions:

        if chrom not in chr_pts:
            new_regs.append((chrom,start,end))
            continue

        lo = bisect_left(
            chr_pos[chrom],
            start
        )

        hi = bisect_right(
            chr_pos[chrom],
            end
        )

        pts = chr_pts[chrom][lo:hi]

        if not pts:
            new_regs.append((chrom,start,end))
            continue

        clusters = []
        curr = [pts[0]]

        for pos,mq in pts[1:]:

            if pos - curr[-1][0] <= window:
                curr.append((pos,mq))
            else:
                clusters.append(curr)
                curr = [(pos,mq)]

        clusters.append(curr)

        reps = sorted(
            max(cl,key=lambda x:x[1])[0]
            for cl in clusters
        )

        cuts = [start] + reps + [end]

        for a,b in zip(cuts,cuts[1:]):
            if a < b:
                new_regs.append(
                    (chrom,a,b)
                )

    fn = os.path.join(
        outdir,
        f"04_{sample}_refined_regions.txt"
    )

    with open(fn,"w") as f:
        f.write("chrom\tstart\tend\n")

        for c,s,e in new_regs:
            f.write(f"{c}\t{s}\t{e}\n")

    logging.info(
        f"Wrote {len(new_regs)} refined regions → {fn}"
    )

    toc(
        t0,
        "Refining regions by high-Q breakpoints"
    )

    return new_regs