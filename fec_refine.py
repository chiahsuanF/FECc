#!/usr/bin/env python3

import os, logging
from collections import defaultdict
from bisect import bisect_left, bisect_right
from fec_login import *
from fec_amplified import merge_regions

# -----------------------
# 3b. Extend regions to anchored breakpoints
# -----------------------


def extend_regions_to_anchored_breakpoints( regions, bps, outdir, sample, min_reads=3, min_mq=30, partner_mq=0, max_extension=50000, window=1000, partner_min_len=0):

    t0 = tic("Extending regions to anchored breakpoints")

    regions = sorted(regions)
    per_chr = defaultdict(list)

    for idx, (chrom, start, end) in enumerate(regions):
        per_chr[chrom].append((start, end, idx))

    starts = {chrom: [x[0] for x in lst] for chrom, lst in per_chr.items()}

    def inside(chrom, pos):
        lst = per_chr.get(chrom)
        if not lst:
            return False
        i = bisect_right(starts[chrom], pos + window) - 1
        return any(
            lst[k][0] - window <= pos <= lst[k][1] + window
            for k in (i - 1, i)
            if 0 <= k < len(lst)
        )

    def nearest_edge(chrom, pos):
        lst = per_chr.get(chrom)
        if not lst:
            return None
        i = bisect_right(starts[chrom], pos)
        best = None
        if i > 0 and 0 < pos - lst[i - 1][1] <= max_extension:
            best = (pos - lst[i - 1][1], lst[i - 1][2], "R")
        if i < len(lst) and 0 < lst[i][0] - pos <= max_extension:
            if best is None or lst[i][0] - pos < best[0]:
                best = (lst[i][0] - pos, lst[i][2], "L")
        return best[1:] if best else None

    def edge_at(chrom, pos):
        lst = per_chr.get(chrom)
        if not lst:
            return None
        i = bisect_right(starts[chrom], pos + window) - 1
        for k in (i - 1, i):
            if 0 <= k < len(lst):
                if abs(pos - lst[k][0]) <= window:
                    return (lst[k][2], "L")
                if abs(pos - lst[k][1]) <= window:
                    return (lst[k][2], "R")
        return None

    candidates = defaultdict(list)
    anchored = defaultdict(set)

    for q,c1,b1,m1,strand1,c2,b2,m2,strand2 in bps:

        if not junction_mq_passes(m1, m2, min_mq, partner_mq):
            continue

        for chrom, pos in ((c1, b1), (c2, b2)):
            edge = edge_at(chrom, pos)
            if edge is not None:
                anchored[edge].add(q)

        in1, in2 = inside(c1, b1), inside(c2, b2)

        if in1 == in2:
            continue

        chrom, pos = (c2, b2) if in1 else (c1, b1)
        edge = nearest_edge(chrom, pos)

        if edge is None:
            continue

        _, start, end = regions[edge[0]]

        if not partner_side_allowed(m1, m2, min_mq, end - start, end - start, partner_min_len):
            continue

        candidates[edge].append((pos, q))

    new_edges = {}

    for (idx, side), pts in candidates.items():
        if len(anchored[(idx, side)]) >= min_reads:
            continue

        pts.sort()
        clusters = []
        curr = [pts[0]]

        for pt in pts[1:]:
            if pt[0] - curr[-1][0] <= window:
                curr.append(pt)
            else:
                clusters.append(curr)
                curr = [pt]

        clusters.append(curr)

        supported = [
            cl for cl in clusters
            if len({x[1] for x in cl}) >= min_reads
        ]

        if not supported:
            continue

        if side == "R":
            cl = max(supported, key=lambda c: c[-1][0])
            new_edges[(idx, side)] = (cl[-1][0], len({x[1] for x in cl}))
        else:
            cl = min(supported, key=lambda c: c[0][0])
            new_edges[(idx, side)] = (cl[0][0], len({x[1] for x in cl}))

    extended = []

    for idx, (chrom, start, end) in enumerate(regions):
        new_start, new_end = start, end

        prev_end = max(
            (e for c, s_, e in regions[:idx] if c == chrom),
            default=None
        )
        next_start = regions[idx + 1][1] if (
            idx + 1 < len(regions) and regions[idx + 1][0] == chrom
        ) else None

        if (idx, "L") in new_edges and prev_end is not None:
            pos, reads = new_edges[(idx, "L")]
            new_edges[(idx, "L")] = (max(pos, prev_end + 1), reads)

        if (idx, "R") in new_edges and next_start is not None:
            pos, reads = new_edges[(idx, "R")]
            new_edges[(idx, "R")] = (min(pos, next_start - 1), reads)

        if (idx, "L") in new_edges:
            new_start, reads = new_edges[(idx, "L")]
            logging.info(
                f"Extended {chrom}:{start}-{end} left to {new_start} "
                f"({start - new_start} bp, {reads} anchored reads)"
            )

        if (idx, "R") in new_edges:
            new_end, reads = new_edges[(idx, "R")]
            logging.info(
                f"Extended {chrom}:{start}-{end} right to {new_end} "
                f"({new_end - end} bp, {reads} anchored reads)"
            )

        extended.append((chrom, new_start, new_end))

    extended = merge_regions(extended)

    fn = os.path.join(
        outdir,
        f"03b_{sample}_extended_regions.txt"
    )

    with open(fn,"w") as f:
        f.write("chrom\tstart\tend\n")

        for c,s,e in extended:
            f.write(f"{c}\t{s}\t{e}\n")

    logging.info(
        f"Extended {len(new_edges)} region edge(s); wrote "
        f"{len(extended)} regions → {fn}"
    )

    toc(t0, "Extending regions to anchored breakpoints")

    return extended

# -----------------------
# 3a. Split bridged gaps that carry a junction
# -----------------------

def split_regions_at_junction_gaps( regions, gaps, bps, outdir, sample, min_reads=3, min_mq=30, partner_mq=0, window=1000):

    t0 = tic("Splitting bridged gaps at junctions")

    pts = defaultdict(list)

    for q,c1,b1,m1,strand1,c2,b2,m2,strand2 in bps:
        if not junction_mq_passes(m1, m2, min_mq, partner_mq):
            continue
        pts[c1].append((b1, q))
        pts[c2].append((b2, q))

    for chrom in pts:
        pts[chrom].sort()

    pos = {chrom: [x[0] for x in lst] for chrom, lst in pts.items()}

    cuts = defaultdict(list)

    for chrom, gap_start, gap_end in gaps:
        if chrom not in pts:
            continue
        lo = bisect_left(pos[chrom], gap_start - window)
        hi = bisect_right(pos[chrom], gap_end + window)
        reads = {q for _, q in pts[chrom][lo:hi]}
        if len(reads) >= min_reads:
            cuts[chrom].append((gap_start, gap_end))
            logging.info(
                f"Splitting at bridged gap {chrom}:{gap_start}-{gap_end} "
                f"({len(reads)} junction reads)"
            )

    out = []

    for chrom, start, end in regions:
        pieces = [(start, end)]
        for gap_start, gap_end in sorted(cuts.get(chrom, [])):
            nxt = []
            for a, b in pieces:
                if a < gap_start and gap_end < b:
                    nxt += [(a, gap_start), (gap_end, b)]
                else:
                    nxt.append((a, b))
            pieces = nxt
        out += [(chrom, a, b) for a, b in pieces]

    fn = os.path.join(
        outdir,
        f"03a_{sample}_gap_split_regions.txt"
    )

    with open(fn,"w") as f:
        f.write("chrom\tstart\tend\n")

        for c,s,e in out:
            f.write(f"{c}\t{s}\t{e}\n")

    logging.info(
        f"Split {sum(map(len, cuts.values()))} bridged gap(s); wrote "
        f"{len(out)} regions → {fn}"
    )

    toc(t0, "Splitting bridged gaps at junctions")

    return out

# -----------------------
# 4. Refine Regions by High-Q Clusters // 
# -----------------------

def refine_regions_by_highQ_breakpoints( regions, bps, outdir, sample, window=1000, min_cluster_reads=3, min_piece=1000):

    t0 = tic("Refining regions by high-Q breakpoints")

    chr_pts = defaultdict(list)

    for q,c1,b1,m1,strand1,c2,b2,m2,strand2 in bps:
        chr_pts[c1].append((b1,m1,q))
        chr_pts[c2].append((b2,m2,q))

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

        for pt in pts[1:]:

            if pt[0] - curr[-1][0] <= window:
                curr.append(pt)
            else:
                clusters.append(curr)
                curr = [pt]

        clusters.append(curr)

        reps = sorted(
            rep
            for rep in (
                max(cl,key=lambda x:x[1])[0]
                for cl in clusters
                if len({x[2] for x in cl}) >= min_cluster_reads
            )
            if start + min_piece < rep < end - min_piece
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
