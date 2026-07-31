#!/usr/bin/env python3

import os, logging, csv,time
from collections import defaultdict
from bisect import bisect_left, bisect_right
from fec_login import *

# -----------------------
# 5. Build Junction Graph  
# -----------------------
def build_graph(regs, bps, amp_thresh, outdir, sample):
    t0 = tic("Building junction graph")
    per_chr = defaultdict(list)
    for idx,(c,s,e) in enumerate(regs):
        per_chr[c].append((s,e,idx))
    for cl in per_chr.values():
        cl.sort()

    def find_idx(c,pos):
        lst=per_chr.get(c,[])
        lo,hi=0,len(lst)-1
        while lo<=hi:
            mid=(lo+hi)//2
            s,e,i = lst[mid]
            if s<=pos<=e: return i # if s<=pos<e: return i
            if pos < s: hi=mid-1
            else: lo=mid+1
        return None

    edge_support = defaultdict(list)
    for q,c1,b1,m1,c2,b2,m2 in bps:
        i1=find_idx(c1,b1)
        i2=find_idx(c2,b2)
        if i1 is None or i2 is None or i1==i2: continue
        #edge_support[(i1,i2)].append((q,b1,b2))
        edge_support[(i1,i2)].append((q, c1, b1, c2, b2))

    read_thr = max(2, int(amp_thresh/3))
    #read_thr = 2
    fn = os.path.join(outdir, f"05_{sample}_graph_edges.txt")
    with open(fn,'w') as f:
        f.write("from\tto\tcount\treads_with_coords\n")
        for (i1,i2), recs in edge_support.items():
            cnt=len(recs)
            if cnt<read_thr: continue
            #coords=[f"{q}@{b1}->{b2}" for q,b1,b2 in recs]
            coords=[f"{q}@{c1}:{b1}->{c2}:{b2}" for q,c1,b1,c2,b2 in recs]
            f.write(f"{i1}\t{i2}\t{cnt}\t{','.join(coords)}\n")
    logging.info(f"Wrote graph edges → {fn} (threshold={read_thr})")
    toc(t0, "Building junction graph")

    graph=defaultdict(dict)
    for (i1,i2), recs in edge_support.items():
        if len(recs)>=read_thr:
            graph[i1][i2]=recs
    return graph
