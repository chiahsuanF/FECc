#!/usr/bin/env python3

import argparse, os, sys, time,random,logging, subprocess, pysam, csv
from collections import defaultdict
from bisect import bisect_left, bisect_right
from concurrent.futures import ThreadPoolExecutor
import numpy as np
from fec_login import *

# -----------------------
# Logger setup
# -----------------------
def setup_logger(output_dir):
    logfile = os.path.join(output_dir, "pipeline.log")
    logging.basicConfig(
        filename=logfile,
        filemode="w",
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s"
    )
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s"))
    logging.getLogger().addHandler(ch)
    logging.info("Logger initialized")

def tic(msg):
    logging.info(f"START: {msg}")
    return time.time()

def toc(start, msg):
    logging.info(f" END: {msg} ({time.time()-start:.1f}s)")

# -----------------------
# merge overlapping regions
# -----------------------
def merge_regions(regions):
    merged = []
    for chrom, start, end in sorted(regions):
        if not merged or chrom != merged[-1][0] or start > merged[-1][2]:
            merged.append([chrom, start, end])
        else:
            merged[-1][2] = max(merged[-1][2], end)
    return [(c, s, e) for c, s, e in merged]

# -----------------------
# 1. Depth & MinSupport
# -----------------------
"""
def estimate_depth_and_minsupp(bam_path, min_support, depth_override,
                               outdir, sample):
    t0 = tic("Estimating depth & minsupp")
    if min_support is None:
        if depth_override is not None:
            depth = depth_override
            minsupp = round(depth/10) + 2
        else:
            total_bases = 0
            ref_len = sum(pysam.AlignmentFile(bam_path,'rb').lengths)
            with pysam.AlignmentFile(bam_path,'rb') as bam:
                for r in bam.fetch(until_eof=True):
                    if not r.is_unmapped:
                        total_bases += r.query_length
            depth = total_bases / ref_len
            minsupp = round(depth/10) + 2
    else:
        minsupp = min_support
        if depth_override is None:
            total_bases = 0
            ref_len = sum(pysam.AlignmentFile(bam_path,'rb').lengths)
            with pysam.AlignmentFile(bam_path,'rb') as bam:
                for r in bam.fetch(until_eof=True):
                    if not r.is_unmapped:
                        total_bases += r.query_length
            depth = total_bases / ref_len
        else:
            depth = depth_override

    fn = os.path.join(outdir, f"01_{sample}_depth_minsupp.txt")
    with open(fn,'w') as f:
        f.write(f"Estimated Depth: {depth:.2f}\n")
        f.write(f"MinSupport: {minsupp}\n")
    logging.info(f"Wrote depth & minsupp → {fn}")
    toc(t0, "Estimating depth & minsupp")
    return depth, minsupp
"""
# -----------------------------
# worker: compute depth for one contig
# -----------------------------
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
    fn = os.path.join(outdir, f"01_{sample}_depth_minsupp.txt")

    with open(fn, "w") as f:
        f.write(f"Estimated Depth: {depth:.2f}\n")
        f.write(f"MinSupport: {minsupp}\n")

    logging.info(f"Wrote depth & minsupp → {fn}")

    return depth, minsupp


# -----------------------
# 2. Identify Amplified Regions // find hight coverage region(correct)
# -----------------------
def identify_amplified_regions(bam_path, outdir, sample,
                               min_len=1000, num_samples=1000):
    t0 = tic("Identifying amplified regions")
    bam = pysam.AlignmentFile(bam_path,'rb')
    contigs = list(zip(bam.references, bam.lengths))

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

    fn = os.path.join(outdir, f"02_{sample}_amplified_regions.txt")
    with open(fn,'w') as f:
        f.write("chrom\tstart\tend\n")
        for c,s,e in regs:
            f.write(f"{c}\t{s}\t{e}\n")
    logging.info(f"Wrote {len(regs)} amplified regions → {fn}")
    toc(t0, "Identifying amplified regions")
    return regs, amp_thresh

# -----------------------
# 3. Extract Breakpoints with MapQ  // c1, bp1 connect c2, bp2 
# -----------------------
def extract_breakpoints(bam_path, outdir, sample):
    t0 = tic("Collecting split-read breakpoints")
    bam = pysam.AlignmentFile(bam_path,'rb')
    bps = []  # (qname, c1, bp1, mq1, c2, bp2, mq2)
    for r in bam.fetch(until_eof=True):
        if not r.has_tag('SA'): continue
        blocks = r.get_blocks()
        if not blocks: continue
        bp1 = blocks[-1][1]; mq1 = r.mapping_quality; c1 = r.reference_name
        qname = r.query_name
        for sa in r.get_tag('SA').strip(';').split(';'):
            if not sa: continue
            c2, bp2, *_ = sa.split(',')[:3]
            mq2 = int(sa.split(',')[4])
            bps.append((qname, c1, bp1, mq1, c2, int(bp2), mq2))
    bam.close()

    fn = os.path.join(outdir, f"03_{sample}_breakpoints.txt")
    with open(fn,'w') as f:
        f.write("qname\tc1\tbp1\tmq1\tc2\tbp2\tmq2\n")
        for q,c1,b1,m1,c2,b2,m2 in bps:
            f.write(f"{q}\t{c1}\t{b1}\t{m1}\t{c2}\t{b2}\t{m2}\n")
    logging.info(f"Wrote {len(bps)} breakpoints → {fn}")
    toc(t0, "Collecting split-read breakpoints")
    return bps

# -----------------------
# 4. Refine Regions by High-Q Clusters // 
# -----------------------
def refine_regions_by_highQ_breakpoints(regions, bps, outdir, sample, window=1000):
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
            if s<=pos<=e: return (s,e)
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
            if s<=pos<=e: return i
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

# -----------------------
# 6. Cycle Detection
# -----------------------
def dfs_visit(node, graph, path, visited, cycles):
    if node in path:
        idx=path.index(node)
        cycles.append(path[idx:]+[node])
        return
    if node in visited:
        return
    visited.add(node)
    for nbr in graph.get(node,{}):
        dfs_visit(nbr, graph, path+[node], visited, cycles)

def find_cycles(graph, outdir, sample):
    t0 = tic("Detecting cycles")
    visited, cycles = set(), []
    for n in graph:
        if n not in visited:
            dfs_visit(n, graph, [], visited, cycles)
    fn = os.path.join(outdir, f"06_{sample}_raw_cycles.txt")
    with open(fn,'w') as f:
        f.write("cycle_nodes\n")
        for c in cycles:
            f.write(",".join(map(str,c))+"\n")
    logging.info(f"Wrote {len(cycles)} raw cycles → {fn}")
    toc(t0, "Detecting cycles")
    return cycles

# -----------------------
# Write cycles with real coords
# -----------------------
def cycle_size(cycle, regions):
    """Return total bp size of a cycle."""
    size = 0
    # ignore last node because it repeats the first node
    for idx in cycle[:-1]:
        chrom, start, end = regions[idx]
        size += (end - start + 1)
    return size


def filter_maximal_cycles(cycles):
    """
    Remove cycles that are strict subsets of another cycle.
    """
    keep = []

    node_sets = [set(c[:-1]) for c in cycles]

    for i, s1 in enumerate(node_sets):
        is_subset = False

        for j, s2 in enumerate(node_sets):
            if i == j:
                continue

            if s1 < s2:      # strict subset
                is_subset = True
                break

        if not is_subset:
            keep.append(cycles[i])

    return keep

def write_cycles_with_coords(cycles, regions, outdir, sample, label):
    cycles = filter_maximal_cycles(cycles) #
    fn = os.path.join( outdir, f"07_{sample}_{label}_cycles_coords.txt")

    with open(fn, 'w') as f:
        f.write("cycle_id\tsize\tpath\n")
        for cid, cycle in enumerate(cycles):
            size = cycle_size(cycle, regions)
            path = ";".join(
                f"{regions[i][0]}:{regions[i][1]}-{regions[i][2]}"
                for i in cycle
            )
            f.write( f"{cid}\t{size}\t{path}\n")

    logging.info(
        f"Wrote {len(cycles)} {label} cycles with coordinates → {fn}"
    )

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

def filter_and_write_cycles(cycles, regs, bam_path, outdir, sample):
    t0 = tic("Filtering cycles")
    filtered=[c for c in cycles if all(check_region_quality(regs[i], bam_path) for i in c)]
    fn = os.path.join(outdir, f"08_{sample}_filtered_cycles.txt")
    with open(fn,'w') as f:
        f.write("cycle_nodes\n")
        for c in filtered:
            f.write(",".join(map(str,c))+"\n")
    logging.info(f"Wrote {len(filtered)} filtered cycles → {fn}")
    toc(t0, "Filtering cycles")
    return filtered

# -----------------------
# 8. POA & Refined Connections
# -----------------------
def build_unique_junctions(cycles, regs):
    """
    Convert cycles → unique junctions (deduplicates POA calls)
    """
    junctions = set()

    for cyc in cycles:
        coords = [regs[i] for i in cyc]

        for j in range(len(coords) - 1):
            c1, s1, _ = coords[j]
            c2, s2, _ = coords[j + 1]

            # canonical key (direction-aware)
            key = (c1, s1, c2, s2)
            junctions.add(key)

    return list(junctions)

def fetch_reads_cached(bam_path, chrom, pos, window=150):
    bam = pysam.AlignmentFile(bam_path, "rb")

    lo = max(0, pos - window)
    hi = pos + window

    reads = []

    for r in bam.fetch(chrom, lo, hi):

        if r.is_unmapped:
            continue

        if r.query_sequence is None:
            continue

        # keep only informative reads
        if r.mapping_quality < 10:
            continue

        reads.append(r)

    bam.close()

    # keep only top reads (VERY IMPORTANT for speed)
    reads = sorted(reads, key=lambda x: x.mapping_quality, reverse=True)[:50]

    return reads

def run_bsalign_poa(junction_id, reads, ref, outdir):
    """
    Run BSAlign POA for one junction
    """

    prefix = os.path.join(outdir, f"poa_{junction_id}")

    fasta = prefix + "_reads.fa"

    with open(fasta, "w") as f:
        for i, r in enumerate(reads):
            f.write(f">{r.query_name}_{i}\n")
            f.write(r.query_sequence + "\n")

    cons = prefix + "_consensus.fa"

    cmd = [ "bsalign", "poa", "-o", cons, "-m", "overlap",
            "-W", "64", "-G", "realn=1", fasta]

    subprocess.check_call(cmd)

    return cons

def process_junctions_fast(junctions, bam_path, ref, outdir, threads=8):

    os.makedirs(outdir, exist_ok=True)

    logging.info(f"Running POA on {len(junctions)} junctions")

    def worker(item):

        jid, junction = item

        c1, s1, c2, s2 = junction

        logging.info(
            f"Junction {jid}: {c1}:{s1} -> {c2}:{s2}"
        )

        reads1 = fetch_reads_cached( bam_path, c1, s1)
        

        reads2 = fetch_reads_cached( bam_path, c2, s2)

        reads = reads1 + reads2

        logging.info(
            f"Junction {jid}: reads1={len(reads1)} "
            f"reads2={len(reads2)} "
            f"total={len(reads)}"
        )

        if len(reads) < 5:
            logging.info(
                f"Junction {jid}: skipped (<5 reads)"
            )
            return None

        try:

            cons = run_bsalign_poa(
                jid,
                reads,
                ref,
                outdir
            )

            logging.info(
                f"Junction {jid}: consensus={cons}"
            )

            return (jid, cons)

        except Exception as e:

            logging.exception(
                f"Junction {jid}: POA failed"
            )

            return None

    results = []

    with ThreadPoolExecutor(max_workers=threads) as pool:

        for res in pool.map(
            worker,
            enumerate(junctions)
        ):
            if res is not None:
                results.append(res)

    logging.info(
        f"Generated {len(results)} consensus sequences"
    )

    return results

def extract_reads_for_poa(bam_path, chrom, pos, window=100, min_sc=20):
    bam=pysam.AlignmentFile(bam_path,'rb')
    reads=[]
    lo,hi=max(0,pos-window),pos+window
    for r in bam.fetch(chrom, lo, hi):
        if r.is_unmapped: continue
        if any(bs<=pos<be for bs,be in r.get_blocks()):
            reads.append(r)
        else:
            sc=sum(l for op,l in r.cigartuples if op==4)
            if sc>=min_sc:
                reads.append(r)
    bam.close()
    return reads
"""
def map_consensus_sequence(seq_id, seq_lines, ref, prefix):

    fa = f"{prefix}_{seq_id}.fa".replace(" ", "_")

    with open(fa, "w") as f:
        f.write(f">{seq_id}\n")
        f.write("".join(seq_lines) + "\n")

    outbam = fa.replace(".fa", "_sorted.bam")

    cmd1 = [
        "minimap2",
        "-t", "4",
        "-ax", "asm5",
        ref,
        fa
    ]

    cmd2 = [
        "samtools",
        "sort",
        "-@",
        "4",
        "-o",
        outbam
    ]

    cmd3 = [
        "samtools",
        "index",
        outbam
    ]

    with subprocess.Popen(cmd1, stdout=subprocess.PIPE) as p1:
        with subprocess.Popen(cmd2, stdin=p1.stdout) as p2:
            p1.stdout.close()
            p2.communicate()

    subprocess.check_call(cmd3)

    return outbam
"""

def map_consensus_sequence(seq_id, seq_lines, ref, prefix):

    fa = f"{prefix}_{seq_id}.fa".replace(" ", "_")

    # -----------------------
    # write consensus FASTA
    # -----------------------
    with open(fa, "w") as f:
        f.write(f">{seq_id}\n")
        f.write("".join(seq_lines) + "\n")

    samfile = fa.replace(".fa", ".sam")
    unsorted_bam = fa.replace(".fa", ".bam")
    outbam = fa.replace(".fa", "_sorted.bam")

    # -----------------------
    # minimap2 alignment
    # -----------------------
    cmd = ["minimap2", "-t", "4", "-ax", "asm5", ref, fa]

    with open(samfile, "w") as sam_out:
        subprocess.check_call(
            cmd,
            stdout=sam_out
        )

    # -----------------------
    # SAM -> BAM
    # -----------------------
    pysam.view("-bS", "-o", unsorted_bam, samfile, catch_stdout=False)

    # -----------------------
    # sort BAM
    # -----------------------
    pysam.sort("-@", "4", "-o", outbam, unsorted_bam)

    # -----------------------
    # index BAM
    # -----------------------
    pysam.index(outbam)

    # -----------------------
    # cleanup
    # -----------------------
    for tmp in (samfile, unsorted_bam):
        if os.path.exists(tmp):
            os.remove(tmp)

    return outbam



def poa_and_mapping(prefix, ref, reads):

    t0 = tic(f"BSAlign POA {prefix}")

    fasta = prefix + "_reads.fasta"

    with open(fasta, "w") as f:
        for i, r in enumerate(reads):

            if len(r.query_sequence) < 50:
                continue

            f.write(f">{r.query_name}_{i}\n")
            f.write(r.query_sequence + "\n")

    cons = prefix + "_consensus.fa"

    cmd = [
        "bsalign",
        "poa",
        "-o", cons,
        "-m", "overlap",
        "-G", "realn=2,trigger=1",
        fasta
    ]

    logging.info("Running: %s", " ".join(cmd))

    subprocess.check_call(cmd)

    mapped = []

    sid = None
    seq = []

    with open(cons) as f:

        for line in f:

            if line.startswith(">"):

                if sid is not None:

                    mapped.append(
                        map_consensus_sequence(
                            sid,
                            seq,
                            ref,
                            prefix
                        )
                    )

                sid = line[1:].strip()
                seq = []

            else:
                seq.append(line.strip())

    if sid is not None:

        mapped.append(
            map_consensus_sequence(
                sid,
                seq,
                ref,
                prefix
            )
        )

    toc(t0, f"BSAlign POA {prefix}")

    return mapped


def extract_refined_coords(bam_path, prefix):
    bam=pysam.AlignmentFile(bam_path,'rb')
    segs=[r for r in bam.fetch(until_eof=True) if not r.is_unmapped and not r.is_secondary]
    bam.close()
    segs.sort(key=lambda r:r.query_alignment_start)
    out=[]
    for a,b in zip(segs, segs[1:]):
        ba=a.get_blocks()[-1][1]; bb=b.get_blocks()[0][0]
        sa='-' if a.is_reverse else '+'; sb='-' if b.is_reverse else '+'
        out.append((a.reference_name,ba,sa,b.reference_name,bb,sb))
    sup=sum(1 for _ in open(prefix+"_reads.fasta") if _.startswith(">"))
    return out, sup

def process_all_connections(bam_path, ref, regs, cycles, outdir, sample):
    t0 = tic("Processing POA connections")
    ws=os.path.join(outdir,"poa_workspace"); os.makedirs(ws,exist_ok=True)
    fn=os.path.join(outdir,f"09_{sample}_refined_connections.txt")
    with open(fn,'w') as out:
        out.write("junction\treads_support\n")
        for cid,cyc in enumerate(cycles):
            coords=[regs[i] for i in cyc]
            for j in range(len(coords)-1):
                c1,s1,_=coords[j]; c2,s2,_=coords[j+1]
                reads=extract_reads_for_poa(bam_path,c1,s1)+extract_reads_for_poa(bam_path,c2,s2)
                if not reads: continue
                prefix=os.path.join(ws,f"cycle{cid}_edge{j}")
                for mb in poa_and_mapping(prefix, ref, reads):
                    jcts,sup=extract_refined_coords(mb,prefix)
                    for a,ea,sa,b,eb,sb in jcts:
                        out.write(f"{a}:{ea}({sa})->{b}:{eb}({sb})\t{sup}\n")
    logging.info(f"Wrote refined connections → {fn}")
    toc(t0, "Processing POA connections")

# -----------------------
# 9. Write cycle FASTA
# -----------------------

def wrap_seq(seq, width=80):
    return "\n".join(
        seq[i:i+width]
        for i in range(0, len(seq), width)
    )


def build_cycle_fasta(filtered_cycles,
                      regions,
                      reference_fasta,
                      outdir,
                      sample):

    t0 = tic("Building cycle FASTA")

    ref = pysam.FastaFile(reference_fasta)

    outfile = os.path.join(
        outdir,
        f"summary.fa"
    )

    with open(outfile, "w") as out:

        for cycle in filter_maximal_cycles(filtered_cycles):

            size = cycle_size(cycle, regions)

            coords = []

            seq_parts = []

            #
            # skip final repeated node
            #
            for idx in cycle[:-1]:

                chrom, start, end = regions[idx]

                coords.append(
                    f"{chrom}:{start}-{end}"
                )

                seq = ref.fetch(
                    chrom,
                    start,
                    end
                )

                seq_parts.append(seq)

            header = (
                f">{size}_"
                + ";".join(coords)
            )

            out.write(header + "\n")

            out.write(
                wrap_seq(
                    "".join(seq_parts)
                )
                + "\n"
            )

    ref.close()

    logging.info(
        f"Wrote cycle FASTA → {outfile}"
    )

    toc(t0, "Building cycle FASTA")

# -----------------------
# Main
# -----------------------
if __name__=="__main__":
    p=argparse.ArgumentParser(description="Amplified-region → graph → cycles → POA pipeline")
    p.add_argument('-b',"--bam", required=True, type=str, help='[STR] Input sorted bam, index required')
    p.add_argument('-o',"--output_dir", required=True, type=str, help='[STR] Name of output directory')
    p.add_argument('-r',"--reference_genome", required=True, type=str, help='[STR] Import reference chromosome genome')
    p.add_argument('-ms',"--min_support", type=int, help='[INT] Minimum number of support reads for insertion reads')
    p.add_argument('-d',"--depth", type=float, help='[FLOAT]Approximately total depth')
    p.add_argument('-t',"--threads", type=int, default=4, help="[INT] Number of CPU threads to use")

    args=p.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    setup_logger(args.output_dir)

    sample=os.path.splitext(os.path.basename(args.bam))[0]
    logging.info(f"Starting pipeline for sample {sample}")

    # 1) Depth & MinSupport
    depth, minsupp = estimate_depth_and_minsupp( args.bam, args.min_support, args.depth, args.output_dir, sample)
    
    # 2) Coarse amplified regions
    regs, amp_thresh = identify_amplified_regions( args.bam, args.output_dir, sample)
    if not regs:
        logging.error("No amplified regions found; exiting"); sys.exit(1)
        
    # 2b) Merge overlapping regions
    merged = merge_regions(regs)
    fnm = os.path.join(args.output_dir, f"10_{sample}_merged_regions.txt")
    with open(fnm,'w') as f:
        f.write("chrom\tstart\tend\n")
        for c,s,e in merged:
            f.write(f"{c}\t{s}\t{e}\n")
    logging.info(f"Wrote {len(merged)} merged regions → {fnm}")

    # 3) Breakpoints
    bps = extract_breakpoints( args.bam, args.output_dir, sample)
    
    # 4) Refine on merged regions
    refined = refine_regions_by_highQ_breakpoints( merged, bps, args.output_dir, sample, window=100000) #1000

    # 5) Build Junction Graph
    graph = build_graph(refined, bps, amp_thresh, args.output_dir, sample)
    
    # 6) Detect cycles & write raw coords
    raw_cycles = find_cycles( graph, args.output_dir, sample)
    write_cycles_with_coords(raw_cycles, refined, args.output_dir, sample, "raw")
    if not raw_cycles:
        logging.error("No cycles detected; exiting"); sys.exit(1)
        
    # 7) Filter & write filtered coords
    filtered = filter_and_write_cycles( raw_cycles, refined, args.bam, args.output_dir, sample)
    write_cycles_with_coords(filtered, refined, args.output_dir, sample, "filtered")
    if not filtered:
        logging.error("No cycles pass quality filter; exiting"); sys.exit(1)
        
    # 8) POA & refined connections
    junctions = build_unique_junctions(filtered, refined)

    results = process_junctions_fast(
        junctions,
        args.bam,
        args.reference_genome,
        os.path.join(args.output_dir, "poa"),
        threads=args.threads)
    
    build_cycle_fasta( filtered, refined, args.reference_genome, args.output_dir, sample)
    
    logging.info("Pipeline completed successfully")

