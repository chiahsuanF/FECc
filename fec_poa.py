#!/usr/bin/env python3

import os, logging, subprocess, pysam, csv
from bisect import bisect_left, bisect_right
from concurrent.futures import ThreadPoolExecutor

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