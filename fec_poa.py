#!/usr/bin/env python3

import os, logging, subprocess, pysam
from concurrent.futures import ThreadPoolExecutor

# -----------------------
# 8. POA & Refined Connections
# -----------------------
def parse_bigcycle_node(node):
    node = node.strip()

    strand = node[-1]

    region = node[:-1]

    chrom, coords = region.rsplit(":", 1)

    start_text, end_text = coords.split("-", 1)

    return (
        chrom,
        int(start_text),
        int(end_text),
        strand
    )

def build_junctions_from_bigcycles(summary_bigcycle_file):

    junctions = []

    with open(summary_bigcycle_file, "r") as handle:

        header = handle.readline().rstrip("\n").split("\t")

        column_index = {
            name: index
            for index, name in enumerate(header)
        }

        path_column = column_index["path"]

        for line in handle:

            line = line.rstrip("\n")

            if not line:
                continue

            fields = line.split("\t")

            path_text = fields[path_column]

            nodes = [
                node.strip()
                for node in path_text.split("->")
                if node.strip()
            ]

            for i in range(len(nodes) - 1):

                c1, s1, e1, strand1 = parse_bigcycle_node(
                    nodes[i]
                )

                c2, s2, e2, strand2 = parse_bigcycle_node(
                    nodes[i + 1]
                )

                junctions.append(
                    (
                        c1,
                        s1,
                        strand1,
                        c2,
                        s2,
                        strand2
                    )
                )

    logging.info(
        "Extracted %d junction occurrences from %s",
        len(junctions),
        summary_bigcycle_file
    )

    return junctions

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
    reads = sorted(reads, key=lambda x: x.mapping_quality, reverse=True)[:50]

    return reads

def run_bsalign_poa(junction_id, reads, ref, outdir):
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

        c1, s1, strand1, c2, s2, strand2 = junction

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
                f"Junction {jid}: "
                f"{c1}:{s1}{strand1} -> "
                f"{c2}:{s2}{strand2}"
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

        except Exception:

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
