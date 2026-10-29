#!/usr/bin/env python3

import argparse, os, sys, shutil, logging, subprocess, pysam
from concurrent.futures import ThreadPoolExecutor
from fec_login import *

# -----------------------
# Consensus sequence per big cycle (run after FEC)
# Reads are cut at reference coordinates from their existing BAM
# alignments (primary + supplementary) in windows along the circle,
# and each window is collapsed with bsalign POA. Window cuts sit inside
# segments, so every junction (also the closing one) is inside a window.
# Only external tool: bsalign.
# -----------------------

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")

REF_OPS = {0, 2, 3, 7, 8}      # M D N = X consume reference
QUERY_OPS = {0, 1, 4, 7, 8}    # M I S = X consume query


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def parse_node(node):
    """
    Convert

        chr2:100-200+

    into

        ("chr2", 100, 200, "+")
    """

    node = node.strip()

    chrom, coords = node[:-1].rsplit(":", 1)
    start_text, end_text = coords.split("-", 1)

    return chrom, int(start_text), int(end_text), node[-1]


def parse_breakpoint(text):
    """
    "chrom:position" (1-based, as alignment_breakpoint writes it) ->
    (chrom, 0-based coordinate of the cut), or None for "NA".

    A left-clipped alignment starts at the breakpoint, a right-clipped
    one ends just before it, so both cut at position - 1 in 0-based,
    end-exclusive coordinates.
    """

    text = text.strip()

    if not text or text == "NA" or ":" not in text:
        return None

    chrom, position = text.rsplit(":", 1)

    try:
        return chrom, int(position) - 1
    except ValueError:
        return None


def trim_to_breakpoints(bigcycle_id, nodes, a_side, b_side, max_trim):
    """
    Move each segment's entry and exit to the exact junction breakpoints.

    Junction reads stop at the true breakpoint, which can sit a few bp
    from the region end; a read that does not reach the region end does
    not fit the window and is dropped. A breakpoint is used only when it
    lies within max_trim of its own end of the region and nearer that
    end than the other one, and the segment keeps a positive length;
    otherwise the region end is kept.

    a_side[i] is the breakpoint where node i is left, b_side[i] where
    node i + 1 is entered (the last step closes the circle).
    """

    if len(a_side) != len(nodes) or len(b_side) != len(nodes):
        return nodes

    trimmed = []

    for i, node in enumerate(nodes):

        chrom, start, end, strand = parse_node(node)

        entry = parse_breakpoint(b_side[i - 1])
        exit_ = parse_breakpoint(a_side[i])

        lo, hi = start, end

        # + is entered at its left end and left at its right end; - the
        # other way round
        left, right = (entry, exit_) if strand == "+" else (exit_, entry)

        for side, point in (("left", left), ("right", right)):

            if point is None:
                continue

            point_chrom, cut = point

            own, other = (start, end) if side == "left" else (end, start)

            if (
                point_chrom != chrom
                or abs(cut - own) > max_trim
                or abs(cut - own) > abs(cut - other)
            ):
                logging.info(
                    "bigcycle %s: breakpoint %s not used for the %s end of "
                    "%s (more than %d bp from it, or nearer the other end)",
                    bigcycle_id, f"{point_chrom}:{cut + 1}", side, node,
                    max_trim
                )
                continue

            if side == "left":
                lo = cut
            else:
                hi = cut

        if hi <= lo:
            lo, hi = start, end

        trimmed.append(f"{chrom}:{lo}-{hi}{strand}")

    return trimmed


def read_cycles(summary_bigcycle_file, ids=None, max_trim=500):
    """
    Return [(bigcycle_id, [node, ...])] in file order. The path keeps
    order and strand; a final node repeating the first is dropped. With
    the A_breakpoint and B_breakpoint columns present, each node is cut
    at its exact breakpoints (trim_to_breakpoints); max_trim 0 keeps the
    region ends.
    """

    cycles = []

    with open(summary_bigcycle_file) as handle:

        header = handle.readline().rstrip("\n").split("\t")

        id_column = header.index("bigcycle_id")
        path_column = header.index("path")
        a_column = header.index("A_breakpoint") if "A_breakpoint" in header else None
        b_column = header.index("B_breakpoint") if "B_breakpoint" in header else None

        for line in handle:

            line = line.rstrip("\n")

            if not line:
                continue

            fields = line.split("\t")

            bigcycle_id = fields[id_column]

            if ids and bigcycle_id not in ids:
                continue

            nodes = [
                node.strip()
                for node in fields[path_column].split("->")
                if node.strip()
            ]

            if len(nodes) > 1 and nodes[0] == nodes[-1]:
                nodes = nodes[:-1]

            if (
                max_trim > 0
                and a_column is not None
                and b_column is not None
                and b_column < len(fields)
            ):
                nodes = trim_to_breakpoints(
                    bigcycle_id,
                    nodes,
                    fields[a_column].split(","),
                    fields[b_column].split(","),
                    max_trim
                )

            cycles.append((bigcycle_id, nodes))

    logging.info(
        "Read %d big cycles from %s",
        len(cycles),
        summary_bigcycle_file
    )

    return cycles


# -----------------------
# Circle coordinates and windows
# -----------------------

def build_segments(nodes):
    """
    [(chrom, ref_start, ref_end, strand, t_start, t_end)] where t is the
    position on the circular template (length = sum of segment lengths).
    """

    segments = []
    t = 0

    for node in nodes:

        chrom, start, end, strand = parse_node(node)

        segments.append((chrom, start, end, strand, t, t + end - start))
        t += end - start

    return segments


def template_to_pieces(t_start, length, segments):
    """
    Reference pieces [(chrom, lo, hi, strand)] in template order for the
    circular template interval [t_start, t_start + length).
    """

    total = segments[-1][5]
    pieces = []
    pos = t_start % total
    remaining = length

    while remaining > 0:

        for chrom, start, end, strand, t0, t1 in segments:
            if t0 <= pos < t1:
                break

        take = min(remaining, t1 - pos)

        a = pos - t0
        b = a + take

        if strand == "+":
            pieces.append((chrom, start + a, start + b, strand))
        else:
            pieces.append((chrom, end - b, end - a, strand))

        remaining -= take
        pos = (pos + take) % total

    return pieces


def make_cuts(segments, window):
    """
    Window cut points on the template, only inside segments: a segment
    >= window gets cuts window/2 from each end and evenly in between
    (at most window apart); a shorter one gets a cut at its middle.
    """

    cuts = []
    half = window // 2

    for _, _, _, _, t0, t1 in segments:

        n = t1 - t0

        if n < 2:
            continue

        if n < window:
            cuts.append(t0 + n // 2)
            continue

        first = t0 + half
        last = t1 - half
        span = last - first

        k = -(-span // window)   # ceil

        if k == 0:
            cuts.append(first)
        else:
            cuts.extend(first + round(i * span / k) for i in range(k + 1))

    return sorted(set(cuts))


def make_windows(segments, window):
    """
    [(t_start, length)] between consecutive cuts, the last one wrapping
    around the end of the template.
    """

    total = segments[-1][5]
    cuts = make_cuts(segments, window)

    if not cuts:
        return [(0, total)]

    windows = []

    for i, cut in enumerate(cuts):

        nxt = cuts[i + 1] if i + 1 < len(cuts) else cuts[0] + total

        windows.append((cut, nxt - cut))

    return windows


def junction_windows(segments, windows):
    """
    [(junction label, window index)]: the window holding each junction
    (segment end -> next segment start; the last one closes the circle).
    """

    total = segments[-1][5]
    result = []

    for i, segment in enumerate(segments):

        nxt = segments[(i + 1) % len(segments)]

        label = (
            f"{segment[0]}:{segment[1]}-{segment[2]}{segment[3]}->"
            f"{nxt[0]}:{nxt[1]}-{nxt[2]}{nxt[3]}"
        )

        pos = segment[5] % total

        for w, (start, length) in enumerate(windows):

            if 0 < (pos - start) % total < length:
                result.append((label, w))
                break
        else:
            result.append((label, None))

    return result


# -----------------------
# Cutting reads at reference coordinates
# -----------------------

def query_cut(rec, ref_x):
    """
    Cut point in the forward read (hard clips included) at reference
    position ref_x, or None if ref_x is outside the alignment.
    """

    if ref_x < rec.reference_start or ref_x > rec.reference_end:
        return None

    ref = rec.reference_start
    q = 0
    found = None

    for op, length in rec.cigartuples:

        if op == 5:                     # hard clip
            q += length
            continue

        if found is None and op in REF_OPS and ref + length >= ref_x:
            found = q + (ref_x - ref) if op in QUERY_OPS else q

        if op in REF_OPS:
            ref += length

        if op in QUERY_OPS:
            q += length

    if found is None:
        found = q

    # q is now the full read length
    return q - found if rec.is_reverse else found


def forward_slice(records, lo, hi):
    """
    Forward-read bases [lo, hi) taken from any record whose stored
    sequence (hard clips excluded) covers them, or None.
    """

    for rec in records:

        if rec.query_sequence is None:
            continue

        cigar = rec.cigartuples
        left = cigar[0][1] if cigar[0][0] == 5 else 0
        right = cigar[-1][1] if cigar[-1][0] == 5 else 0

        # hard clip at the start of the forward read
        offset = right if rec.is_reverse else left

        seq = rec.get_forward_sequence()

        if offset <= lo and hi <= offset + len(seq):
            return seq[lo - offset:hi - offset]

    return None


def read_piece(records, pieces, max_gap):
    """
    The part of one read that runs through all pieces in template order,
    oriented like the template, or None if the read does not fit.
    records: all alignments of this read overlapping the pieces.
    """

    cuts = []
    orientation = None

    for chrom, lo, hi, strand in pieces:

        best = None

        for rec in records:

            if (
                rec.reference_name != chrom
                or rec.reference_start > lo
                or rec.reference_end < hi
            ):
                continue

            o = rec.is_reverse != (strand == "-")

            if orientation is not None and o != orientation:
                continue

            entry_x, exit_x = (lo, hi) if strand == "+" else (hi, lo)

            entry = query_cut(rec, entry_x)
            exit_ = query_cut(rec, exit_x)

            if entry is None or exit_ is None:
                continue

            if cuts:
                step = (entry - cuts[-1]) if not o else (cuts[-1] - entry)

                if step < 0 or step > max_gap:
                    continue
            else:
                step = 0

            if best is None or step < best[0]:
                best = (step, o, entry, exit_)

        if best is None:
            return None

        _, orientation, entry, exit_ = best
        cuts.extend([entry, exit_])

    lo, hi = min(cuts), max(cuts)

    if hi <= lo:
        return None

    piece = forward_slice(records, lo, hi)

    if piece is None:
        return None

    return revcomp(piece) if orientation else piece


# -----------------------
# One window: reads -> bsalign POA
# -----------------------

def reference_sequence(ref, pieces):

    parts = []

    for chrom, lo, hi, strand in pieces:

        seq = ref.fetch(chrom, lo, hi).upper()
        parts.append(seq if strand == "+" else revcomp(seq))

    return "".join(parts)


def run_poa(reads_fa, out_fa, mode):
    """
    bsalign POA consensus of the reads in reads_fa; returns the sequence.
    """

    cmd = [
        "bsalign", "poa",
        "-o", out_fa,
        "-m", mode,
        "-W", "64",
        "-G", "realn=1",
        reads_fa
    ]

    with open(out_fa + ".log", "w") as log:
        subprocess.check_call(cmd, stdout=log, stderr=subprocess.STDOUT)

    with pysam.FastxFile(out_fa) as fasta:
        for entry in fasta:
            return entry.sequence.upper()

    return ""


def consensus_window(index, t_start, length, segments, bam_path, ref_fasta,
                     workdir, opts):
    """
    Consensus of one window. Returns a dict for windows.tsv.
    """

    pieces = template_to_pieces(t_start, length, segments)

    by_name = {}

    with pysam.AlignmentFile(bam_path, "rb") as bam:

        for chrom, lo, hi, _ in pieces:

            for rec in bam.fetch(chrom, lo, hi):

                if (
                    rec.is_unmapped
                    or rec.is_secondary
                    or rec.is_qcfail
                    or rec.is_duplicate
                    or rec.mapping_quality < opts.min_mapq
                ):
                    continue

                by_name.setdefault(rec.query_name, {})[
                    (rec.reference_name, rec.reference_start,
                     rec.is_supplementary, rec.is_reverse)
                ] = rec

    fits = []

    for name, recs in by_name.items():

        recs = list(recs.values())

        piece = read_piece(recs, pieces, opts.max_junction_gap)

        if piece:
            mapq = min(r.mapping_quality for r in recs)
            fits.append((mapq, name, piece))

    fits.sort(key=lambda f: (-f[0], f[1]))
    used = fits[:opts.max_reads]

    result = {
        "window": index,
        "template_start": t_start,
        "template_end": t_start + length,
        "pieces": ",".join(f"{c}:{lo}-{hi}{s}" for c, lo, hi, s in pieces),
        "n_reads_fit": len(fits),
        "n_reads_used": len(used),
        "source": "reference",
        "seq": ""
    }

    if len(used) >= opts.min_reads:

        reads_fa = os.path.join(workdir, f"w{index}.fa")

        with open(reads_fa, "w") as out:
            for _, name, piece in used:
                out.write(f">{name}\n{piece}\n")

        try:
            seq = run_poa(reads_fa, os.path.join(workdir, f"w{index}.cns.fa"),
                          opts.poa_mode)

            if seq:
                result["seq"] = seq
                result["source"] = "poa"

        except Exception:
            logging.exception("window %d: bsalign poa failed", index)

    if not result["seq"]:

        with pysam.FastaFile(ref_fasta) as ref:
            result["seq"] = reference_sequence(ref, pieces)

    result["consensus_len"] = len(result["seq"])

    return result


# -----------------------
# One big cycle
# -----------------------

def consensus_cycle(bigcycle_id, nodes, bam_path, ref_fasta, outdir, opts):

    t0 = tic(f"bigcycle {bigcycle_id}: consensus")

    workdir = os.path.join(outdir, f"bigcycle_{bigcycle_id}")
    windir = os.path.join(workdir, "windows")
    os.makedirs(windir, exist_ok=True)

    segments = build_segments(nodes)
    total = segments[-1][5]
    windows = make_windows(segments, opts.window)

    logging.info(
        "bigcycle %s: %d bp, %d segments, %d windows",
        bigcycle_id, total, len(segments), len(windows)
    )

    with ThreadPoolExecutor(max_workers=opts.threads) as pool:

        results = list(pool.map(
            lambda item: consensus_window(
                item[0], item[1][0], item[1][1], segments,
                bam_path, ref_fasta, windir, opts
            ),
            enumerate(windows)
        ))

    consensus = "".join(r["seq"] for r in results)

    with open(os.path.join(workdir, "windows.tsv"), "w") as out:

        columns = ["window", "template_start", "template_end", "pieces",
                   "n_reads_fit", "n_reads_used", "source", "consensus_len"]

        out.write("\t".join(columns) + "\n")

        for r in results:
            out.write("\t".join(str(r[c]) for c in columns) + "\n")

    junctions = junction_windows(segments, windows)

    rows = [
        (label, w, results[w]["n_reads_fit"] if w is not None else 0)
        for label, w in junctions
    ]

    with open(os.path.join(workdir, "junctions.tsv"), "w") as out:

        out.write("junction\twindow\tspanning_reads\n")

        for label, w, n in rows:
            out.write(f"{label}\t{w}\t{n}\n")

    weakest = min(rows, key=lambda row: row[2])
    n_reference = sum(1 for r in results if r["source"] == "reference")

    path = " -> ".join(nodes + [nodes[0]])
    header = (
        f"bigcycle_{bigcycle_id} length={len(consensus)} circular=Y "
        f"start={nodes[0]}@{windows[0][0]} path={path}"
    )

    with open(os.path.join(workdir, "consensus.fa"), "w") as out:
        out.write(f">{header}\n{consensus}\n")

    logging.info(
        "bigcycle %s: consensus %d bp; %d/%d windows from reference; "
        "weakest junction %s (%d reads)",
        bigcycle_id, len(consensus), n_reference, len(windows),
        weakest[0], weakest[2]
    )

    toc(t0, f"bigcycle {bigcycle_id}: consensus")

    return {
        "id": bigcycle_id,
        "template_len": total,
        "consensus_len": len(consensus),
        "windows": len(windows),
        "reference_windows": n_reference,
        "min_support": weakest[2],
        "weakest_junction": weakest[0],
        "header": header,
        "seq": consensus
    }


# -----------------------
# All big cycles of one summary
# -----------------------

def run_consensus(summary_file, bam_path, ref_fasta, outdir, opts, ids=None):
    """
    Consensus of every big cycle in summary_file, written to
    consensus_bigcycle.fa and consensus_summary.tsv in outdir.
    opts needs window, max_reads, min_reads, min_mapq, max_junction_gap,
    poa_mode, threads and max_breakpoint_trim. Returns the results.
    """

    t0 = tic("Consensus sequence of big cycles")

    os.makedirs(outdir, exist_ok=True)

    cycles = read_cycles(summary_file, ids, opts.max_breakpoint_trim)

    results = []

    for bigcycle_id, nodes in cycles:

        try:
            results.append(consensus_cycle(
                bigcycle_id, nodes, bam_path, ref_fasta, outdir, opts
            ))

        except Exception:
            logging.exception("bigcycle %s: consensus failed", bigcycle_id)

    with open(os.path.join(outdir, "consensus_bigcycle.fa"), "w") as fa, \
         open(os.path.join(outdir, "consensus_summary.tsv"), "w") as tsv:

        tsv.write("bigcycle_id\ttemplate_len\tconsensus_len\twindows\treference_windows\tmin_junction_support\tweakest_junction\n")

        for r in results:
            fa.write(f">{r['header']}\n{r['seq']}\n")
            tsv.write(
                f"{r['id']}\t{r['template_len']}\t{r['consensus_len']}\t"
                f"{r['windows']}\t{r['reference_windows']}\t"
                f"{r['min_support']}\t{r['weakest_junction']}\n"
            )

    logging.info(
        "Wrote consensus for %d of %d big cycles",
        len(results), len(cycles)
    )

    toc(t0, "Consensus sequence of big cycles")

    return results


# Fixed consensus parameters (not exposed on the command line)
WINDOW = 2000               # consensus window length in bp
MAX_READS = 30              # reads per window given to bsalign poa (highest MAPQ first)
MIN_READS = 3               # fewer fitting reads -> the window keeps the reference sequence
MIN_MAPQ = 0                # ignore alignments below this MAPQ in the consensus
MAX_JUNCTION_GAP = 2000     # largest read stretch allowed between two pieces of a window
MAX_BREAKPOINT_TRIM = 500   # cut a segment at its junction breakpoint within this many bp of the region end; 0 keeps the region ends


def add_consensus_arguments(parser, prefix="consensus_"):
    """
    The consensus options, shared by this script and fec_main.py. Only
    --poa_mode is user-adjustable; the fixed parameters above are set as
    parser defaults so opts can still be the parsed arguments.
    """
    parser.add_argument("--poa_mode", type=str, default="global", help="[STR] bsalign poa -m mode (global or overlap)")
    parser.set_defaults(window=WINDOW, max_reads=MAX_READS, min_reads=MIN_READS, min_mapq=MIN_MAPQ,
        max_junction_gap=MAX_JUNCTION_GAP, max_breakpoint_trim=MAX_BREAKPOINT_TRIM)


# -----------------------
# Main
# -----------------------
if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Consensus sequence of each FEC big cycle from the BAM reads, with bsalign POA in windows along the circle")
    p.add_argument('-s', "--summary", required=True, type=str, help="[STR] summary_bigcycle.txt from FEC")
    p.add_argument('-b', "--bam", required=True, type=str, help="[STR] Sorted, indexed BAM used for FEC")
    p.add_argument('-r', "--reference_genome", required=True, type=str, help="[STR] Reference FASTA used for FEC (with .fai); used for coordinates and for windows without reads")
    p.add_argument('-o', "--output_dir", default="./consensus_out/", type=str, help="[STR] Output directory")
    p.add_argument('-t', "--threads", type=int, default=4, help="[INT] Windows run at the same time (one bsalign each)")
    p.add_argument("--ids", type=str, default="", help="[STR] Comma-separated bigcycle_id values to run (default: all)")
    add_consensus_arguments(p, prefix="")

    args = p.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    setup_logger(args.output_dir)

    if shutil.which("bsalign") is None:
        logging.error("bsalign not found on PATH")
        sys.exit(1)

    ids = {i.strip() for i in args.ids.split(",") if i.strip()}

    results = run_consensus(
        args.summary, args.bam, args.reference_genome, args.output_dir,
        args, ids
    )

    if not results:
        logging.error("No big cycles to run; exiting")
        sys.exit(1)
