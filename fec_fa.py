#!/usr/bin/env python3

import os, logging, pysam, csv
from bisect import bisect_left, bisect_right

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