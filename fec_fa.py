#!/usr/bin/env python3

import os, logging, pysam
from fec_login import *
from fec_cycled import *
from Bio.Seq import Seq

# -----------------------
# 9. Write cycle FASTA
# -----------------------


def build_cycle_fasta(filtered_cycles,
                      regions,
                      reference_fasta,
                      outdir,
                      sample):

    t0 = tic("Building cycle FASTA")

    ref = pysam.FastaFile(reference_fasta)

    outfile = os.path.join(
        outdir,
        "summary.fa"
    )

    with open(outfile, "w") as out:

        for cycle in filter_maximal_cycles(filtered_cycles):

            size = cycle_size(cycle, regions)

            coords = []

            seq_parts = []

            #
            # skip final repeated node
            #
            for node in cycle[:-1]:
                strand = node[1]
                chrom, start, end = regions[node[0]]

                coords.append(
                    f"{chrom}:{start}-{end}{strand}"
                )

                seq = ref.fetch(
                    chrom,
                    start,
                    end
                )
                if strand == "-":seq = str(Seq(seq).reverse_complement())

                seq_parts.append(seq)

            header = (
                f">{size}_"
                + ";".join(coords)
            )

            out.write(header + "\n")

            out.write("".join(seq_parts) + "\n")

    ref.close()

    logging.info(
        f"Wrote cycle FASTA → {outfile}"
    )

    toc(t0, "Building cycle FASTA")