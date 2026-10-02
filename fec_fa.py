#!/usr/bin/env python3

import os, logging, pysam
from fec_login import *
from fec_cycled import *
from Bio.Seq import Seq

# -----------------------
# 9. Write cycle FASTA
# -----------------------

def build_bigcycle_fasta(
    summary_bigcycle_file,
    reference_genome,
    output_dir,
    sample=None
):
    t0 = tic("Building big-cycle FASTA")

    ref = pysam.FastaFile(reference_genome)

    outfile = os.path.join(
        output_dir,
        "summary_bigcycle.fa"
    )

    with open(summary_bigcycle_file) as handle, \
         open(outfile, "w") as out:

        header = handle.readline().rstrip("\n").split("\t")

        column_index = {
            name: index
            for index, name in enumerate(header)
        }

        id_col = column_index["bigcycle_id"]
        path_col = column_index["path"]

        for line in handle:

            line = line.rstrip("\n")

            if not line:
                continue

            fields = line.split("\t")

            bigcycle_id = fields[id_col]

            path_text = fields[path_col]
            path_label = path_text.replace(" -> ", ";")

            nodes = [
                x.strip()
                for x in path_text.split("->")
                if x.strip()
            ]

            if (
                len(nodes) > 1
                and nodes[0] == nodes[-1]
            ):
                nodes = nodes[:-1]

            coords = []
            seq_parts = []

            for node in nodes:

                strand = node[-1]

                region = node[:-1]

                chrom, position_text = region.rsplit(":", 1)

                start_text, end_text = position_text.split("-", 1)

                start = int(start_text)
                end = int(end_text)

                coords.append(node)

                seq = ref.fetch(
                    chrom,
                    start,
                    end
                )

                if strand == "-":
                    seq = str(
                        Seq(seq).reverse_complement()
                    )

                seq_parts.append(seq)

            out.write(
                f">bigcycle_{bigcycle_id}_{path_label}\n"
            )

            out.write(
                "".join(seq_parts) + "\n"
            )

    ref.close()

    logging.info(
        f"Wrote big-cycle FASTA -> {outfile}"
    )

    toc(t0, "Building big-cycle FASTA")
