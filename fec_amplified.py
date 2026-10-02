#!/usr/bin/env python3

import argparse, os, sys, random, logging,  pysam, csv, time
from bisect import bisect_left, bisect_right
import numpy as np
from fec_login import *

def merge_regions(regions):

    merged = []

    for chrom, start, end in sorted(regions):
        start = int(start)
        end = int(end)

        if start >= end:
            continue

        if (
            not merged
            or chrom != merged[-1][0]
            or start >= merged[-1][2]
        ):
            merged.append(
                [chrom, start, end]
            )
        else:
            merged[-1][2] = max(
                merged[-1][2],
                end
            )

    return [
        (chrom, start, end)
        for chrom, start, end in merged
    ]


def identify_amplified_regions(
    bam_path,
    outdir,
    sample,
    min_len=100,
    num_samples=100, max_low_depth_gap=1000
):

    t0 = tic("Identifying amplified regions")

    bam = pysam.AlignmentFile(
        bam_path,
        "rb"
    )

    allowed = {
        f"chr{i}"
        for i in range(1, 23)
    }

    allowed.update({
        "chrX",
        "chrY"
    })

    contigs = [
        (chrom, length)
        for chrom, length in zip(
            bam.references,
            bam.lengths
        )
        if chrom in allowed
    ]

    if not contigs:
        bam.close()
        logging.error(
            "No recognized chromosomes found in BAM"
        )
        sys.exit(1)
        
    covs = []

    target_valid_samples = num_samples

    max_sampling_attempts = max(
        10000,
        target_valid_samples * 100
    )

    minimum_sampling_mapq = 30
    minimum_high_quality_ratio = 0.80

    sampling_attempts = 0
    covered_positions = 0
    quality_rejected_positions = 0
    uncovered_positions = 0

    while (
        len(covs) < target_valid_samples
        and sampling_attempts < max_sampling_attempts
    ):
        sampling_attempts += 1

        chrom, chromosome_length = random.choice(
            contigs
        )

        position = random.randrange(
            chromosome_length
        )

        found_column = False

        for column in bam.pileup(
            chrom,
            position,
            position + 1,
            truncate=True,
            stepper="samtools",
            min_base_quality=0,
            max_depth=1000000
        ):
            found_column = True

            usable_reads = [
                pileup_read
                for pileup_read in column.pileups
                if not pileup_read.is_del
                and not pileup_read.is_refskip
                and not pileup_read.alignment.is_unmapped
                and not pileup_read.alignment.is_secondary
            ]

            usable_depth = len(
                usable_reads
            )

            if usable_depth <= 0:
                break

            covered_positions += 1

            high_quality_count = sum(
                pileup_read.alignment.mapping_quality
                >= minimum_sampling_mapq
                for pileup_read in usable_reads
            )

            high_quality_ratio = (
                high_quality_count
                / float(usable_depth)
            )

            if (
                high_quality_ratio
                >= minimum_high_quality_ratio
            ):
                covs.append(
                    usable_depth
                )
            else:
                quality_rejected_positions += 1

            break

        if not found_column:
            uncovered_positions += 1

    bam.close()

    logging.info(
        "Coverage sampling: attempts=%d accepted=%d/%d "
        "covered=%d uncovered=%d quality_rejected=%d "
        "minimum_MAPQ=%d minimum_quality_ratio=%.2f",
        sampling_attempts,
        len(covs),
        target_valid_samples,
        covered_positions,
        uncovered_positions,
        quality_rejected_positions,
        minimum_sampling_mapq,
        minimum_high_quality_ratio
    )

    minimum_required_samples = min(
        20,
        target_valid_samples
    )

    if len(covs) < minimum_required_samples:
        logging.error(
            "Only %d usable coverage samples were collected; "
            "at least %d are required. Check BAM coverage and "
            "mapping-quality distribution.",
            len(covs),
            minimum_required_samples
        )
        sys.exit(1)

    if len(covs) < target_valid_samples:
        logging.warning(
            "Proceeding with %d coverage samples instead of "
            "the requested %d",
            len(covs),
            target_valid_samples
        )

    coverage_array = np.asarray(
        covs,
        dtype=float
    )

    amp_thresh = (
        coverage_array.mean()
        + 3.0 * coverage_array.std()
    )

    logging.info(
        "Amplification threshold = %.1f",
        amp_thresh
    )

    # Scan the genome.
    regs = []

    bam = pysam.AlignmentFile(
        bam_path,
        "rb"
    )

    try:
        for chrom, chromosome_length in contigs:
            amplified_start = None
            previous_position = None

            for column in bam.pileup(
                chrom,
                0,
                chromosome_length,
                truncate=True
            ):
                position = column.pos

                # Close the interval if pileup skips genomic positions.
                if (
                    previous_position is not None
                    and position > previous_position + 1
                ):
                    if amplified_start is not None:
                        amplified_end = previous_position + 1

                        if (
                            amplified_end - amplified_start
                            >= min_len
                        ):
                            regs.append(
                                (
                                    chrom,
                                    amplified_start,
                                    amplified_end
                                )
                            )

                    amplified_start = None

                if column.n >= amp_thresh:
                    if amplified_start is None:
                        amplified_start = position
                else:
                    if amplified_start is not None:
                        # Half-open endpoint. Do not subtract one.
                        amplified_end = position

                        if (
                            amplified_end - amplified_start
                            >= min_len
                        ):
                            regs.append(
                                (
                                    chrom,
                                    amplified_start,
                                    amplified_end
                                )
                            )

                    amplified_start = None

                previous_position = position

            # Close an amplified interval at the final covered position.
            if amplified_start is not None:
                amplified_end = (
                    previous_position + 1
                    if previous_position is not None
                    else chromosome_length
                )

                if (
                    amplified_end - amplified_start
                    >= min_len
                ):
                    regs.append(
                        (
                            chrom,
                            amplified_start,
                            amplified_end
                        )
                    )

    finally:
        bam.close()

    output_file = os.path.join(
        outdir,
        f"{sample}_amplified_regions.txt"
    )

    with open(output_file, "w") as handle:
        handle.write(
            "chrom\tstart\tend\n"
        )

        for chrom, start, end in regs:
            handle.write(
                f"{chrom}\t{start}\t{end}\n"
            )

    logging.info(
        "Wrote %d amplified regions -> %s",
        len(regs),
        output_file
    )

    toc(
        t0,
        "Identifying amplified regions"
    )

    return regs, amp_thresh
