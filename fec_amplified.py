#!/usr/bin/env python3

import argparse, os, sys, random, logging,  pysam, csv, time
import multiprocessing as mp
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


SCAN_CHUNK_SIZE = 10_000_000


def _scan_chunk(job):

    bam_path, chrom, start, end, amp_thresh = job

    runs = []
    amplified_start = None
    previous_position = None

    with pysam.AlignmentFile(bam_path, "rb") as bam:

        for column in bam.pileup(
            chrom,
            start,
            end,
            truncate=True
        ):
            position = column.pos

            if (
                previous_position is not None
                and position > previous_position + 1
            ):
                if amplified_start is not None:
                    runs.append(
                        (amplified_start, previous_position + 1)
                    )

                amplified_start = None

            if column.n >= amp_thresh:
                if amplified_start is None:
                    amplified_start = position
            else:
                if amplified_start is not None:
                    runs.append(
                        (amplified_start, position)
                    )

                amplified_start = None

            previous_position = position

    if amplified_start is not None:
        runs.append(
            (amplified_start, previous_position + 1)
        )

    return runs


def identify_amplified_regions(
    bam_path,
    outdir,
    sample,
    min_len=100,
    num_samples=1000, max_low_depth_gap=0, seed=1,
    threads=1
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

    rng = random.Random(seed)

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

        chrom, chromosome_length = rng.choice(
            contigs
        )

        position = rng.randrange(
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

    median_depth = float(
        np.median(coverage_array)
    )

    mad_depth = 1.4826 * float(
        np.median(np.abs(coverage_array - median_depth))
    )

    amp_thresh = max(
        median_depth + 3.0 * mad_depth,
        2.0 * median_depth
    )

    logging.info(
        "Amplification threshold = %.1f (median=%.1f, MAD=%.2f, n=%d; "
        "max(median+3*MAD, 2*median))",
        amp_thresh,
        median_depth,
        mad_depth,
        len(coverage_array)
    )

    jobs = [
        (bam_path, chrom, start, min(start + SCAN_CHUNK_SIZE, chromosome_length), amp_thresh)
        for chrom, chromosome_length in contigs
        for start in range(0, chromosome_length, SCAN_CHUNK_SIZE)
    ]

    logging.info(
        "Scanning %d windows with %d process(es)",
        len(jobs),
        threads
    )

    if threads > 1:
        with mp.get_context("fork").Pool(threads) as pool:
            chunk_runs = pool.map(_scan_chunk, jobs, chunksize=1)
    else:
        chunk_runs = [_scan_chunk(job) for job in jobs]

    regs = []
    gaps = []
    open_run = None

    for (_, chrom, start, _, _), runs in zip(jobs, chunk_runs):

        for run_start, run_end in runs:

            if (
                open_run is not None
                and open_run[0] == chrom
                and run_start - open_run[2] <= max_low_depth_gap
            ):
                if run_start > open_run[2]:
                    gaps.append((chrom, open_run[2], run_start))
                open_run[2] = run_end
                continue

            if open_run is not None and open_run[2] - open_run[1] >= min_len:
                regs.append(tuple(open_run))

            open_run = [chrom, run_start, run_end]

    if open_run is not None and open_run[2] - open_run[1] >= min_len:
        regs.append(tuple(open_run))

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

    return regs, amp_thresh, gaps
