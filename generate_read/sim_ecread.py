import random
import os
import argparse
import math
import sys


def ensure_directory(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)


def read_genome(fasta_file, chromosome):
    sequence = []
    capture = False

    with open(fasta_file, "r") as f:
        for line in f:
            line = line.strip()

            if line.startswith(">"):
                header_chr = line[1:].split()[0]

                if capture and header_chr != chromosome:
                    break

                capture = (header_chr == chromosome)
                continue

            if capture:
                sequence.append(line)

    return "".join(sequence)


def read_bed_file(bed_file, chromosome):
    excluded_regions = []

    with open(bed_file, "r") as file:
        for line in file:
            if not line.strip() or line.startswith("#"):
                continue

            parts = line.strip().split()

            if len(parts) < 3:
                continue

            if parts[0] == chromosome:
                start = int(parts[1])
                end = int(parts[2])
                excluded_regions.append((start, end))

    return excluded_regions


def is_valid_fragment(start, length, excluded_regions,
                      used_regions, genome_length):

    end = start + length

    if end > genome_length or length < 1000:
        return False

    for ex_start, ex_end in excluded_regions + used_regions:
        if start < ex_end and end > ex_start:
            return False

    return True


def pick_random_fragments(chromosome,
                          genome_sequence,
                          num_fragments,
                          excluded_regions):

    fragments = []
    used_regions = []

    genome_length = len(genome_sequence)

    if genome_length < 1000:
        return []

    attempts = 0
    max_attempts = num_fragments * 200

    while len(fragments) < num_fragments and attempts < max_attempts:

        length = max(1000, int(random.gauss(75000, 25000)))

        if genome_length <= length:
            attempts += 1
            continue

        start = random.randint(0, genome_length - length)

        if is_valid_fragment(start,
                             length,
                             excluded_regions,
                             used_regions,
                             genome_length):

            seq = genome_sequence[start:start + length]

            fragments.append(
                (
                    chromosome,
                    start,
                    start + length - 1,
                    length,
                    seq
                )
            )

            used_regions.append(
                (
                    start,
                    start + length - 1
                )
            )

        attempts += 1

    if len(fragments) < num_fragments:
        print(
            "WARNING: Only found {} fragments for {}".format(
                len(fragments),
                chromosome
            )
        )

    return fragments


def save_to_fasta(sequence, header, directory, filename):

    path = os.path.join(directory, filename)

    with open(path, "w") as f:
        f.write(">{}\n".format(header))

        for i in range(0, len(sequence), 80):
            f.write(sequence[i:i + 80] + "\n")


def main():

    parser = argparse.ArgumentParser(
        description="Generate eccDNA sequences from random non-overlapping fragments."
    )

    parser.add_argument(
        "--genome_file",
        required=True
    )

    parser.add_argument(
        "--bed_file",
        required=True
    )

    parser.add_argument(
        "--fragments_per_chr",
        "--num_fragments_per_chromosome",
        dest="fragments_per_chr",
        type=int,
        default=5
    )

    parser.add_argument(
        "--total_eccDNA",
        type=int,
        default=10
    )

    parser.add_argument(
        "--output_directory",
        default="eccDNA_sequences"
    )

    args = parser.parse_args()

    chromosomes = ["chr{}".format(i) for i in range(1, 23)]
    chromosomes.append("chrX")

    num_chromosomes = len(chromosomes)
    max_breakpoints = 22

    required_total = args.total_eccDNA * max_breakpoints

    min_frag_per_chr = int(
        math.ceil(float(required_total) / num_chromosomes)
    )

    if args.fragments_per_chr < min_frag_per_chr:

        print(
            "Adjusting fragments_per_chr from {} to {}".format(
                args.fragments_per_chr,
                min_frag_per_chr
            )
        )

        args.fragments_per_chr = min_frag_per_chr

    ensure_directory(args.output_directory)

    all_fragments = []

    for chrom in chromosomes:

        print("Processing {}...".format(chrom))

        seq = read_genome(args.genome_file, chrom)

        if len(seq) == 0:
            print(
                "WARNING: {} not found in FASTA, skipping.".format(
                    chrom
                )
            )
            continue

        print(
            "{} length: {}".format(
                chrom,
                len(seq)
            )
        )

        excluded = read_bed_file(
            args.bed_file,
            chrom
        )

        pool = pick_random_fragments(
            chrom,
            seq,
            args.fragments_per_chr,
            excluded
        )

        all_fragments.extend(pool)

    if len(all_fragments) < 3:
        print(
            "ERROR: Insufficient fragments generated. "
            "Check FASTA and BED files."
        )
        sys.exit(1)

    coord_path = os.path.join(
        args.output_directory,
        "eccDNA_coordinates.txt"
    )

    with open(coord_path, "w") as coord_file:

        for i in range(args.total_eccDNA):

            if len(all_fragments) < 3:
                print(
                    "WARNING: Not enough fragments left "
                    "to build more eccDNA."
                )
                break

            k = random.randint(
                3,
                min(22, len(all_fragments))
            )

            ecc_frags = random.sample(
                all_fragments,
                k
            )

            for frag in ecc_frags:
                all_fragments.remove(frag)

            total_bp = sum(
                frag[3]
                for frag in ecc_frags
            )

            header = "eccDNA_B{}_{} Total_Size:{}bp".format(
                k,
                i + 1,
                total_bp
            )

            coord_file.write(
                "# {} Coordinates:\n".format(header)
            )

            for chrom, start, end, length, seq in ecc_frags:

                coord_file.write(
                    "{}\t{}\t{}\tFragment_Size:{}bp\n".format(
                        chrom,
                        start,
                        end,
                        length
                    )
                )

            coord_file.write("\n")

            full_sequence = "".join(
                frag[4]
                for frag in ecc_frags
            )

            fasta_name = "eccDNA_B{}_{}.fasta".format(
                k,
                i + 1
            )

            save_to_fasta(
                full_sequence,
                header,
                args.output_directory,
                fasta_name
            )

    print("eccDNA generation complete.")


if __name__ == "__main__":
    main()
