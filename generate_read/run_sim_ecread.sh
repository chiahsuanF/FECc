#!/bin/bash

#SBATCH --job-name=eccdna_sim
#SBATCH --output=_eccdna_sim.out
#SBATCH --error=_eccdna_sim.err
#SBATCH --partition=medium
#SBATCH --cpus-per-task=1
#SBATCH --time=20:00:00
#SBATCH --mem-per-cpu=20G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=fanc@uab.edu

# python ./sim44.py --genome_file ref.fa --num_fragments_per_chromosome [int] --total_eccDNA [int]  --bed_file hg38.largegaps.bed  --output_directory simulated_eccDNAsm  
python ./sim44.py --genome_file GRCh38.p13.genome.fa --num_fragments_per_chromosome 10 --total_eccDNA 50  --bed_file hg38.largegaps.bed  --output_directory simulated_eccDNAsm  

