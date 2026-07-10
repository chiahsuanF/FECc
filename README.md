# FEC

Finding ECCDNA from long-read sequencing data.

Author: Chia Hsuan(Josh) Fan

Email: rge4br@virginia.edu

## Quick Start

```sh
git clone https://github.com/chiahsuanF/fec.git
cd fec/
```
```
./fec_main.py -h
```

```sh
# Finding NUMTs with indexed bam. The reference using human hg38.
./fec_main.py -b <index.bam> -o <output_dir> -r <ref.fa>

# Example usage:
./fec_main.py -b index_aligned_human_hg38.bam -o fec_output -r hg38.fa
```

```sh
# Finding NUMTs with indexed bam. The reference using other species' chrM indexed by BWA.
./fec_main.py -b <index.bam> -o <output_dir> -r <ref.fa> -ms <int> -d <flost> -t <int>

# Example usage:
./fec_main.py -b index_aligned_human_hg38.bam -o fec_output -r hg38.fa -m 15 -d 30 -t 4
```

## Description
fec is a tool for finding ECCDNA(Extrachromosomal circular DNA) in nuclear genome.      

Amplified-region → graph → cycles → POA pipeline

The input file should be a long-read indexed BAM(PacBio Hifi and Oxford Nanopore) file. The output is a standard VCF file which passed filter. This program was test on a x86_64 Linux system with 40GB physical memory.

## Depencency
Dependencies for fec:

* python3(version 3.10.9 or later)
* pysam(version 0.22.0 or later)
* argparse(version 1.4.0 or later)
* pandas(version 2.1.1 or later)
* minimap2(version 2.24-r1122 or later)
* bsalign(version 1.2.1 or later)
* numpy(version 1.26.0 or later)

  
## Installation

### Git
```
git clone https://github.com/chiahsuanF/fec.git
```
Next, add the directory to PATH:
```
export PATH=$PWD/fec/:$PATH
```

## General usage
| Option |  |Description| Required|
|:---------|:-------------|:-------------|:-------------|
|-h |--help|Show this help message and exit||
|-b |--bam|[[STR] Input sorted bam, index required |True|
|-o |--output_dir|[STR] Name of output directory|True|
|-r |--reference_genome|[STR] Import reference chromosome genome|True|
|-ms |--min_support|[INT] Minimum number of support reads for insertion reads|Fales|
|-d |--depth|[FLOAT]Approximately total depth|Fales|
|-t |--threads|[INT] Number of CPU threads to use|Fales|

## Use cases
fec required input as an aligned BAM file with header. Mapping raw data suggests using minimap2 and samtools. Input reference only needs reference chrM.fa, duplicated sequence and required indexed by using BWA.

### Mapping long-read raw data with minimap2 and samtools:
Download the reference FASTA file and sample's FASTQ file, then:
```
minimap2 -ax map-hifi ref.fa sample.fq > aligned.sam
samtools sort aligned.sam -o aligned.bam
samtools index -b aligned.bam
```
after this process, user will get <aligned.bam> file and <aligned.bam.bai> file.    

For more detail, please look at [Minimap2 GitHub pages](https://github.com/lh3/minimap2) and [SAMtools GitHub pages](https://github.com/samtools/samtools).

## Output files
The output file directory includes:

### summary.txt
| cycle_id |	size	| path|
|:---------|:-------------|:-------------|
|ID|eccDNA size|eccDNA path|
```
cycle_id	size	path
0	94077	chr22:24446938-24518142;chr14:42123483-42146354;chr22:24446938-24518142
1	111332	chr2:15829592-15854254;chr6:116974394-117061062;chr2:15829592-15854254
2	312359	chr14:40103961-40230096;chr21:15051654-15079079;chr22:31441698-31504777;chr18:15179505-15275221;chr14:40103961-40230096
```

### summary.fa
 &gt; size_path <br> seq

```
>94077_chr22:24446938-24518142;chr14:42123483-42146354;chr22:24446938-24518142
ACTGGGTAATTTATAAAGGAGAGAGGTTTAATTACTCACAGTTACACATGGCTGGAGAGGCCTCACAATCATGATGGAAG...
>111332_chr2:15829592-15854254;chr6:116974394-117061062;chr2:15829592-15854254
TGTCTTCTCTGTTTGTGAACAACCCTGGTATCTCTGTGTATCCAAATTTCTCTTTTTATAATGACTTGGTCAGATTGAAT...
>312359_chr14:40103961-40230096;chr21:15051654-15079079;chr22:31441698-31504777;chr18:15179505-15275221;chr14:40103961-40230096
CTCTGACTATATGAGAGAAGAATTTGGACAAAAGATAATAAAGAAATCTTTGAGATAAGTTACCTGAAACTGCCAAGAAA...
```


### login_not.txt
```
2026-06-15 09:56:05,095 [INFO] Logger initialized
2026-06-15 09:56:05,096 [INFO] Starting pipeline for sample index_aligned_human_hg38
2026-06-15 09:56:05,096 [INFO] Estimating depth & minsupp (multi-threaded)
...
2026-06-15 10:24:11,809 [INFO] Wrote cycle FASTA → 08_L5_v4_v2_output/summary.fa
2026-06-15 10:24:11,809 [INFO]  END: Building cycle FASTA (1.6s)
2026-06-15 10:24:11,809 [INFO] Pipeline completed successfully
```

## Additional information
1. [pysam](https://pysam.readthedocs.io/en/latest/api.html)
2. [argparse](https://docs.python.org/3/library/argparse.html)
3. [pandas](https://devdocs.io/pandas~0.25/)
4. [minimap2](https://github.com/lh3/minimap2)
5. [bsalign](https://github.com/ruanjue/bsalign)


