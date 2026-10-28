
# FEC
test
Finding ECCDNA(FEC) from long-read sequencing data.

Author: Chia Hsuan(Josh) Fan

Email: rge4br@virginia.edu

## Quick Start

```sh
git clone https://github.com/chiahsuanF/FECc.git
cd FECc/
```
```
./fec_main.py -h
```

```sh
# Finding ECCDNA with indexed bam using required functions. The reference using human hg38.
./fec_main.py -b <index.bam> -r <ref.fa>

# Example usage:
./fec_main.py -b index_aligned_human_hg38.bam -r hg38.fa
```

```sh
# Finding ECCDNA with indexed bam using full function.
./fec_main.py -b <index.bam> -o <output_dir> -r <ref.fa> -ms <int> -d <flost> -t <int>

# Example usage:
./fec_main.py -b index_aligned_human_hg38.bam -o fec_output -r hg38.fa -ms 15 -d 30 -t 4
```

## Description
fec is a tool for finding ECCDNA(Extrachromosomal circular DNA) in nuclear genome.      

Amplified-region → graph → cycles → POA pipeline

The input file should be a long-read indexed BAM(PacBio Hifi and Oxford Nanopore) file. The output is a standard VCF file which passed filter. This program was test on a x86_64 Linux system with 40GB physical memory.

## Dependency
Dependencies for fec:

* python3(version 3.10 or later)
* pysam(version 0.22.0 or later)
* numpy(version 1.23 or later)
* scipy(version 1.9 or later, for `scipy.optimize.milp`)
* matplotlib(version 3.5 or later)
* bsalign(version 1.2.1 or later), must be on PATH; only needed for the consensus step (skip it with `--no_consensus`)

minimap2 and samtools are only needed to prepare the input BAM (see [Use cases](#use-cases)); fec does not call them.

## Installation

### Git
```
git clone https://github.com/chiahsuanF/FECc.git
cd FECc/
pip install -r requirements.txt
```

### bsalign
```
git clone https://github.com/ruanjue/bsalign.git
cd bsalign && make
export PATH=$PWD:$PATH
```

Next, add the fec directory to PATH:
```
export PATH=/path/to/FECc/:$PATH
```

## General usage
| Option |  |Description| Required|
|:---------|:-------------|:-------------|:-------------|
|-h |--help|Show this help message and exit||
|-b |--bam|[STR] Input sorted bam, index required |True|
|-o |--output_dir|[STR] Name of output directory. Default: fec_output|Fales|
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

### summary_bigcycle.txt
| column name |		| 
|:---------|:-------------|
|bigcycle_id	| the circle's number|
|total_size	| the lengths of all the pieces|
|fragment_count | number of fregment make up the circle |	
|region_depths	| number of coverage that cover each fregment on average |
|region_copy_numbers	| depth divided by the normal genome depth |
|junction_supports	| for each place where one fregment connects to the next | 
|path	| the order of the fregment going around the circle |
|A_breakpoint	| for each join, the exact position where the first fregment ends |
|B_breakpoint | the exact position where the next fregment begins |

If no exact position could be found for a join, that entry shows NA.
```
bigcycle_id	total_size	fragment_count	region_depths	region_copy_numbers	junction_supports	path	A_breakpoint	B_breakpoint
1	924140	2	chr2:15694017-16145119=220.45,chr2:16186561-16659597=221.72	chr2:15694017-16145119=223,chr2:16186561-16659597=224	chr2:16186561-16659597-->chr2:15694017-16145119-=80,chr2:15694017-16145119-->chr2:16186561-16659597-=187	chr2:16186561-16659597- -> chr2:15694017-16145119- -> chr2:16186561-16659597-	chr2:16224459,chr2:15694017	chr2:16035532,chr2:16659597

```

### summary_bigcycle.fa
 &gt; bigcycle_id_path <br> seq

```
>bigcycle_1_chr2:16186561-16659597-;chr2:15694017-16145119-;chr2:16186561-16659597-
TCATTTTCCCCTTCAATGATGGGGCCCCAAAATTTAAAGTACTTGTTCGGGCCCTTAAACCACAAAAAACAAAATTTTTTATT...
```


### login_not.txt
```
2026-10-19 09:49:09,972 [INFO] Logger initialized
2026-10-19 09:49:09,974 [INFO] Starting pipeline for sample SRR19746237
2026-10-19 09:49:09,974 [INFO] Estimating depth & minsupp from 100 random regions
...
2026-10-19 10:09:28,278 [INFO] Wrote 2 contigs -> SRR19746237_fec_v5_out14_mq30/summary_bigcycle.fa
2026-10-19 10:09:28,278 [INFO]  END: Building big-cycle FASTA (0.0s)
2026-10-19 10:09:28,278 [INFO] Pipeline completed successfully
```

## Additional information
1. [pysam](https://pysam.readthedocs.io/en/latest/api.html)
2. [bsalign](https://github.com/ruanjue/bsalign)


