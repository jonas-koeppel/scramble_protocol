# Supplementary Information 4: Example dataset and reproducibility run

This document reproduces the computational tutorial (Supplementary Information 4) end-to-end.

## Data inputs

Two SRA runs are used:

- **HAP1 parental (SRR28738516)**: `.sra` download is ~55.6 GB; extracted FASTQ is larger. Plan for **>3×** the `.sra` size for the FASTQ plus additional temporary space during extraction (see `-t` below).
- **HAP1 scrambled clone (SRR29347399)**: `.sra` download is ~11.9 GB; FASTQ is larger for the same reason.

### Download SRA runs

Repeat `prefetch` if it reports incomplete files.

```bash
prefetch SRR28738516
prefetch SRR29347399
```

### Convert to FASTQ

Set `-t` to a large local scratch temp folder.

```bash
fasterq-dump SRR28738516 -e 8 -t <tmp_folder> -p
fasterq-dump SRR29347399 -e 8 -t <tmp_folder> -p
```

Expected outputs:
- `SRR28738516.fastq`
- `SRR29347399.fastq`

## Reference files

See `resources/reference/README.md` for checksum verification and indexing details.

### Download hg38 + tandem repeats (TRF)

```bash
# Human reference genome (hg38)
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz

# Tandem repeats
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.trf.bed.gz
gunzip hg38.trf.bed.gz

# Organize
mkdir -p genome
mv -f hg38.fa genome/hg38.fa
mv -f hg38.trf.bed genome/hg38.trf.bed
```

## Software requirements

A dedicated conda/mamba environment is recommended to avoid R/Bioconductor conflicts.

```bash
mamba create -n scrambling_analysis -c bioconda -c conda-forge   python=3.10 sra-tools r-base=4.3 r-argparse r-tidyverse r-igraph   bioconductor-genomicranges bioconductor-plyranges   bioconductor-variantannotation bioconductor-structuralvariantannotation   samtools minimap2 sniffles mosdepth nanomonsv mafft
conda activate scrambling_analysis
```

Install `gintools` in R:
```r
install.packages("remotes")
remotes::install_github("cnobles/gintools")
```

## Tested compute environment

Validated on two independent systems:
- Intel Xeon Gold 6312U @ 2.40GHz, Ubuntu 22.04.3 LTS, 20 GB RAM
- 30-core system with T4 GPU, Ubuntu 20.04, 120 GB RAM

Timings below reflect the Intel Xeon setup and are for guidance only.

## Procedure

### 1) Align ONT reads to hg38 (minimap2)

Scrambled clone:
```bash
minimap2 -ax map-ont -t 8 genome/hg38.fa SRR29347399.fastq > clone17.sam
# (8846 seconds)
```

Parental:
```bash
minimap2 -ax map-ont -t 8 genome/hg38.fa SRR28738516.fastq > parental.sam
# (73162 seconds)
```

### 2) Convert, sort and index (samtools)

```bash
samtools sort -@ 10 -o clone17.bam clone17.sam
# (7m35.642s)

samtools sort -@ 10 -o parental.bam parental.sam
# (72m10.083s)

samtools index -@ 10 clone17.bam
samtools index -@ 10 parental.bam
```

### 3) Call SVs in the parental for loxP mapping (Sniffles2)

```bash
mkdir -p vcf

sniffles --input parental.bam   --vcf vcf/parental.vcf.gz   --reference genome/hg38.fa   --tandem-repeats genome/hg38.trf.bed   --threads 8   --mosaic   --minsupport 3   --minsvlen 25   --qc-output-all   --mapq 5

# Output: vcf/parental.vcf.gz
# (8m11.405s)
```

### 4) Map loxP insertions in the parental (R script)

> This step is only necessary for the parental.

```bash
find_insertions.R --path vcf/parental.vcf.gz --site loxPsym --ploidy 1
```

Outputs:
- `parental_all.tsv`
- `parental_clonal.tsv`

Reference run result:
- 1012 loxPsym insertions detected in the parental clone, of which 294 are clonal.
- (37 seconds on a Macbook Pro M3 Pro, 18 GB RAM)

### 5) Variant calling with nanomonsv (using parental as control)

```bash
nanomonsv parse clone17.bam ./nanomonsv/clone17
# (6m37.292s)

nanomonsv parse parental.bam ./nanomonsv/parental
# (30m2s)

nanomonsv get nanomonsv/clone17 clone17.bam genome/hg38.fa   --control_prefix nanomonsv/parental --control_bam parental.bam   --min_tumor_variant_read_num 2 --min_indel_size 1000   --processes 16 --min_tumor_VAF 0.02   --max_control_variant_read_num 0 --check_read_max_num 200

# Output: clone17.nanomonsv.result.txt
# (47m38s)
```

### 6) Filter to Cre-induced variants (R script)

```bash
filter_rearrangements.R   --rearrangements clone17.nanomonsv.result.txt   --insertion_sites LINE1_nick_3mm.tsv   --liftover_sites parental_all.tsv
```

Output:
- `clone17_cre_induced.tsv`

Reference run result:
- two Cre-induced inversions and one deletion
- (0.85 seconds on a Macbook Pro M3 Pro, 18 GB RAM)

### 7) Compute coverages (mosdepth)

```bash
mosdepth -n --by 50000 -t 8 clone17 clone17.bam
# (0m58.132s)

mosdepth -n --by 50000 -t 8 parental parental.bam
# (2m52.737s)
```