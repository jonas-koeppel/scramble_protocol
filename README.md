# scramble_protocol

Computational tutorial for the scrambling protocol (Nature Protocols Supplementary Information 4).

This repository provides:
- a **tutorial** (docs/tutorial.md) that reproduces the example analysis

The example dataset uses public ONT long-read sequencing from SRA:
- HAP1 parental: **SRR28738516**
- HAP1 scrambled clone: **SRR29347399**

## Quickstart (full tutorial)

1) Clone the repo
```bash
git clone git@github.com:jonas-koeppel/scramble_protocol.git
cd scramble_protocol
```

2) Create and activate the conda/mamba environment
```bash
mamba create -n scrambling_analysis -c bioconda -c conda-forge   python=3.10 sra-tools r-base=4.3 r-argparse r-tidyverse r-igraph   bioconductor-genomicranges bioconductor-plyranges   bioconductor-variantannotation bioconductor-structuralvariantannotation   samtools minimap2 sniffles mosdepth nanomonsv mafft
conda activate scrambling_analysis
```

3) Install `gintools` in R
```r
install.packages("remotes")
remotes::install_github("cnobles/gintools")
```

4) Follow the tutorial
- **Main tutorial:** `docs/tutorial.md`
- **Reference setup details:** `resources/README.md`

## Repository layout

- `docs/tutorial.md` — step-by-step commands corresponding to Supplementary Information 4
- `scripts/` — R scripts used by the tutorial (e.g., `find_insertions.R`, `filter_rearrangements.R`)
- `resources/` — reference download + checksum verification instructions
- `expected/` — expected outputs

## Tested compute environments

This tutorial has been validated on:
- Intel Xeon Gold 6312U @ 2.40GHz, Ubuntu 22.04.3 LTS, 20 GB RAM
- 30-core system with T4 GPU, Ubuntu 20.04, 120 GB RAM
