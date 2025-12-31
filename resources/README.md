# Reference resources (hg38 + tandem repeats)

This folder contains input files used by the tutorial. The hg38.fa reference genome is too large to include here but the correctness of the downloaded file can be verified by comparing the checksum to hg38.fa.sha256 within this folder

The tutorial uses:
- UCSC hg38 FASTA: `hg38.fa`

## Download

From the repository root:

```bash
mkdir -p genome
cd genome

wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz

wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.trf.bed.gz
gunzip hg38.trf.bed.gz
```

This should produce:
- `genome/hg38.fa`
- `genome/hg38.trf.bed`

## Verify with checksums

Generate checksums after download (Linux):
```bash
sha256sum genome/hg38.fa genome/hg38.trf.bed > genome.reference.sha256
```

macOS:
```bash
shasum -a 256 genome/hg38.fa genome/hg38.trf.bed > genome.reference.sha256
```
