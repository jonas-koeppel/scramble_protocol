#!/usr/bin/env Rscript

# ==============================================================================
# SCRIPT: filter_vcf_for_loxP.R
#
# DESCRIPTION:
#   This script processes a VCF file (e.g., from Sniffles) to identify loxP or
#   loxPsym sites. It takes the VCF path, site type, and sample ploidy as
#   input. It automatically generates two output files in the current directory:
#   1. {samplename}_all.tsv: All insertions matching the site sequence.
#   2. {samplename}_clonal.tsv: A subset of insertions that pass clonal filters.
#
# REQUIREMENTS:
#   R packages: argparse, VariantAnnotation, StructuralVariantAnnotation,
#               tidyverse, gintools.
# ==============================================================================

# --- 1. Load Libraries ---
suppressPackageStartupMessages({
  library(argparse)
  library(VariantAnnotation)
  library(StructuralVariantAnnotation)
  library(tidyverse)
  library(gintools)
})

# --- 2. Define Command-Line Arguments ---
parser <- ArgumentParser(description="Filter a Sniffles VCF for loxP/loxPsym insertions.")

parser$add_argument("-p", "--path", required=TRUE,
                    help="Path to the input VCF file (.vcf or .vcf.gz).")
parser$add_argument("-s", "--site", required=TRUE, choices=c("loxP", "loxPsym"),
                    help="The specific loxP site to filter for: 'loxP' or 'loxPsym'.")
parser$add_argument("-d", "--ploidy", required=TRUE, type="integer", choices=c(1, 2, 3),
                    help="Ploidy of the sample for setting the allele frequency threshold (1, 2, or 3).")

args <- parser$parse_args()


# --- 3. Define Global Variables and Functions ---

chr_list <- sprintf("chr%s", c(seq(1, 22, 1), "X", "Y"))
annotations_sniffles <- tibble(
  SVTYPE = c("DEL", "INV", "DUP", "BND", "INS"),
  category = c("Deletion", "Inversion", "Circle", "Translocation", "Insertion")
)

#' Load and process a Sniffles VCF file into a tidy tibble.
load_sniffles <- function(path) {
  vcf <- VariantAnnotation::readVcf(path, "GRCh38")
  gr <- GRanges(
    vcf@rowRanges, chr_2 = vcf@info$CHR2, END = vcf@info$END, SVLEN = vcf@info$SVLEN,
    SVTYPE = vcf@info$SVTYPE, STRANDS = vcf@info$STRAND, REF = vcf@fixed$REF,
    ALT = vcf@fixed$ALT, FILTER = vcf@fixed$FILTER, AF = vcf@info$VAF,
    supp_reads = vcf@info$SUPPORT, read_names = vcf@info$RNAMES,
    coverage = vcf@info$COVERAGE
  ) %>% gintools::unique_granges()
  
  df <- BiocGenerics::as.data.frame(gr) %>%
    dplyr::rename("chr" = "seqnames") %>% dplyr::select(-strand)
  
  read_names_flat <- sapply(df$read_names, paste, collapse = "_")
  
  tibble(df) %>%
    mutate(
      chr_2 = ifelse(is.na(chr_2), as.character(chr), chr_2),
      read_names = read_names_flat, ALT = as.character(unlist(ALT))
    ) %>%
    mutate(ALT = str_remove_all(ALT, "\\[|\\]|N")) %>%
    separate(ALT, into = c("ALT_chr", "ALT_start"), sep = ":", remove = FALSE, fill = "right") %>%
    mutate(
      ALT_start = as.numeric(ALT_start),
      END = ifelse(is.na(END), ALT_start, END),
      chr_2 = ifelse(is.na(chr_2), ALT_chr, chr_2)
    ) %>%
    dplyr::select(
      chr, start, chr_2, "end" = END, SVLEN, SVTYPE, STRANDS,
      ALT, AF, supp_reads, read_names, coverage
    ) %>%
    left_join(annotations_sniffles, by = "SVTYPE") %>%
    distinct() %>%
    filter(chr %in% chr_list, chr_2 %in% chr_list)
}

#' Filter for insertions containing the loxPsym sequence.
filter_vcf_loxPsym <- function(x) {
  x[agrep("TAACTTCGTATAATGTACATTATACGAAGTTA", x$ALT, max.distance = 2), ] %>%
    filter(SVTYPE == "INS")
}

#' Filter for insertions containing the loxP sequence.
filter_vcf_loxP <- function(x) {
  bind_rows(
    x[agrep("TAACTTCGTATAGCATACATTATACGAAGTTA", x$ALT, max.distance = 2), ],
    x[agrep("TAACTTCGTATAATGTATGCTATACGAAGTTA", x$ALT, max.distance = 2), ]
  ) %>% filter(SVTYPE == "INS")
}


# --- 4. Main Execution Block ---

# Derive sample name and define output file paths
samplename <- basename(args$path) %>% str_remove("\\.vcf(\\.gz)?$")
output_all_path <- paste0(samplename, "_all.tsv")
output_clonal_path <- paste0(samplename, "_clonal.tsv")

cat(paste("--- Starting Analysis for Sample:", samplename, "---\n"))
cat(paste("Input VCF:", args$path, "\n"))
cat(paste("Site Type:", args$site, "\n"))
cat(paste("Ploidy:", args$ploidy, "\n\n"))

# Step 1: Load and process the VCF file
cat("Step 1: Loading and parsing VCF file...\n")
variants <- load_sniffles(args$path) %>%
  mutate(sample = samplename)

# Step 2: Filter for the specified loxP site type
cat(paste("Step 2: Filtering for", args$site, "insertions...\n"))
loxP_insertions <- if (args$site == "loxP") {
  filter_vcf_loxP(variants)
} else {
  filter_vcf_loxPsym(variants)
}

# Step 3: Write ALL identified insertions to a file
cat(paste("   -> Found", nrow(loxP_insertions), "total insertions. Writing to", output_all_path, "\n"))
write_tsv(loxP_insertions, output_all_path)

# Step 4: Determine allele frequency (AF) threshold based on ploidy
thr <- switch(as.character(args$ploidy), "1" = 0.8, "2" = 0.3, "3" = 0.1)

# Step 5: Apply filters for clonal insertions
cat(paste("Step 3: Filtering for clonal events (AF >", thr, "and reads > 3)...\n"))
clonal_insertions <- dplyr::filter(
  loxP_insertions,
  AF > thr,
  supp_reads > 3
)

# Step 6: Write CLONAL insertions to a separate file
cat(paste("   -> Found", nrow(clonal_insertions), "clonal insertions. Writing to", output_clonal_path, "\n"))
write_tsv(clonal_insertions, output_clonal_path)

cat("\n--- Script finished successfully. ---\n")