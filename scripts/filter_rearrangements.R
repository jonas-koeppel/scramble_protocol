#!/usr/bin/env Rscript

# ==============================================================================
# SCRIPT: find_cre_rearrangements.R
#
# DESCRIPTION:
#   This script processes structural variant (SV) calls from NanoMonSV to
#   identify Cre-recombinase induced rearrangements. It works by finding SVs
#   whose breakpoints overlap with a given set of insertion sites (e.g., loxP
#   or LINE1 elements) and then uses a set of "liftover" sites to refine the
#   breakpoint coordinates.
#
# REQUIREMENTS:
#   R packages: argparse, tidyverse, GenomicRanges, plyranges.
#
#   To install packages, run this in an R session:
#   install.packages(c("argparse", "tidyverse"))
#   if (!require("BiocManager", quietly = TRUE))
#       install.packages("BiocManager")
#   BiocManager::install(c("GenomicRanges", "plyranges"))
# ==============================================================================

# --- 1. Load Libraries ---
suppressPackageStartupMessages({
  library(argparse)
  library(tidyverse)
  library(GenomicRanges)
  library(plyranges)
})

# --- 2. Define Command-Line Arguments ---
parser <- ArgumentParser(description="Find Cre-induced rearrangements from NanoMonSV output.")

parser$add_argument("-r", "--rearrangements", required=TRUE,
                    help="Path to the raw rearrangements TSV file from NanoMonSV.")
parser$add_argument("-i", "--insertion_sites", required=TRUE,
                    help="Path to the insertion sites TSV file (e.g., LINE1). Requires 'chr', 'start', 'end' columns.")
parser$add_argument("-l", "--liftover_sites", required=TRUE,
                    help="Path to the clonal insertion sites for coordinate liftover. Requires 'chr', 'start', 'end' columns.")

args <- parser$parse_args()


# --- 3. Define Global Variables and Functions ---

# List of standard chromosomes to filter against
chr_list <- sprintf("chr%s", c(seq(1, 22, 1), "X", "Y"))

# Annotations to classify NanoMonSV outputs into standard SV types
sv_annotations <- tibble(
  orientation = c("+-c", "++c", "--c", "-+c", "+-t", "++t", "--t", "-+t", "+-l", "++l", "--l"),
  SVTYPE = c("DEL", "INVh2h", "INVt2t", "CIRC", "TRAh2t", "TRAh2h", "TRAt2t", "TRAt2h", "INS", "FBh2h", "FBt2t"),
  category = c("Deletion", "Inversion", "Inversion", "Circle", "Translocation", "Translocation", "Translocation", "Translocation", "Insertion", "Fold back", "Fold back")
)

#' Read and process a NanoMonSV output file.
#' @param path File path to the NanoMonSV TSV.
#' @return A tibble of processed structural variants.
read_nanomonsv <- function(path) {
  read_tsv(path, col_types = cols()) %>%
    mutate(
      type = ifelse(Chr_1 != Chr_2, "t", ifelse(abs(Pos_2 - Pos_1) < 50, "l", "c")),
      orientation = paste0(Dir_1, Dir_2, type)
    ) %>%
    left_join(sv_annotations, by = "orientation") %>%
    mutate(
      SVLEN = ifelse(!str_detect(SVTYPE, "TRA"), Pos_2 - Pos_1, 0)
    ) %>%
    filter(
      Chr_1 %in% chr_list,
      Chr_2 %in% chr_list,
      type %in% c("l", "t") | SVLEN > 3000
    ) %>%
    dplyr::select(
      "chr" = "Chr_1", "start" = "Pos_1", "chr_2" = "Chr_2", "end" = "Pos_2",
      SVLEN, category, SVTYPE, orientation, "ins_seq" = "Inserted_Seq",
      "supp_reads" = "Supporting_Read_Num_Tumor"
    )
}

#' Harmonize SV coordinates by "lifting over" to precise insertion sites.
#' This corrects for potential inaccuracies in breakpoint detection.
#' @param x A tibble of SVs to be corrected.
#' @param liftover_gr A GRanges object of high-confidence insertion sites.
#' @return A tibble with corrected start and end coordinates.
liftover_coordinates <- function(x, liftover_gr) {
  message("Harmonizing coordinates...")
  x_filtered <- filter(x, SVTYPE != "INS", chr %in% chr_list, chr_2 %in% chr_list) %>%
    mutate(id = 1:nrow(.))

  if (nrow(x_filtered) == 0) return(tibble())

  # Liftover start coordinates
  sv_starts_gr <- dplyr::select(x_filtered, chr, start) %>%
    mutate(end = start, sv_start = start) %>%
    GRanges()
  sv_liftover_start <- join_overlap_left(liftover_gr, sv_starts_gr) %>%
    as_tibble() %>%
    filter(!is.na(sv_start)) %>%
    dplyr::select(chr = seqnames, start, sv_start) %>% distinct() %>%
    left_join(dplyr::select(x_filtered, chr, "sv_start" = "start", chr_2, end, SVTYPE, category, id),
              by = c("chr", "sv_start"), relationship = "many-to-many")

  # Liftover end coordinates
  sv_ends_gr <- dplyr::select(x_filtered, "chr" = "chr_2", "start" = "end") %>%
    mutate(end = start, sv_end = start) %>%
    GRanges()
  sv_liftover_end <- join_overlap_left(liftover_gr, sv_ends_gr) %>%
    as_tibble() %>%
    filter(!is.na(sv_end)) %>%
    dplyr::select(chr_2 = seqnames, end = start, sv_end) %>% distinct()

  # Combine lifted start and end points
  sv_liftover <- left_join(sv_liftover_end,
                           dplyr::select(sv_liftover_start, chr, start, "sv_end" = "end", chr_2, end, SVTYPE, category, id),
                           by = c("chr_2", "sv_end"), relationship = "many-to-many") %>%
    filter(!is.na(id)) %>%
    distinct()

  return(sv_liftover)
}

#' Main function to find rearrangements associated with insertion sites.
#' @param sv_data A tibble of SVs from read_nanomonsv.
#' @param insertion_gr A GRanges object of sites to check for overlaps.
#' @param liftover_gr A GRanges object for coordinate correction.
#' @param output_path Path for the final output TSV file.
find_rearrangements_nano <- function(sv_data, insertion_gr, liftover_gr, output_path) {
  message("Detecting rearrangements linked to insertion sites...")

  # Find SVs with inserted sequences (evidence of mediation) at breakpoints
  breakpoint_sites <- sv_data %>%
    filter(str_detect(ins_seq, "TAACTTCGTAT|ATACGAAGTTA|AATGTACATTAT|GTATAATGTAC"))

  # Identify start/end coordinates that overlap with the provided insertion sites
  starts_gr <- breakpoint_sites %>% dplyr::select(chr, start) %>% mutate(end = start) %>% GRanges()
  ends_gr <- breakpoint_sites %>% dplyr::select("chr" = "chr_2", "start" = "end") %>% mutate(end = start) %>% GRanges()

  overlaps_start <- subsetByOverlaps(starts_gr, insertion_gr) %>% as.data.frame() %>% mutate(region = paste0(seqnames, ":", start)) %>% `$`(region)
  overlaps_end <- subsetByOverlaps(ends_gr, insertion_gr) %>% as.data.frame() %>% mutate(region = paste0(seqnames, ":", start)) %>% `$`(region)

  message("Filtering for variants with both breakpoints at insertion sites...")
  passing_rearrangements <- sv_data %>%
    mutate(region_start = paste0(chr, ":", start), region_end = paste0(chr_2, ":", end)) %>%
    filter(region_start %in% overlaps_start & region_end %in% overlaps_end, SVTYPE != "INS")

  # If no matching rearrangements, write an empty file and exit
  if (nrow(passing_rearrangements) == 0) {
    message("No Cre-induced rearrangements found.")
    write_tsv(passing_rearrangements, output_path)
    return()
  }

  # Use high-confidence sites to correct ("liftover") breakpoint coordinates
  passing_rearrangements_lifted <- liftover_coordinates(passing_rearrangements, liftover_gr = liftover_gr)

  # Update coordinates if a liftover was successful
  if (nrow(passing_rearrangements_lifted) > 0) {
      final_rearrangements <- passing_rearrangements %>%
        mutate(id = 1:n()) %>%
        left_join(dplyr::select(passing_rearrangements_lifted, id, "start_lift" = "start", "end_lift" = "end"), by = "id") %>%
        mutate(
          start = ifelse(!is.na(start_lift), start_lift, start),
          end = ifelse(!is.na(end_lift), end_lift, end)
        ) %>%
        dplyr::select(-id, -start_lift, -end_lift, -region_start, -region_end)
  } else {
      final_rearrangements <- passing_rearrangements %>% dplyr::select(-region_start, -region_end)
  }
  
  # Summarize reads for events that become identical after liftover
  final_rearrangements_sum <- final_rearrangements %>%
      group_by(chr, start, chr_2, end, SVTYPE, category, orientation) %>%
      summarise(supp_reads = sum(supp_reads), ins_seq = ins_seq[1], .groups = "drop")


  message(paste(nrow(final_rearrangements_sum), "Cre-induced rearrangements detected."))
  write_tsv(final_rearrangements_sum, output_path)
}


# --- 4. Main Execution Block ---

# Derive a sample name and define the output file path
samplename <- basename(args$rearrangements) %>% str_remove("\\.nanomonsv\\.result\\.txt$")
output_file <- paste0(samplename, "_cre_induced.tsv")

cat(paste("--- Starting Analysis for Sample:", samplename, "---\n"))
cat(paste("Input Rearrangements:", args$rearrangements, "\n"))
cat(paste("Insertion Sites:", args$insertion_sites, "\n"))
cat(paste("Liftover Sites:", args$liftover_sites, "\n"))
cat(paste("Output File:", output_file, "\n\n"))

# Step 1: Read all input files
rearrangements_data <- read_nanomonsv(args$rearrangements) %>% mutate(sample = samplename)
insertion_sites_gr <- read_tsv(args$insertion_sites, col_types = cols()) %>% GRanges()
liftover_sites_gr <- read_tsv(args$liftover_sites, col_types = cols()) %>% GRanges()

# Step 2: Run the main analysis function
find_rearrangements_nano(
  sv_data = rearrangements_data,
  insertion_gr = insertion_sites_gr,
  liftover_gr = liftover_sites_gr,
  output_path = output_file
)

cat("\n--- Script finished successfully. ---\n")