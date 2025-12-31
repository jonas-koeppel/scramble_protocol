library(tidyverse)
library(circlize)

setwd("/path/to/scramble_protocol")

# Insertion sites
loxPsym_insertions <- read_tsv("./expected/parental_all.tsv")

chr_sizes <- read_tsv("./resources/GRCh38.chrom_sizes.txt", col_names = c("tmp", "chr", "end")) %>% mutate(start = 0, supp_reads = 0) %>% filter(chr != "chrY") %>%
  pivot_longer(c(start, end), names_to = "label", values_to = "start") %>% dplyr::select(chr, start, supp_reads) 

p <- loxPsym_insertions %>% bind_rows(mutate(chr_sizes, clone = "HAP1 loxPsym")) %>% mutate(chr = str_remove(chr, "chr")) %>%
  ggplot(aes(x = start/1000000, y = supp_reads)) +
  geom_col(size = 0.15, position = "stack", show.legend = F, color = "steelblue", fill = "steelblue") +
  facet_grid(cols = vars(factor(chr, levels = c(1:22, "X"))), scales = ("free"), space = ("free_x")) +
  #coord_cartesian(ylim = c(0, 60)) +
  scale_y_continuous("Number of supporting reads", breaks = c(0, 20, 40, 60), expand = c(0,0)) +
  scale_x_continuous("Chromosome position [Mbp]") +
  theme_sv +
  theme(axis.text = element_text(color = "black"),
        panel.background = element_rect(fill = "white", color = "black", size = 0.2),
        strip.background = element_rect(color = "white", fill = "white"),
        legend.position = "bottom",
        text = element_text(color = "black"),
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        panel.spacing = unit(0.05, "lines"))
ggsave("./expected/Insertions.pdf", p, width = 20, height = 5, units = "cm", dpi = 300)

# Circos track
clone17_deletions = read_tsv("./expected/clone17_cre_induced.tsv")

deletions_1 <- clone17_deletions %>% filter(category == "Deletion") %>% dplyr::select(chr, start) %>% mutate(end = start+1)
deletions_2 <- clone17_deletions %>% filter(category == "Deletion") %>% dplyr::select("chr" = "chr_2", "start" = "end") %>% mutate(end = start+1)
inversions_1 <- clone17_deletions %>% filter(category == "Inversion") %>% dplyr::select(chr, start) %>% mutate(end = start+1)
inversions_2 <- clone17_deletions %>% filter(category == "Inversion") %>% dplyr::select("chr" = "chr_2", "start" = "end") %>% mutate(end = start+1)

insertion_reads <- loxPsym_insertions %>% dplyr::select(chr, start, end, supp_reads)

pdf("./expected/Circos_clone17.pdf", width = unit(3, "cm"), height = unit(3, "cm"))
circos.clear()
circos.par(track.height = 0.2)
circos.initializeWithIdeogram(species = "hg38", plotType = c("ideogram", "labels"))
circos.genomicTrack(insertion_reads, panel.fun = function(region, value, ...) {
  circos.genomicLines(region, value, type = "h", lwd = 0.5, col = "#555555")
})
circos.genomicLink(deletions_1, deletions_2, col = "#0074A8", lwd = 0.5)
circos.genomicLink(inversions_1, inversions_2, col = "#E8753D", lwd = 0.5)
dev.off()


# Deletion coverage track
parental_coverage = read_tsv("./expected/parental.regions.bed", col_names = c("chr", "start", "end", "coverage")) %>% mutate(sample = "parental")
clone17_coverage = read_tsv("./expected/clone17.regions.bed", col_names = c("chr", "start", "end", "coverage")) %>% mutate(sample = "clone17")

p <- bind_rows(parental_coverage, clone17_coverage) %>% filter(chr == "chr3", start > 107510791, start < 111510791) %>%
  group_by(sample) %>%
  mutate(relative_coverage = coverage/mean(coverage)) %>%
  ggplot(aes(x = start/1000000, y = relative_coverage, col = sample)) +
  geom_point(size = 0.1) +
  annotate(geom = "rect", xmin = 109.510791, xmax = 109.762718, ymin = 0, ymax = 2.5, alpha = 0.1) +
  scale_color_manual(values = c("steelblue", "darkgray")) +
  labs(x = "Position on chr3 [Mb]", y = "Relative read coverage") +
  theme_sv
ggsave("./expected/Deletion_coverage.pdf", p, width = 10, height = 4, units = "cm", dpi = 300)

