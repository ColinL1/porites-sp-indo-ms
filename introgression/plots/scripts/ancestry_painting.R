#!/usr/bin/env Rscript
# Ancestry painting along the genome: at the sites where the two parents differ, the share of sites
# where Porites sp. Indonesia is homozygous for the P. lobata allele, heterozygous with one allele of
# each parent (the F1 pattern) or homozygous for the P. harrisoni allele, in 1 Mb bins on the
# concatenated scaffold axis. Built from the 100 kb windows of 04_ancestry_painting (counts summed).
# Uses the site patterns, not the per-haplotype tracks: phase blocks (NG50 ~37 kb) are shorter
# than the bins, so "haplotype 1" is not the same haplotype from one block to the next.
# usage: ancestry_painting.R <results dir> <out prefix>
# writes <prefix>.pdf/.png, <prefix>.tsv (plotted values), <prefix>.rds (panel, for the composite)
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(script_dir, "setup_introgression.R"))

args <- commandArgs(trailingOnly = TRUE)
res <- args[1]; prefix <- args[2]
bin_size <- 1e6
min_sites <- 50   # bins with fewer patterned sites are left empty

ax <- genome_axis(file.path(res, "00_inputs", "scaffolds.bed"))
w <- read_tsv(file.path(res, "04_ancestry_painting", "window_ancestry.tsv"), show_col_types = FALSE)

d <- w |>
  mutate(bin = floor(win_start / bin_size) * bin_size) |>
  group_by(chrom, bin) |>
  summarise(hom_P1 = sum(hom_P1), het_opposite = sum(het_opposite), hom_P2 = sum(hom_P2), .groups = "drop") |>
  mutate(total = hom_P1 + het_opposite + hom_P2) |>
  filter(total >= min_sites) |>
  inner_join(ax, by = "chrom") |>
  mutate(bin_end = pmin(bin + bin_size, len), xmin = offset + bin, xmax = offset + bin_end)
write_tsv(select(d, chrom, bin_start = bin, bin_end, hom_P1, het_opposite, hom_P2, total),
          paste0(prefix, ".tsv"))

pattern_lab <- c(hom_P1 = "Homozygous P. lobata", het_opposite = "Heterozygous",
                 hom_P2 = "Homozygous P. harrisoni")
long <- d |>
  pivot_longer(c(hom_P1, het_opposite, hom_P2), names_to = "pattern", values_to = "n") |>
  mutate(pattern = factor(pattern, levels = names(pattern_lab))) |>
  arrange(chrom, bin, pattern) |>
  group_by(chrom, bin) |>
  mutate(ymax = cumsum(n) / total, ymin = lag(ymax, default = 0)) |> ungroup()

p <- ggplot(long) +
  # outline in the fill colour: closes the anti-aliasing seams between adjacent bins in PDF viewers
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = pattern, colour = pattern),
            linewidth = 0.1) +
  scale_fill_manual(values = c(hom_P1 = col_p1, het_opposite = col_both, hom_P2 = col_p2),
                    labels = label_expr(pattern_lab), name = NULL, aesthetics = c("fill", "colour")) +
  scale_x_genome_axis(ax) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1"), expand = c(0, 0)) +
  labs(x = "Scaffold", y = "Share of\ndiagnostic sites") +
  theme_fig() +
  theme(panel.grid = element_blank(), legend.position = "bottom", legend.justification = "left",
        legend.location = "panel", legend.key.size = unit(3, "mm"),
        legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(1, "mm"))

saveRDS(p, paste0(prefix, ".rds"))
save_fig(p, prefix, width = 7.5, height = 1.9)
