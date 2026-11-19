#!/usr/bin/env Rscript
# Protein domains whose count in the new assembly differs most from the other Porites
# (funannotate compare pfam/pfam.results.csv or interpro/interproscan.results.csv: domain x genome counts).
# Score = log2((count in Porites sp. + 1) / (median count in the other Porites + 1)); domains with a max count
# >= --min-count across Porites are ranked and the top N in each direction shown as a count heatmap.
# Domains whose description matches transposable-element keywords are marked with a dagger: their counts
# depend strongly on how each annotation pipeline treated repeats.
# usage: domain_counts.R <results csv> <out prefix> [top N per direction, default 15] [min count, default 10]
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(script_dir, "..", "..", "..", "shared", "plots", "setup.R"))   # names, colours, genome order, shared functions

args <- commandArgs(trailingOnly = TRUE)
in_csv <- args[1]; prefix <- args[2]
n_top <- if (length(args) >= 3) as.integer(args[3]) else 15L
min_count <- if (length(args) >= 4) as.numeric(args[4]) else 10

te_regex <- regex(paste("transpos", "reverse transcriptase", "integrase", "retrotrans", "RNase H",
                        "H2C2 zinc finger", "Tc5", "Tc1", "piggyBac", "\\bDDE\\b", "\\bhAT\\b", "gag",
                        "CENP-?B", "Jerky", "Pol polyprotein", "Psq", "integrase-like", "Pao", sep = "|"),
                  ignore_case = TRUE)

d <- read_fc_table(in_csv)
others <- setdiff(porites, new_genome)
scored <- d |> select(id, descriptions, all_of(porites)) |>
  mutate(across(all_of(porites), as.numeric)) |>
  rowwise() |> mutate(median_other = median(c_across(all_of(others))),
                      max_count = max(c_across(all_of(porites)))) |> ungroup() |>
  mutate(log2r = log2((.data[[new_genome]] + 1) / (median_other + 1)),
         te = str_detect(descriptions, te_regex) &
           !str_detect(descriptions, regex("Half-A-TPR|Exuperantia", ignore_case = TRUE))) |>
  filter(max_count >= min_count)
write_tsv(scored |> arrange(log2r), paste0(prefix, ".scores.tsv"))

sel <- bind_rows(scored |> filter(log2r > 0) |> slice_max(log2r, n = n_top, with_ties = FALSE),
                 scored |> filter(log2r < 0) |> slice_min(log2r, n = n_top, with_ties = FALSE)) |>
  mutate(direction = ifelse(log2r > 0, paste("More in", new_genome_label), paste("Fewer in", new_genome_label)),
         label = paste0(str_trunc(str_squish(descriptions), 55), ifelse(te, " †", ""), "  ", id),
         label = fct_reorder(label, log2r))
message(sprintf("%d domains pass min count %g; showing %d (%d TE-related)", nrow(scored), min_count,
                nrow(sel), sum(sel$te)))

long <- sel |> select(label, direction, all_of(porites)) |>
  pivot_longer(all_of(porites), names_to = "genome", values_to = "count") |>
  mutate(genome = factor(genome, levels = porites))

fig <- ggplot(long, aes(x = genome, y = label, fill = log10(count + 1))) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = count, colour = log10(count + 1) > 1.6), size = 2.3) +
  facet_grid(direction ~ ., scales = "free_y", space = "free_y", labeller = labeller_italic) +
  scale_fill_gradientn(colours = seq_ramp, name = "domain count", breaks = log10(c(1, 11, 101, 1001)),
                       labels = c(0, 10, 100, 1000)) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = ink), guide = "none") +
  scale_x_genome() +
  labs(x = NULL, y = NULL) +
  theme_fig() + theme(axis.text.x = element_text(angle = 30, hjust = 1),
                      axis.text.y = element_text(size = 7), panel.grid = element_blank(),
                      strip.text.y = element_text(angle = -90))

save_fig(fig, prefix, width = 7.5, height = 2 + nrow(sel) * 0.17)
