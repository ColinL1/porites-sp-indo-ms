#!/usr/bin/env Rscript
# Orthogroup sharing across genomes (funannotate compare orthology/orthology_groups.txt, ProteinOrtho),
# drawn as an UpSet plot: bars = orthogroups with exactly that set of genomes, matrix = which genomes,
# side bars = orthogroups per genome. The Acropora outgroup is kept here: it separates families shared
# with Acropora from Porites-specific ones. Proteins in no orthogroup are not in this file (see
# genome_stats, "Proteins without orthologs").
# usage: orthology_upset.R <funannotate_compare dir> <out prefix> [number of intersections, default 20]
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(script_dir, "..", "..", "..", "shared", "plots", "setup.R"))   # names, colours, genome order, shared functions

args <- commandArgs(trailingOnly = TRUE)
fc_dir <- args[1]; prefix <- args[2]
n_top <- if (length(args) >= 3) as.integer(args[3]) else 20L

og <- read_tsv(file.path(fc_dir, "orthology", "orthology_groups.txt"), col_names = FALSE,
               col_types = cols(.default = "c"))
members <- og |> select(og = X1, member = X5) |>
  separate_rows(member, sep = ",\\s*") |>
  # "Porites.lobata-gene-PLOB.00028857-T1" -> Porites_lobata; the new assembly's IDs carry no prefix
  mutate(genome = ifelse(str_detect(member, "^[A-Z][a-z]+\\.[a-z]+-"),
                         str_replace(str_extract(member, "^[A-Z][a-z]+\\.[a-z]+"), "\\.", "_"), new_genome))
stopifnot(all(members$genome %in% all_genomes))

pattern <- members |> distinct(og, genome) |>
  mutate(genome = factor(genome, levels = all_genomes)) |>
  arrange(og, genome) |>
  group_by(og) |> summarise(set = paste(genome, collapse = "&"), n_genomes = n(),
                            has_new = new_genome %in% genome, .groups = "drop")

inter <- pattern |> count(set, n_genomes, has_new, name = "n") |> arrange(desc(n)) |>
  mutate(rank = row_number())
write_tsv(inter |> select(-rank), paste0(prefix, ".intersections.tsv"))
# bar colour classes (legend labels)
k_all <- "Shared by all genomes"; k_porites <- "Shared by all Porites"
k_with <- paste("Including", new_genome_label); k_without <- paste("Without", new_genome_label)
top <- inter |> slice_head(n = n_top) |>
  mutate(x = factor(rank),
         kind = factor(case_when(n_genomes == length(all_genomes) ~ k_all,
                                 set == paste(porites, collapse = "&") ~ k_porites,
                                 has_new ~ k_with, TRUE ~ k_without),
                       levels = c(k_all, k_porites, k_with, k_without)))
sizes <- pattern |> separate_rows(set, sep = "&") |> count(genome = set, name = "n")
message(sprintf("%d orthogroups, %d intersections; top %d cover %.1f%%", nrow(pattern), nrow(inter),
                n_top, 100 * sum(top$n) / nrow(pattern)))

kind_col <- setNames(c(ink, "#5a5a55", col_new, col_other), c(k_all, k_porites, k_with, k_without))
x_scale <- scale_x_discrete(limits = levels(top$x))

p_bars <- ggplot(top, aes(x = x, y = n, fill = kind)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = format(n, big.mark = ",")), vjust = -0.4, size = 2.2, colour = ink_muted) +
  # legend: only the classes the matrix does not already make obvious
  scale_fill_manual(values = kind_col, breaks = c(k_porites, k_with), labels = function(x) label_expr(x),
                    name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12)), labels = scales::label_comma()) +
  x_scale + labs(x = NULL, y = "Orthogroups") +
  theme_fig() + theme(axis.text.x = element_blank(), panel.grid.major.x = element_blank(),
                      # legend inside the bar panel, top right, above the small bars
                      legend.position = "inside", legend.position.inside = c(1, 1),
                      legend.justification = c(1, 1), legend.key.size = unit(3.5, "mm"),
                      legend.text = element_text(size = 8))

grid_df <- top |> select(x, set) |> crossing(genome = all_genomes) |>
  mutate(present = str_detect(paste0("&", set, "&"), paste0("&", genome, "&")))
lines_df <- grid_df |> filter(present) |> group_by(x) |>
  summarise(lo = min(match(genome, rev(all_genomes))), hi = max(match(genome, rev(all_genomes))))
p_matrix <- ggplot(grid_df, aes(x = x, y = factor(genome, levels = rev(all_genomes)))) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = length(all_genomes) - 0.5,
           ymax = length(all_genomes) + 0.5, fill = band) +
  geom_point(data = filter(grid_df, !present), colour = grid_col, size = 2.2) +
  geom_segment(data = lines_df, aes(x = x, xend = x, y = lo, yend = hi), inherit.aes = FALSE,
               colour = ink, linewidth = 0.6) +
  geom_point(data = filter(grid_df, present), colour = ink, size = 2.2) +
  x_scale + scale_y_genome() + labs(x = NULL, y = NULL) +
  theme_fig() + theme(axis.text.x = element_blank(), panel.grid = element_blank())

p_sizes <- ggplot(sizes, aes(y = factor(genome, levels = rev(all_genomes)), x = n,
                             fill = genome == new_genome)) +
  geom_col(width = 0.6) +
  scale_fill_manual(values = c(`TRUE` = col_new, `FALSE` = col_other), guide = "none") +
  scale_x_reverse(labels = scales::label_number(scale_cut = scales::cut_short_scale()),
                  breaks = c(0, 20000)) +
  labs(x = "Total orthogroups", y = NULL) +
  theme_fig() + theme(axis.text.y = element_blank(), panel.grid.major.y = element_blank())

fig <- (plot_spacer() + p_bars + p_sizes + p_matrix) +
  plot_layout(ncol = 2, widths = c(1, 4.5), heights = c(3, 1.6))

save_fig(fig, prefix, width = 7.5, height = 4.4)
saveRDS(fig, paste0(prefix, ".rds"))   # for composite.R
