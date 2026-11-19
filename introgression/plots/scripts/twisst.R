#!/usr/bin/env Rscript
# Twisst topology weights of the quartet (Porites sp. Indonesia, P. lobata, P. harrisoni, outgroup):
#   along: T1 - T2 weight difference along the genome for the along_outgroup (P. cylindrica: most windows, and the
#          outgroup of the robustness grid), 1 Mb bins on the concatenated scaffold axis
#   mean weights (one panel, stacked bars): per outgroup at default settings (06_twisst/twisst_summary.tsv)
#   above the robustness grid (07_twisst_robustness, P. cylindrica outgroup)
# Topologies: T1 = target with P. lobata, T2 = target with P. harrisoni, T3 = the species tree
# (P. lobata with P. harrisoni). Colours as in the ancestry painting (P. lobata blue, P. harrisoni red).
# "along" shows the directional signal of introgression, T1 - T2 per bin, rather than the three weights,
# so it does not repeat the form of the ancestry painting.
# usage: twisst.R <results dir> <out prefix>
# writes <prefix>.pdf/.png (the two panels), <prefix>.along.tsv / .mean_weights.tsv (plotted values),
# <prefix>.panels.rds (list of panels, for the composite)
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(script_dir, "setup_introgression.R"))

args <- commandArgs(trailingOnly = TRUE)
res <- args[1]; prefix <- args[2]
bin_size <- 1e6
min_windows <- 3   # bins with fewer Twisst windows are left empty
along_outgroup <- "Pcylindrica1"

topo_lab <- c(T1 = "P. sp. Indonesia + P. lobata", T3 = "P. lobata + P. harrisoni",
              T2 = "P. sp. Indonesia + P. harrisoni")
topo_col <- c(T1 = col_p1, T3 = col_both, T2 = col_p2)
og_lab <- sample_lab[outgroups]
scale_fill_topo <- function(...) scale_fill_manual(values = topo_col, labels = label_expr(topo_lab),
                                                   breaks = names(topo_lab), name = "Topology", ...)

# a weights header topology -> T1/T2/T3, read from the pair joined first in the unrooted quartet:
# (P1,P2) = T3; (P1,target) or (P2,outgroup) = T1; (P1,outgroup) or (P2,target) = T2
classify_topo <- function(nwk, og) {
  pair <- sort(str_split(str_match(nwk, "^\\(\\(([^)]+)\\)")[, 2], ",")[[1]])
  is <- function(a, b) identical(pair, sort(c(a, b)))
  if (is(p1, p2)) "T3" else if (is(p1, target) || is(p2, og)) "T1"
  else if (is(p1, og) || is(p2, target)) "T2" else stop("unknown topology: ", nwk)
}

ax <- genome_axis(file.path(res, "00_inputs", "scaffolds.bed"))

## ---- along the genome
read_og <- function(og) {
  dir <- file.path(res, "06_twisst", og)
  f <- file.path(dir, paste0(og, "_weights.csv.gz"))
  hdr <- grep("^#topo", readLines(f, n = 10), value = TRUE)
  map <- setNames(vapply(str_remove(hdr, "^#\\S+ "), classify_topo, "", og = og), str_match(hdr, "^#(\\S+)")[, 2])
  wt <- read_tsv(f, comment = "#", show_col_types = FALSE)
  names(wt) <- map[names(wt)]
  win <- read_tsv(file.path(dir, paste0(og, "_windows.tsv")), show_col_types = FALSE)
  stopifnot(nrow(win) == nrow(wt))
  bind_cols(win, wt) |>
    mutate(tot = T1 + T2 + T3) |> filter(tot > 0) |>
    mutate(across(c(T1, T2, T3), ~ .x / tot), outgroup = og,
           bin = floor((start + end) / 2 / bin_size) * bin_size)
}
along <- read_og(along_outgroup) |>
  group_by(outgroup, chrom, bin) |>
  summarise(n_windows = n(), across(c(T1, T3, T2), mean), .groups = "drop") |>
  filter(n_windows >= min_windows) |>
  inner_join(ax, by = "chrom") |>
  mutate(bin_end = pmin(bin + bin_size, len), xmin = offset + bin, xmax = offset + bin_end)
write_tsv(select(along, outgroup, chrom, bin_start = bin, bin_end, n_windows, T1, T3, T2), paste0(prefix, ".along.tsv"))

# T1 - T2 per bin: above 0 the new genome groups more often with P. lobata, below 0 with P. harrisoni
diff <- along |> mutate(d = T1 - T2, side = factor(ifelse(d >= 0, "T1", "T2"), levels = c("T1", "T2")))
lim <- max(abs(diff$d))
# legend: which parent the new genome groups with more often in the bin (T1 > T2: P. lobata); no title
scale_fill_closer <- scale_fill_manual(values = topo_col[c("T1", "T2")], labels = label_expr(c("P. lobata", "P. harrisoni")),
                                       name = NULL, aesthetics = c("fill", "colour"))
pAlong <- ggplot(diff) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = pmin(d, 0), ymax = pmax(d, 0), fill = side, colour = side),
            linewidth = 0.1) +   # outline in the fill colour: no seams between adjacent bins in PDF viewers
  geom_hline(yintercept = 0, colour = ink_muted, linewidth = 0.25) +
  scale_fill_closer +
  scale_x_genome_axis(ax) +
  scale_y_continuous(limits = c(-lim, lim), expand = expansion(mult = 0.05)) +
  labs(x = "Scaffold", y = "Topology weight\ndifference") +
  theme_fig() +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        legend.position = "bottom", legend.justification = "left", legend.location = "panel",
        legend.key.size = unit(3, "mm"), legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(1, "mm"))
# the same, turned 90 degrees (genome down the y axis, scaffold 1 at the top), for a narrow column;
# legend below, one entry per line
pAlongV <- ggplot(diff) +
  geom_rect(aes(ymin = xmin, ymax = xmax, xmin = pmin(d, 0), xmax = pmax(d, 0), fill = side, colour = side),
            linewidth = 0.1) +
  geom_vline(xintercept = 0, colour = ink_muted, linewidth = 0.25) +
  scale_fill_closer +
  scale_y_reverse(breaks = ax$mid, labels = ax$label, expand = expansion(mult = 0.005)) +
  scale_x_continuous(limits = c(-lim, lim), breaks = c(-1, -0.5, 0, 0.5, 1), expand = expansion(mult = 0.05)) +
  labs(y = "Scaffold", x = "Topology weight difference") +
  theme_fig() +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        legend.position = "bottom", legend.justification = "left", legend.location = "panel",
        legend.direction = "horizontal", legend.title.position = "top",
        legend.key.size = unit(3, "mm"), legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(1.5, "mm"), legend.spacing.y = unit(0.5, "mm"))

## ---- mean topology weights as stacked bars: the three weights of each run as shares of 100 %
# (the weights sum to ~1; dividing by the sum only absorbs rounding). Top: one bar per outgroup at the
# default settings (06_twisst/twisst_summary.tsv). Bottom: the robustness grid (07_twisst_robustness,
# P. cylindrica outgroup), one column per representation of the new genome, rows by window size with NJ
# above ML. Segments blue | grey | red, labelled with their rounded percentage; dashed lines at 1/3 and
# 2/3 = where the boundaries would be if the three topologies were equally frequent.
default_lab <- "Outgroups (default settings)"
rep_lab <- c(rand2 = "Pseudo-haplotypes", phased2 = "Phased haplotypes", cons1 = "IUPAC consensus")
to_shares <- function(d, ...) d |>
  pivot_longer(c(T1, T3, T2), names_to = "topo", values_to = "w") |>
  group_by(...) |> mutate(share = w / sum(w)) |> ungroup() |>
  mutate(topo = factor(topo, levels = names(topo_lab)))

gw <- read_tsv(file.path(res, "06_twisst", "twisst_summary.tsv"), show_col_types = FALSE)
names(gw) <- str_replace(names(gw), "^(T[123])_.*", "\\1")
gw_sh <- to_shares(gw, outgroup) |>
  mutate(outgroup = factor(outgroup, levels = rev(outgroups)), panel = default_lab)
rob <- read_tsv(file.path(res, "07_twisst_robustness", "robustness_summary.tsv"), show_col_types = FALSE)
rob_sh <- to_shares(rob, representation, window_snps, method) |>
  mutate(representation = factor(rep_lab[representation], levels = rep_lab),
         window = factor(paste(window_snps, "SNPs"), levels = paste(sort(unique(window_snps)), "SNPs")),
         method = factor(method, levels = c("ML", "NJ")))                 # NJ on top
write_tsv(bind_rows(mutate(gw_sh, run = "default"), mutate(rob_sh, run = "robustness")) |>
            select(run, outgroup, representation, window_snps, method, topo, weight = w, share),
          paste0(prefix, ".mean_weights.tsv"))

bar_layers <- function() list(
  geom_col(aes(fill = topo), width = 0.8, colour = "white", linewidth = 0.3, position = position_stack(reverse = TRUE)),
  geom_text(aes(label = round(100 * share), colour = topo), size = 6 / .pt, family = fig_font,
            position = position_stack(vjust = 0.5, reverse = TRUE)),
  geom_vline(xintercept = c(1, 2) / 3, linetype = "22", colour = ink, linewidth = 0.25),
  scale_fill_topo(),
  # percentages: white on blue and red, ink on grey
  scale_colour_manual(values = c(T1 = "white", T3 = ink, T2 = "white"), guide = "none"),
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c(0, 50, 100), expand = c(0, 0)),
  theme_fig(),
  theme(panel.grid = element_blank(), strip.text = element_text(face = "plain", hjust = 0.5)))

pOutBar <- ggplot(gw_sh, aes(x = share, y = outgroup)) + bar_layers() +
  facet_wrap(~ panel) +
  scale_y_discrete(labels = function(x) sample_expr(x)) +
  labs(x = NULL, y = NULL)
pGridBar <- ggplot(rob_sh, aes(x = share, y = method)) + bar_layers() +
  facet_grid(window ~ representation, switch = "y") +
  labs(x = "Share of topology weight (%)", y = NULL) +
  theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0, hjust = 1),
        panel.spacing.y = unit(1.2, "mm"), panel.spacing.x = unit(5, "mm"))
# relative heights of the two parts: number of bars (3 outgroups; 4 window sizes x 2 methods) plus strips
bars_heights <- c(3 + 1.2, 8 + 2.6)

# the two parts stacked; the topology key collected (`legend_theme`: its position; `extra`: text size)
mean_weight_bars <- function(p_out, p_grid, heights, legend_theme, extra = theme())
  (p_out / p_grid) + plot_layout(heights = heights, guides = "collect") & legend_theme & extra
# key below, title on top, entries in one row
key_bottom <- theme(legend.position = "bottom", legend.justification = "left", legend.title.position = "top",
                    legend.key.size = unit(3, "mm"), legend.margin = margin(0, 0, 0, 0),
                    legend.box.spacing = unit(1.5, "mm"))

panels <- list(along = pAlong, along_vertical = pAlongV, outgroup_bars = pOutBar, grid_bars = pGridBar,
               bars_heights = bars_heights)
saveRDS(panels, paste0(prefix, ".panels.rds"))

fig <- (pAlong / wrap_elements(full = mean_weight_bars(pOutBar, pGridBar, bars_heights, key_bottom))) +
  plot_layout(heights = c(1.3, 3.6)) +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 11, face = "bold"))
save_fig(fig, prefix, width = 7.5, height = 5.1)
