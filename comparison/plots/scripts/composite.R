#!/usr/bin/env Rscript
# Composite figure of the funannotate compare analysis, sized for a portrait US Letter page
# (7.5 in wide, <= 10 in tall; heights follow the content), in two versions that differ only by the table:
#   top:          A. orthogroup sharing (UpSet, largest combinations), full width
#   bottom left:  B. GO enrichment dot matrix (go_enrichment panel B)
#   bottom right: C. COG composition (new genome against the other Porites), then (<prefix>_with_table
#                 only) D. the assembly/annotation table. GO legend: below B with the table, right of B without.
# The UpSet and GO panels are the plot objects saved by orthology_upset.R and go_enrichment.R (.rds).
# usage: composite.R <funannotate_compare dir> <out prefix> <orthology .rds> <go .panels.rds> <stats table .tsv>
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(script_dir, "..", "..", "..", "shared", "plots", "setup.R"))   # names, colours, genome order, shared functions
suppressPackageStartupMessages(library(gridExtra))

args <- commandArgs(trailingOnly = TRUE)
fc_dir <- args[1]; prefix <- args[2]; ortho_rds <- args[3]; go_rds <- args[4]; stats_tsv <- args[5]

page_w <- 7.5          # in; portrait US Letter with margins
base   <- 7            # pt; text size across the composite

upset <- readRDS(ortho_rds)
go <- readRDS(go_rds)

## GO dot matrix at its natural size. Its legend: one line centred below the whole panel, labels included
## (version with the table); one column to the right of the panel (version without)
pGO <- go$B + theme(legend.title = element_text(size = base - 0.5),
                    legend.text = element_text(size = base - 1),
                    legend.key.height = unit(2.5, "mm"),
                    axis.text.y = element_text(size = base - 1, lineheight = 0.85))
pGO_below <- pGO + theme(legend.position = "bottom", legend.location = "plot", legend.justification = "center",
                         legend.direction = "horizontal", legend.box = "horizontal",
                         legend.box.just = "bottom", legend.title.position = "top",
                         legend.key.width = unit(4, "mm"))
pGO_right <- pGO + theme(legend.position = "right", legend.direction = "vertical", legend.box = "vertical",
                         legend.box.just = "left", legend.key.height = unit(3.5, "mm"),
                         legend.key.width = unit(2.5, "mm"))

## COG composition: % of COG-assigned proteins; grey bar = range of the other Porites
cog <- read_fc_table(file.path(fc_dir, "cogs", "COGS.all.results.csv")) |>
  select(id, all_of(porites)) |>
  pivot_longer(-id, names_to = "genome", values_to = "n") |>
  group_by(genome) |> mutate(pct = 100 * n / sum(n)) |> ungroup() |>
  mutate(label = str_replace(id, "^\\((\\w)\\) ", "\\1  "), is_new = genome == new_genome)
# COG panel; `wrap` = label width in characters (narrower labels leave more room for the plot)
cog_panel <- function(wrap) {
  d <- cog |> mutate(category = str_wrap(label, wrap),
                     category = fct_rev(factor(category, levels = sort(unique(category)))))
  rng <- d |> filter(!is_new) |> group_by(category) |> summarise(lo = min(pct), hi = max(pct))
  p <- ggplot(d, aes(x = pct, y = category)) +
    geom_segment(data = rng, aes(x = lo, xend = hi, y = category, yend = category),
                 inherit.aes = FALSE, colour = grid_col, linewidth = 1.8, lineend = "round") +
    geom_point(data = filter(d, !is_new), colour = col_other, size = 1) +
    geom_point(data = filter(d, is_new), colour = col_new, size = 1.8) +
    scale_x_continuous(breaks = c(0, 10, 20)) +
    labs(x = "Proteins (%)", y = NULL) +
    theme_fig(base) + theme(axis.text.y = element_text(size = base - 1.5, lineheight = 0.8),
                            panel.grid.major.y = element_blank())
  # not aligned with the GO panel: keeps its own x-axis title close to the axis
  wrap_elements(full = p)
}

## assembly/annotation table, genomes as rows (same order as the UpSet matrix), new genome in the accent
cols <- c("Assembly size (Mb)" = "Assembly\n(Mb)", "Scaffold N50 (Mb)" = "N50\n(Mb)",
          "Protein-coding genes" = "Genes", "Proteins with an ortholog" = "Proteins\nw/ ortholog",
          "Proteins without an ortholog" = "Proteins\nw/o ortholog")
digits <- c("Scaffold N50 (Mb)" = 1)                     # decimals per metric (default 0)
st <- read_tsv(stats_tsv, show_col_types = FALSE) |> filter(Metric %in% names(cols))
m <- t(as.matrix(st[, porites]))
dig <- ifelse(is.na(digits[st$Metric]), 0, digits[st$Metric])
m[] <- mapply(function(v, d) formatC(v, format = "f", digits = d, big.mark = ","), m, rep(dig, each = nrow(m)))
colnames(m) <- cols[st$Metric]; m <- m[, cols]
rownames(m) <- vapply(genome_lab[porites], plotmath_label, "")
tt <- ttheme_minimal(
  base_size = base, base_family = fig_font,
  core    = list(fg_params = list(hjust = 1, x = 0.9, col = ink),
                 bg_params = list(fill = c("white", band), col = NA)),
  colhead = list(fg_params = list(hjust = 1, x = 0.9, col = ink_muted, fontface = "plain", lineheight = 0.9)),
  rowhead = list(fg_params = list(parse = TRUE, hjust = 0, x = 0.05, col = ink, fontface = "plain")),
  padding = unit(c(1.6, 1.8), "mm"))
tab <- tableGrob(m, theme = tt)
for (i in which(tab$layout$t == 2 & tab$layout$name %in% c("core-fg", "rowhead-fg")))
  tab$grobs[[i]]$gp$col <- col_new
pTAB <- wrap_elements(full = tab)

# text size applied before wrapping: a wrapped panel does not receive the composite-wide theme
pUP <- wrap_elements(full = upset & theme(text = element_text(size = base),
                                          axis.title = element_text(size = base)))

# one text size for the whole composite (panel-specific sizes set above stay)
small <- theme(text = element_text(size = base), legend.text = element_text(size = base - 1),
               plot.tag = element_text(size = base + 3, face = "bold"))

# heights come from the content (in per row), not from the page: the figure only has to fit within
# 7.5 x 10 in, not fill it. No panel is stretched: the shorter column gets empty space below it.
h_top <- 2.8
h_go  <- go$n_rows_B * 0.16 + 0.9               # GO rows + genome labels
h_leg   <- 0.5                                   # GO legend: one line below the panel
h_cog <- function(row) n_distinct(cog$label) * row + 0.5   # COG height for a given row height
h_tab <- 1.6
# a column of panels with fixed heights, padded with empty space to height H
column <- function(plots, heights, H) {
  if (H > sum(heights) + 0.05) { plots <- c(plots, list(plot_spacer())); heights <- c(heights, H - sum(heights)) }
  wrap_plots(plots, ncol = 1, heights = heights)
}

make <- function(with_table) {
  if (with_table) {
    H <- max(h_go + h_leg, h_cog(0.17) + h_tab)
    left  <- column(list(pGO_below), h_go + h_leg, H)
    right <- column(list(cog_panel(46), pTAB), c(h_cog(0.17), h_tab), H)
    widths <- c(0.45, 1)
  } else {
    H <- max(h_go, h_cog(0.21))   # narrower, two-line COG labels need taller rows
    left  <- column(list(pGO_right), h_go, H)
    right <- column(list(cog_panel(30)), h_cog(0.21), H)
    widths <- c(0.5, 1)
  }
  stopifnot(h_top + H + 0.2 <= 10)
  fig <- (pUP / ((left | right) + plot_layout(widths = widths))) + plot_layout(heights = c(h_top, H)) +
    plot_annotation(tag_levels = "A") & small
  list(fig = fig, height = h_top + H + 0.2)
}

v1 <- make(FALSE)
save_fig(v1$fig, prefix, width = page_w, height = v1$height)
v2 <- make(TRUE)
save_fig(v2$fig, paste0(prefix, "_with_table"), width = page_w, height = v2$height)
