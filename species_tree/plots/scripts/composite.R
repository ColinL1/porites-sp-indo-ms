#!/usr/bin/env Rscript
# Composite species tree figure for a portrait US Letter page (<= 7.5 x 10 in), panels at their natural size:
#   A. full width (tanglegram)
#   B. below A, left-aligned (not stretched; empty space to its right), or, with a third panel,
#   B and C side by side below A
# Panel letters sit in a narrow gutter left of each panel, level with its top.
# usage: composite.R <out prefix> <A .rds> <B .rds> [C .rds]
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(script_dir, "../../../shared/plots/setup.R"))   # names, colours, shared functions
suppressPackageStartupMessages(library(cowplot))

args <- commandArgs(trailingOnly = TRUE)
prefix <- args[1]
p <- lapply(args[-1], readRDS)
names(p) <- LETTERS[seq_along(p)]

gutter <- 0.2; gap <- 0.15                                 # in
bottom <- p[-1]                                            # panels of the bottom row, left to right
bottom_w <- sum(vapply(bottom, function(x) gutter + x$width, 0))
bottom_h <- max(vapply(bottom, function(x) x$height, 0))
W <- max(gutter + p$A$width, bottom_w)
H <- p$A$height + gap + bottom_h
stopifnot(W <= 7.5, H <= 10)

# cowplot positions are fractions of the page, from the bottom left; each panel top-aligned in its row
place <- function(fig, panel, letter, x, top) {
  fig +
    draw_plot(panel$plot, x = (x + gutter) / W, y = (top - panel$height) / H,
              width = panel$width / W, height = panel$height / H) +
    draw_label(letter, x = (x + 0.02) / W, y = (top - 0.04) / H, hjust = 0, vjust = 1, fontface = "bold",
               size = 10, colour = ink)
}
fig <- place(ggdraw(), p$A, "A", 0, H)
x <- 0
for (k in names(bottom)) {
  fig <- place(fig, p[[k]], k, x, bottom_h)
  x <- x + gutter + p[[k]]$width
}
fig <- fig + theme(plot.background = element_rect(fill = "white", colour = NA))

save_fig(fig, prefix, W, H)
cat(sprintf("%s: %.2f x %.2f in\n", prefix, W, H))
