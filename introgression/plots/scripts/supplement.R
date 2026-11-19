#!/usr/bin/env Rscript
# Supplementary composite, sized for a portrait US Letter page (7.5 in wide; height from the content):
#   A. ADMIXTURE, K = 2 to 5 | B. mitogenome tree
# The PCA is left out for now (pca_admixture.R still draws it in its own figure).
# usage: supplement.R <out prefix> <pca_admixture .panels.rds> <mito .rds>
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(script_dir, "setup_introgression.R"))

args <- commandArgs(trailingOnly = TRUE)
prefix <- args[1]
pa <- readRDS(args[2]); pTree <- readRDS(args[3])

base <- 7
small <- theme(text = element_text(size = base), axis.text = element_text(size = base - 1),
               strip.text = element_text(size = base),
               plot.tag = element_text(size = base + 3, face = "bold"))

fig <- (pa$admixture | pTree) + plot_layout(widths = c(1.2, 1)) +
  plot_annotation(tag_levels = "A") & small
save_fig(fig, prefix, width = 7.5, height = 2.4)
