#!/usr/bin/env Rscript
# Mitogenome ML tree (IQ-TREE, 10_mitogenome/iqtree/mito_tree.treefile), rooted on P. divaricata
# (the pipeline's mito_root). Branches with ultrafast bootstrap < min_support are collapsed into
# polytomies; the remaining support values are printed. Drawn as a cladogram: the P. divaricata branch
# is ~20 times longer than any other, so branch lengths to scale would compress the rest to a point.
# usage: mitogenome.R <results dir> <out prefix>
# writes <prefix>.pdf/.png, <prefix>.rds (panel, for the composite)
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(script_dir, "setup_introgression.R"))
suppressPackageStartupMessages({ library(ape); library(ggtree) })

args <- commandArgs(trailingOnly = TRUE)
res <- args[1]; prefix <- args[2]
root_on <- "Pdivaricata2"
min_support <- 50

tr <- read.tree(file.path(res, "10_mitogenome", "iqtree", "mito_tree.treefile"))
tr <- root(tr, outgroup = root_on, resolve.root = TRUE, edgelabel = TRUE)   # support stays with its branch
# collapse weakly supported internal branches (support is the node label of the branch's child node)
sup <- c(rep(NA, Ntip(tr)), suppressWarnings(as.numeric(tr$node.label)))   # indexed by node number
weak <- which(!is.na(sup[tr$edge[, 2]]) & sup[tr$edge[, 2]] < min_support)
tr$edge.length[weak] <- 0
tr <- di2multi(tr, tol = 1e-12)
tr$edge.length <- NULL

tr$node.label[is.na(suppressWarnings(as.numeric(tr$node.label)))] <- ""   # no label on the root
tips <- tibble(label = tr$tip.label, lab = vapply(sample_lab[tr$tip.label], plotmath_label, ""),
               is_new = tr$tip.label == target)
p <- ggtree(tr, colour = ink, linewidth = 0.35, ladderize = TRUE) %<+% tips +
  geom_tiplab(aes(label = lab, colour = is_new), parse = TRUE, size = 2.4, offset = 0.3, family = fig_font) +
  geom_tippoint(aes(colour = is_new), size = 1.2) +
  geom_nodelab(aes(label = label), size = 2.2, family = fig_font, colour = ink_muted, hjust = 1.25, vjust = -0.5) +
  scale_colour_manual(values = c(`TRUE` = col_new, `FALSE` = ink), guide = "none") +
  coord_cartesian(clip = "off") +
  theme_tree(bgcolor = "white") + theme(text = element_text(family = fig_font), plot.margin = margin(4, 70, 4, 4))

saveRDS(p, paste0(prefix, ".rds"))
save_fig(p, prefix, width = 3, height = 2.2)
