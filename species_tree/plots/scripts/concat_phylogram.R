#!/usr/bin/env Rscript
# Concatenated ML tree (IQ-TREE, partitioned) with branch lengths, rooted on the outgroup. The path between
# the Porites crown and the outgroup is drawn shortened, with a break mark; its real length is in
# <prefix>.tsv. Supports (SH-aLRT / UFBoot) shown only where below 100/100. Scale bar in substitutions/site.
# Writes <prefix>.pdf/.png, <prefix>.tsv and <prefix>.rds (for the composite).
# usage: concat_phylogram.R <concat .treefile> <outgroup tip ID> <out prefix> [tree width, in; default 3.4]
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(script_dir, "../../../shared/plots/setup.R"))   # names, colours, shared functions
suppressPackageStartupMessages(library(ape))

args <- commandArgs(trailingOnly = TRUE)
tree_file <- args[1]; outgroup <- args[2]; prefix <- args[3]
tree_w <- if (length(args) >= 4) as.numeric(args[4]) else 3.4

tr <- shorten_outgroup(read_tree_rooted(tree_file, outgroup), outgroup)
d <- tree_layout(tr)
fig <- tree_plot(d, tree_w)
upi <- attr(fig, "upi"); h <- tree_height(Ntip(tr))

## break mark halfway along the outgroup branch
og <- d |> filter(label == outgroup)
fig <- fig + break_mark((d$x[d$node == og$parent] + og$x) / 2, og$y, upi, 1 / tree_tip_h)

## supports below maximal
sup <- d |> filter(!isTip, nzchar(label), label != "100/100")
fig <- fig + geom_text(data = sup, aes(x, y, label = label), hjust = 1.08, vjust = -0.35,
                       size = node_pt / .pt, colour = ink_muted)

## scale bar: top left inside the panel, where the tree is empty; for a narrow tree (< 2 in) the top left is
## too small for its label, so it goes right-aligned between the two bottom rows, label under the bar (the
## outgroup label is short, the row above is not)
sb <- 0.005
narrow <- tree_w < 2
sb_y <- if (narrow) 1.6 else Ntip(tr) - 0.3
sb_x <- if (narrow) max(layer_scales(fig)$x$limits) - sb - 0.02 * upi else 0
fig <- fig +
  annotate("segment", x = sb_x, xend = sb_x + sb, y = sb_y, yend = sb_y, colour = ink, linewidth = 0.45) +
  annotate("text", x = if (narrow) sb_x + sb else sb_x, y = sb_y, label = paste(sb, "substitutions/site"),
           hjust = if (narrow) 1 else 0, vjust = if (narrow) 1.5 else -0.9, size = node_pt / .pt, colour = ink)

w <- attr(fig, "fig_w")
save_fig(fig, prefix, w, h)
saveRDS(list(plot = fig, width = w, height = h), paste0(prefix, ".rds"))
write_tsv(tibble(outgroup = tree_tip_label(outgroup), crown_to_outgroup_length = attr(tr, "crown_to_outgroup"),
                 drawn_length = attr(tr, "drawn"), max_ingroup_root_to_tip = attr(tr, "ingroup_depth"),
                 scale_bar = sb), paste0(prefix, ".tsv"))
