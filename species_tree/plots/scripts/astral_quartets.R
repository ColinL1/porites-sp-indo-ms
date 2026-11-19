#!/usr/bin/env Rscript
# ASTRAL-IV species tree as a cladogram, with the gene-tree quartet support of each branch as a pie
# (q1 main topology, q2 and q3 alternatives), q1 above each branch and the local PP where it is below 1.
# Writes <prefix>.pdf/.png, <prefix>.tsv (values per clade) and <prefix>.rds (plot and size, for the composite).
# With "pp_only" (for a narrow tree) the q1 numbers are left out, only PP < 1 is written (above the branch
# into the node, ending at the pie) and the legend is a compact q1/q2/q3 key.
# usage: astral_quartets.R <astral .tre> <outgroup tip ID> <out prefix> [tree width, in; default 3.4] [pp_only]
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(script_dir, "../../../shared/plots/setup.R"))   # names, colours, shared functions
suppressPackageStartupMessages(library(ape))

args <- commandArgs(trailingOnly = TRUE)
tree_file <- args[1]; outgroup <- args[2]; prefix <- args[3]
tree_w <- if (length(args) >= 4) as.numeric(args[4]) else 3.4
pp_only <- length(args) >= 5 && args[5] == "pp_only"

tr <- read_tree_rooted(tree_file, outgroup)
d <- tree_layout(tr, branch.length = "none")
nodes <- d |>
  filter(!isTip) |>
  mutate(q1 = astral_value(label, "q1"), q2 = astral_value(label, "q2"), q3 = astral_value(label, "q3"),
         pp = astral_value(label, "localPP"), x_parent = d$x[match(parent, d$node)]) |>
  filter(!is.na(q1))
nodes <- nodes |>
  mutate(lab = ifelse(pp < 0.995, sprintf("PP %.2f\nq1 %.2f", pp, q1), sprintf("q1 %.2f", q1)))

## sizes (in): tree width fixed, label column as wide as the widest tip label
r_mm <- 1.7
fig <- tree_plot(d, tree_w)
upi <- attr(fig, "upi")
fig <- fig + if (pp_only) {
  geom_text(data = filter(nodes, pp < 0.995) |> mutate(xs = x - (r_mm / 25.4 + 0.02) * upi),
            aes(xs, y, label = sprintf("PP %.2f", pp)), hjust = 1, vjust = -0.35, size = node_pt / .pt,
            colour = ink_muted)
} else {
  geom_text(data = nodes, aes((x_parent + x) / 2, y, label = lab), vjust = -0.35, size = node_pt / .pt,
            colour = ink_muted, lineheight = 0.9)
}
# pies: one grob per node, a fixed physical radius so they stay round whatever the axis scales
for (i in seq_len(nrow(nodes)))
  fig <- fig + annotation_custom(pie_grob(unlist(nodes[i, c("q1", "q2", "q3")]), r_mm = r_mm),
                                 xmin = nodes$x[i] - 1, xmax = nodes$x[i] + 1,
                                 ymin = nodes$y[i] - 1, ymax = nodes$y[i] + 1)
# legend: a dummy layer gives the fill key; placed inside the panel, top left, where the cladogram is empty
fig <- fig +
  geom_point(data = tibble(q = factor(names(col_quartet), names(col_quartet)), x = 0, y = 1),
             aes(x, y, fill = q), shape = 21, size = 0, stroke = 0, inherit.aes = FALSE) +
  # narrow tree: compact key (q1-q3 spelled out in the legend text of the figure)
  scale_fill_manual(values = col_quartet, name = if (pp_only) "Quartets" else "Gene-tree quartets",
                    labels = if (pp_only) names(col_quartet) else
                      c("Main topology (q1)", "First alternative (q2)", "Second alternative (q3)")) +
  guides(fill = guide_legend(override.aes = list(size = 3.2, stroke = 0.3))) +
  theme(legend.position = "inside", legend.position.inside = c(0, 1), legend.justification = c(0, 1),
        legend.title = element_text(size = tip_pt), legend.text = element_text(size = tip_pt - 0.5),
        legend.key.size = unit(3.6, "mm"))

h <- tree_height(Ntip(tr))
w <- attr(tree_plot(d, tree_w), "fig_w")
save_fig(fig, prefix, w, h)
saveRDS(list(plot = fig, width = w, height = h), paste0(prefix, ".rds"))
nodes |>
  transmute(clade = tree_clades(tr, outgroup)[node - Ntip(tr)], localPP = pp, q1, q2, q3) |>
  write_tsv(paste0(prefix, ".tsv"))
