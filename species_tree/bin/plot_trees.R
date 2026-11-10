#!/usr/bin/env Rscript
# Plots of the earlier plot_astral_res.R, run on one set's trees:
#   1. ASTRAL species tree as a cladogram with q1/q2/q3 quartet-support pies
#   2. tanglegram ASTRAL vs concatenated ML tree, differing clades in red, normalised RF distance
# usage: plot_trees.R astral.tre concat.treefile outgroup prefix
suppressPackageStartupMessages({ library(ape); library(phytools); library(phangorn) })

args <- commandArgs(trailingOnly = TRUE)
astral_file <- args[1]; concat_file <- args[2]; outgroup <- args[3]; prefix <- args[4]

col_q <- c("#2a78d6", "#eb6834", "#1baf7a")   # q1, q2, q3
col_diff <- "#c8312b"; col_new <- "#6a3d9a"; col_note <- "grey35"

# "Plobata1_GCA_973946955.1" -> "Porites lobata 1 (GCA_973946955.1)" (acc = FALSE for the tanglegram, where
# phytools would turn "GCA_" into "GCA "); a tip without an accession is the new assembly
# ("Porites_sp_nov_indonesia" -> "Porites sp. Indonesia")
genera <- c(P = "Porites", A = "Acropora")
pretty_tip <- function(tip, acc = TRUE) {
  m <- regmatches(tip, regexec("^([A-Z])([a-z]+)([0-9]*)_(GC[AF]_[0-9.]+)$", tip))[[1]]
  if (length(m)) return(sprintf("%s %s%s%s", genera[[m[2]]], m[3], if (m[4] == "") "" else paste0(" ", m[4]),
                                if (acc) sprintf(" (%s)", m[5]) else ""))
  w <- strsplit(tip, "_")[[1]]
  rest <- w[-1][!tolower(w[-1]) %in% c("sp", "nov")]
  sprintf("%s sp. %s", w[1], paste0(toupper(substring(rest, 1, 1)), substring(rest, 2), collapse = " "))
}
is_new <- function(tip) !grepl("_GC[AF]_", tip)
qval <- function(label, key) {
  m <- regmatches(label, regexpr(paste0("\\b", key, "=[0-9.eE+-]+"), label))
  if (length(m) == 0) NA_real_ else as.numeric(sub(".*=", "", m))
}
clades <- function(t) sapply(prop.part(t), function(i) paste(sort(t$tip.label[i]), collapse = ","))
# supports are read off the trees as written and keyed by clade: root() on the already rooted ASTRAL
# tree shifts node labels one node rootwards, and cophylo rotates nodes
supports <- function(t) setNames(c(t$node.label, rep("", Nnode(t)))[seq_len(Nnode(t))], clades(t))
astral_raw <- read.tree(astral_file); concat_raw <- read.tree(concat_file)
sup_astral <- supports(astral_raw); sup_concat <- supports(concat_raw)
reroot <- function(t) { t <- root(t, outgroup, resolve.root = TRUE); t$node.label <- NULL; t }
# every plot is written as PDF and as 300 dpi PNG
save_plot <- function(stem, width, height, draw) {
  cairo_pdf(paste0(stem, ".pdf"), width = width, height = height); draw(); invisible(dev.off())
  png(paste0(stem, ".png"), width = width, height = height, units = "in", res = 300, type = "cairo"); draw(); invisible(dev.off())
}

## 1. quartet support
tr <- ladderize(reroot(astral_raw), right = FALSE)
tr$node.label <- unname(sup_astral[clades(tr)])
q_mat <- t(sapply(tr$node.label, function(l) c(q1 = qval(l, "q1"), q2 = qval(l, "q2"), q3 = qval(l, "q3"))))
pp <- sapply(tr$node.label, qval, key = "localPP")
has_data <- stats::complete.cases(q_mat)
node_ids <- ((Ntip(tr) + 1):(Ntip(tr) + Nnode(tr)))[has_data]
q_lab <- ifelse(pp < 0.995, sprintf("PP %.2f\nq1 %.2f", pp, q_mat[, "q1"]), sprintf("q1 %.2f", q_mat[, "q1"]))[has_data]
clade_ids <- clades(tr)
new_tip <- is_new(tr$tip.label)
tr$tip.label <- sapply(tr$tip.label, pretty_tip)

save_plot(paste0(prefix, "_astral_quartets"), width = 12, height = 0.5 * Ntip(tr) + 1, function() {
  par(mar = c(2.4, 0.5, 0.5, 0.5))
  plot(tr, use.edge.length = FALSE, cex = 0.95, label.offset = 0.15, edge.width = 1.8, font = 3, underscore = TRUE,
       tip.color = ifelse(new_tip, col_new, "black"), align.tip.label = 3)   # 3 = dotted leaders
  nodelabels(node = node_ids, pie = q_mat[has_data, ], piecol = col_q, cex = 0.8)
  # labels sit above the branch entering each node, ending just left of the pie (ape's pie radius)
  lp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  r <- 0.8 * diff(par("usr")[1:2]) / 50
  text(lp$xx[node_ids] - 1.6 * r, lp$yy[node_ids] + 0.12, q_lab, adj = c(1, 0), cex = 0.82, col = "grey15")
  legend("bottom", inset = c(0, -0.3 / par("pin")[2]), xpd = NA, horiz = TRUE, bty = "n", cex = 0.85, pch = 22,
         pt.cex = 1.7, pt.bg = col_q, legend = c("Main topology (q1)", "First alternative (q2)", "Second alternative (q3)"))
})
write.table(data.frame(node = node_ids, clade = clade_ids[has_data], round(q_mat[has_data, , drop = FALSE], 4)),
            paste0(prefix, "_astral_quartets.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## 2. tanglegram and RF distance
tr1 <- reroot(astral_raw)
tr2 <- reroot(concat_raw)
tr1$edge.length <- NULL
tr2$edge.length <- NULL
rf <- RF.dist(tr1, tr2, normalize = TRUE)
rf_splits <- RF.dist(tr1, tr2)
writeLines(sprintf("normalised_RF\t%s", rf), paste0(prefix, "_rf.tsv"))

obj <- cophylo(tr1, tr2, rotate = TRUE, print = FALSE)
r1 <- obj$trees[[1]]; r2 <- obj$trees[[2]]
c1 <- clades(r1); c2 <- clades(r2)
n1 <- (Ntip(r1) + 1):(Ntip(r1) + Nnode(r1)); n2 <- (Ntip(r2) + 1):(Ntip(r2) + Nnode(r2))
d1 <- n1[!(c1 %in% c2)]; d2 <- n2[!(c2 %in% c1)]
s1 <- sapply(sup_astral[c1], qval, key = "localPP"); s2 <- sup_concat[c2]
# supports are printed on every branch that differs between the trees and on any branch below maximum
lab1 <- which(n1 %in% d1 | (!is.na(s1) & s1 < 0.995)); lab2 <- which(n2 %in% d2 | (nzchar(s2) & s2 != "100/100"))
lab2 <- lab2[!is.na(s2[lab2]) & nzchar(s2[lab2])]
for (i in 1:2) obj$trees[[i]]$tip.label <- sapply(obj$trees[[i]]$tip.label, pretty_tip, acc = FALSE)
obj$assoc <- apply(obj$assoc, 2, function(x) sapply(x, pretty_tip, acc = FALSE))
n <- Ntip(r1)
part <- 0.45   # share of the width for each tree and its labels

save_plot(paste0(prefix, "_tanglegram"), width = 14, height = 0.5 * n + 1.2, function() {
  plot(obj, part = part, fsize = 0.9, ftype = "i", pts = FALSE, lwd = 2, mar = c(2.2, 0.5, 1.6, 0.5),
       link.type = "curved", link.lty = "solid", link.lwd = 1, link.col = make.transparent("grey40", 0.6),
       edge.col = list(left = ifelse(r1$edge[, 2] %in% d1, col_diff, "black"),
                       right = ifelse(r2$edge[, 2] %in% d2, col_diff, "black")))
  # values sit at node height between the two child branches, where cladograms leave free space
  if (length(lab1)) nodelabels.cophylo(sprintf("PP %.2f", s1[lab1]), n1[lab1], which = "left",
                                       frame = "none", adj = c(-0.12, 0.5), cex = 0.82, col = "grey15")
  if (length(lab2)) nodelabels.cophylo(s2[lab2], n2[lab2], which = "right",
                                       frame = "none", adj = c(1.12, 0.5), cex = 0.82, col = "grey15")
  # phytools has no per-tip label colour: redraw the new assembly's name in purple where phylogram()
  # puts it (labels end at x = +/-(part - 0.5), on the side facing the links)
  lpc <- get("last_plot.cophylo", envir = .PlotPhyloEnv)
  for (k in 1:2) {
    tip <- which(is_new(c(r1, r2)[[k]]$tip.label))
    text(c(1, -1)[k] * (part - 0.5), lpc[[k]]$yy[tip], obj$trees[[k]]$tip.label[tip], pos = c(2, 4)[k],
         offset = 0, cex = 0.9, font = 3, col = col_new)
  }
  mtext("ASTRAL-IV", side = 3, line = 0.2, adj = 0.02, font = 2, cex = 0.95)
  mtext("Concatenated ML", side = 3, line = 0.2, adj = 0.98, font = 2, cex = 0.95)
  legend("bottom", inset = c(0, -0.3 / par("pin")[2]), xpd = NA, bty = "n", cex = 0.85, lwd = 2,
         col = col_diff, seg.len = 2, legend = "Discordant branch")
})

## notes: what the figures show, for captions
short <- function(key) paste(sapply(strsplit(key, ",")[[1]], pretty_tip, acc = FALSE), collapse = ", ")
diff_rows <- c(sprintf("| ASTRAL-IV | %s | local PP %.2f |", sapply(c1[n1 %in% d1], short), s1[n1 %in% d1]),
               sprintf("| Concatenated ML | %s | SH-aLRT/UFBoot %s |", sapply(c2[n2 %in% d2], short), s2[n2 %in% d2]))
writeLines(c(
  sprintf("# Species-tree plots, %s set", prefix), "",
  sprintf("%d genomes, rooted on %s. The new assembly is %s, shown in purple.",
          n, pretty_tip(outgroup), pretty_tip(tr1$tip.label[is_new(tr1$tip.label)])), "",
  sprintf("## %s_astral_quartets (.pdf, .png)", prefix), "",
  "ASTRAL-IV species tree from the BUSCO single-copy gene trees, drawn as a cladogram (no branch lengths).",
  "Pies: share of gene-tree quartets supporting the main topology of each branch (q1, blue) and the two",
  "alternatives (q2, orange; q3, green). Each node is labelled with q1; the local posterior probability (PP)",
  sprintf("is added where it is below 1. Values per clade are in %s_astral_quartets.tsv.", prefix), "",
  sprintf("## %s_tanglegram (.pdf, .png)", prefix), "",
  "ASTRAL-IV species tree (left) against the partitioned concatenated ML tree from IQ-TREE (right), both as",
  "cladograms, with nodes rotated to minimise crossing links. Red: discordant branches, i.e. branches leading",
  "to clades absent from the other tree. Support values are given on every red branch and on any branch below maximal support:",
  "local PP on the left, SH-aLRT (%) / ultrafast bootstrap (%) on the right; all other branches have",
  "local PP >= 0.99 and 100/100.", "",
  sprintf("Normalised Robinson-Foulds distance: %.2f (%d of %d splits differ; also in %s_rf.tsv).",
          rf, as.integer(rf_splits), 2L * (n - 3L), prefix), "",
  "Clades present in one tree only:", "",
  "| Tree | Clade | Support |", "|---|---|---|", diff_rows),
  paste0(prefix, "_plots_notes.md"))
