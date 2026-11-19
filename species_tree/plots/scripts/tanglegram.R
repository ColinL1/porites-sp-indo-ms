#!/usr/bin/env Rscript
# Tanglegram: ASTRAL-IV species tree (left) against the concatenated ML tree (right), both as cladograms,
# the right tree rotated to follow the left tip order. Discordant branches (leading to a clade absent from
# the other tree) in col_discord; supports on discordant branches and on any branch below maximal support.
# Writes <prefix>.pdf/.png, <prefix>.discordant.tsv, <prefix>.rf.tsv and <prefix>.rds (for the composite).
# Options (any order, after the prefix):
#   phylogram  right tree drawn with branch lengths (dotted leaders to the labels, outgroup path shortened
#              with a break mark, scale bar bottom right)
#   quartets   left tree with the ASTRAL gene-tree quartet support of each branch as pies (q1 main topology,
#              q2/q3 alternatives); their key joins the discordant-branch key in the header row
# usage: tanglegram.R <astral .tre> <concat .treefile> <outgroup tip ID> <out prefix> [phylogram] [quartets]
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(script_dir, "../../../shared/plots/setup.R"))   # names, colours, shared functions
suppressPackageStartupMessages(library(ape))

args <- commandArgs(trailingOnly = TRUE)
astral_file <- args[1]; concat_file <- args[2]; outgroup <- args[3]; prefix <- args[4]
right_phylo <- "phylogram" %in% args[-(1:4)]
left_pies <- "quartets" %in% args[-(1:4)]
r_mm <- 1.6                                                # pie radius

t1 <- read_tree_rooted(astral_file, outgroup)
t2 <- read_tree_rooted(concat_file, outgroup)
d1 <- tree_layout(t1, branch.length = "none")
t2 <- rotateConstr(   # tip order of the left tree, bottom to top
  t2, t1$tip.label[order(d1$y[d1$isTip])])
if (right_phylo) t2 <- shorten_outgroup(t2, outgroup)
d2 <- if (right_phylo) tree_layout(t2) else tree_layout(t2, branch.length = "none")

## discordance and supports, per node of each layout
node_info <- function(t, d, other) {
  key <- tree_clades(t, outgroup)
  d |> mutate(clade = ifelse(isTip, NA, key[pmax(node - Ntip(t), 1)]),
              discordant = !isTip & node != parent & !is.na(clade) & !(clade %in% other) &
                           nchar(clade) > 0 & node != (Ntip(t) + 1))
}
k1 <- tree_clades(t1, outgroup); k2 <- tree_clades(t2, outgroup)
d1 <- node_info(t1, d1, k2) |> mutate(pp = astral_value(label, "localPP"),
                                      sup = ifelse(!isTip & !is.na(pp), sprintf("PP %.2f", pp), ""),
                                      show = !isTip & (discordant | (!is.na(pp) & pp < 0.995)))
d2 <- node_info(t2, d2, k1) |> mutate(sup = ifelse(isTip, "", label),
                                      show = !isTip & nzchar(sup) & (discordant | sup != "100/100"))

## horizontal layout (in): tree | labels | links | labels | tree
# tree widths: equal for two cladograms; a right phylogram gets more room for its short internal branches,
# taken from the left tree unless that carries pies
tw1 <- if (right_phylo && !left_pies) 1.0 else 1.35; tw2 <- if (right_phylo) 1.7 else 1.35
off <- 0.06; link_w <- 0.45
lab_w <- label_width_in(tree_tip_label(t1$tip.label), tip_pt)
W <- tw1 + tw2 + 2 * (off + lab_w) + link_w
to_left  <- function(x, d) x / max(d$x) * tw1
to_right <- function(x, d) W - x / max(d$x) * tw2
e1 <- tree_edges(d1 |> mutate(x = to_left(x, d1)), ifelse(d1$discordant, col_discord, ink))
e2 <- tree_edges(d2 |> mutate(x = to_right(x, d2)), ifelse(d2$discordant, col_discord, ink))
tips <- bind_rows(
  d1 |> filter(isTip) |> mutate(xl = tw1 + off, hj = 0),
  d2 |> filter(isTip) |> mutate(xl = W - tw2 - off, hj = 1)) |>
  mutate(expr = vapply(tree_tip_label(label), plotmath_label, ""),
         col = ifelse(label == new_genome, col_new, ink))
links <- tibble(label = d1$label[d1$isTip], y1 = d1$y[d1$isTip]) |>
  mutate(y2 = d2$y[d2$isTip][match(label, d2$label[d2$isTip])],
         x1 = tw1 + off + lab_w + 0.04, x2 = W - tw2 - off - lab_w - 0.04)   # in the data: the .rds is redrawn elsewhere
# right tree with branch lengths: tips end short of the label column; dotted leaders close the gap
leaders <- d2 |> filter(isTip) |> mutate(x1 = to_right(x, d2) - 0.03, x2 = W - tw2 - off + 0.02) |>
  filter(x1 - x2 > 0.03)
# supports: left tree above the branch into the node, ending at the node; right tree at node height, on the
# tip side of the node (free space between the two child branches)
sups <- bind_rows(d1 |> filter(show) |>
                    mutate(xs = to_left(x, d1) - 0.03 - if (left_pies) r_mm / 25.4 else 0, hj = 1, vj = -0.35),
                  d2 |> filter(show) |> mutate(xs = to_right(x, d2) - 0.04, hj = 1, vj = 0.5))
n <- Ntip(t1)

fig <- ggplot() +
  geom_tree_edges(bind_rows(e1, e2)) +
  geom_segment(data = leaders, aes(x = x1, xend = x2, y = y, yend = y), colour = col_other, linewidth = 0.3,
               linetype = "dotted") +
  geom_segment(data = links, aes(x = x1, xend = x2, y = y1, yend = y2), colour = col_other, linewidth = 0.35) +
  geom_text(data = tips, aes(xl, y, label = expr, colour = col, hjust = hj), parse = TRUE, size = tip_pt / .pt) +
  geom_text(data = sups, aes(xs, y, label = sup, hjust = hj, vjust = vj), size = node_pt / .pt, colour = ink_muted) +
  scale_colour_identity() +
  scale_x_continuous(limits = c(-0.03, W + 0.03), expand = c(0, 0)) +   # root lines not clipped
  scale_y_continuous(expand = expansion(add = c(0.6, 0.5))) +
  theme_void(base_size = tip_pt) +
  theme(plot.background = element_rect(fill = "white", colour = NA), plot.margin = margin(4, 4, 4, 4))

## key: the quartet colours (with pies), then the discordant branch, a wider gap between the two groups. Short (no pies): centred between the two titles.
## With the pies it is too long for that and goes in a row at the bottom, left-aligned under the left tree
## (which carries the pies); a right phylogram's scale bar then sits under its title, which moves up for it
key <- tibble(sym = "line", col = col_discord, text = "Discordant branch")
if (left_pies) key <- bind_rows(tibble(sym = "dot", col = col_quartet,
  text = c("Main topology (q1)", "First alternative (q2)", "Second alternative (q3)")), key)
key_pt <- tip_pt - 0.5
key <- key |> mutate(sym_w = ifelse(sym == "line", 0.2, 0.1),
                     text_w = vapply(text, label_width_in, 0, size_pt = key_pt),
                     item_w = sym_w + 0.05 + text_w,
                     gap = ifelse(lead(sym) != sym, 0.3, 0.18))   # gap after each item; wider between groups
key_w <- sum(key$item_w) + sum(key$gap, na.rm = TRUE)
title_w <- vapply(c("ASTRAL-IV", "Concatenated ML"), label_width_in, 0, size_pt = tip_pt)
key_row <- W / 2 - key_w / 2 < title_w[1] + 0.25 || W / 2 + key_w / 2 > W - title_w[2] - 0.25
key <- key |> mutate(x0 = (if (key_row) 0 else W / 2 - key_w / 2) + cumsum(lag(item_w + gap, default = 0)))
bottom_y <- 0.1                                            # key row below the outgroup row
title_y <- if (key_row && right_phylo) n + 1.35 else n + 1  # room for the scale bar under the right title
fig <- fig + annotate("text", x = c(0, W), y = title_y, label = c("ASTRAL-IV", "Concatenated ML"),
                      hjust = c(0, 1), vjust = 0.5, size = tip_pt / .pt, fontface = "bold", colour = ink)
key <- key |> mutate(ky = if (key_row) bottom_y else n + 1)   # in the data: the .rds is redrawn elsewhere
fig <- fig +
  geom_segment(data = filter(key, sym == "line"), aes(x = x0, xend = x0 + sym_w, y = ky, yend = ky),
               colour = col_discord, linewidth = 0.6) +
  geom_point(data = filter(key, sym == "dot"), aes(x0 + sym_w / 2, ky, fill = col), shape = 21, size = 2.6,
             stroke = 0.3, colour = ink) +
  scale_fill_identity() +
  geom_text(data = key, aes(x0 + sym_w + 0.05, ky, label = text), hjust = 0, size = key_pt / .pt, colour = ink)

## pies on the left tree, a fixed physical radius so they stay round
if (left_pies) {
  pies <- d1 |> filter(!isTip) |>
    mutate(q1 = astral_value(label, "q1"), q2 = astral_value(label, "q2"), q3 = astral_value(label, "q3")) |>
    filter(!is.na(q1))
  for (i in seq_len(nrow(pies))) {
    xi <- to_left(pies$x[i], d1)
    fig <- fig + annotation_custom(pie_grob(unlist(pies[i, c("q1", "q2", "q3")]), r_mm = r_mm),
                                   xmin = xi - 0.5, xmax = xi + 0.5, ymin = pies$y[i] - 1, ymax = pies$y[i] + 1)
  }
}

if (right_phylo) {
  og <- d2 |> filter(label == outgroup)
  fig <- fig + break_mark((to_right(d2$x[d2$node == og$parent], d2) + to_right(og$x, d2)) / 2, og$y,
                          1, 1 / tree_tip_h)
  # scale bar: under the right title, right-aligned, label to its left (with the bottom key row); otherwise
  # between the outgroup row and the row above it
  sb <- 0.005; sb_in <- sb / max(d2$x) * tw2
  xr <- if (key_row) W else W - 0.25; sy <- if (key_row) n + 0.65 else 1.45
  fig <- fig +
    annotate("segment", x = xr - sb_in, xend = xr, y = sy, yend = sy, colour = ink, linewidth = 0.45) +
    if (key_row) annotate("text", x = xr - sb_in - 0.05, y = sy, label = paste(sb, "substitutions/site"),
                          hjust = 1, vjust = 0.5, size = key_pt / .pt, colour = ink)
    else annotate("text", x = xr, y = sy, label = paste(sb, "substitutions/site"), hjust = 1, vjust = -0.9,
                  size = node_pt / .pt, colour = ink)
}
w <- W + 0.06 + 8 / 72; h <- tree_height(n + 1) + (if (key_row) 0.9 * tree_tip_h else 0) + (title_y - n - 1) * tree_tip_h
save_fig(fig, prefix, w, h)
saveRDS(list(plot = fig, width = w, height = h), paste0(prefix, ".rds"))

## values behind the figure: clades in one tree only, and the Robinson-Foulds distance
rf <- dist.topo(unroot(t1), unroot(t2), method = "PH85")[1]   # unrooted Robinson-Foulds: splits differing
rf_norm <- rf / (2 * (n - 3))
bind_rows(d1 |> filter(discordant) |> transmute(tree = "ASTRAL-IV", clade, support = sup),
          d2 |> filter(discordant) |> transmute(tree = "Concatenated ML", clade, support = sup)) |>
  mutate(clade = vapply(strsplit(clade, ","), function(x) paste(tree_tip_label(x), collapse = ", "), "")) |>
  write_tsv(paste0(prefix, ".discordant.tsv"))
write_tsv(tibble(splits_differing = rf, splits_total = 2 * (n - 3), normalised_RF = round(rf_norm, 4)),
          paste0(prefix, ".rf.tsv"))
