# setup.R: customisation and shared functions for the funannotate compare plots.
# Sourced at the top of every plot script; change names, colours and genome order here, nowhere else.
# env: conda env r_env_DESeq2
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
  library(stringr); library(forcats); library(patchwork)
})

## ================================================================ customisation

# the new assembly: its ID in all three pipelines (params `name`: species_tree, introgression; comparison
# `genomes:`), and how it is named in figures and tables
new_genome       <- "Porites_sp_indonesia"
# new_genome_label <- "Porites sp. Indonesia"
new_genome_label <- "P. sp. Indonesia"

# genome order: the new assembly first, then the other Porites, then the Acropora outgroup.
# Count/composition plots show only `porites`; the outgroup is used where it adds information
# (orthology: which families are Porites-specific).
porites     <- c(new_genome, "Porites_lobata", "Porites_lutea", "Porites_evermanni", "Porites_harrisoni")
outgroup    <- "Acropora_millepora"
all_genomes <- c(porites, outgroup)

# text italicised wherever it appears in a label (regular expressions): genus names and abbreviated
# binomials ("P. lobata"), never "sp.": "Porites sp. Indonesia" -> *Porites* sp. Indonesia,
# "P. sp. Indonesia" -> *P.* sp. Indonesia.
italic_terms <- c("[A-Z]\\. (?!sp\\.)[a-z]+", "[A-Z]\\.(?= sp\\.)", "Porites", "Acropora")

# colours (dataviz reference palette): ink/grid neutrals, the accent of the new assembly (purple, as in
# the species tree plots), grey for the other genomes, a blue <-> red diverging pair with a grey
# midpoint, and a single-hue blue sequential ramp
col_new   <- "#6a3d9a"
col_other <- "#a3a29c"
ink <- "#2b2b29"; ink_muted <- "#6b6b66"; grid_col <- "#e6e5e1"; band <- "#f6f5f2"
col_under <- "#2a78d6"; col_over <- "#e34948"; col_mid <- "#f0efec"
seq_ramp  <- c("#f4f8fd", "#cfe0f6", "#94bdec", "#5596e0", "#2a78d6", "#1b56a0", "#123a6e")

# species tree plots (species_tree/plots): colours of the ASTRAL
# quartet support (main topology q1, alternatives q2 and q3) and of branches that differ between two trees
col_quartet <- c(q1 = "#2a78d6", q2 = "#eb6834", q3 = "#1baf7a")
col_discord <- col_over
# their text sizes (pt): tip labels, node values; height per tip (in)
tip_pt <- 8; node_pt <- 6.5; tree_tip_h <- 0.27

# genome assembly and annotation plots (genome-assembly-annotation/plots), from the categorical slots of the
# reference palette in fixed order: mitogenome gene groups; RNA genes, repeats and GC in the circular genome
# overview; BlobToolKit phyla (the host phylum first, then the next largest, then the rest pooled)
col_mito <- c(COX = "#2a78d6", ND = "#eb6834", CYTB = "#1baf7a", ATP = "#eda100", rRNA = "#e87ba4",
              tRNA = "#52514e")
col_genome_rna <- c(tRNA = "#52514e", ncRNA = "#eb6834", rRNA = "#e87ba4")
col_repeat <- c(classified = "#1baf7a", unclassified = col_other)
col_gc <- "#2a78d6"
col_link <- c(unplaced = col_other, between = "#eb6834", within = "#2a78d6")
col_phylum <- c("#2a78d6", "#eb6834", "#1baf7a")
col_kmer <- c(observed = "#94bdec", model = ink, unique = "#1baf7a", errors = "#eb6834")
# composite_2 (annotation): repeat classes in palette order, unknown families in the neutral grey; a single-hue
# ramp for repeat coverage (aqua, beside the blue seq_ramp of gene density); RNA-seq support of gene models
# (full -> none, dark -> light blue); OMArk categories
col_te_class <- c(DNA = "#2a78d6", LINE = "#eb6834", SINE = "#1baf7a", LTR = "#eda100", RC = "#e87ba4",
                  Unknown = col_other)
repeat_ramp <- c("#f2fbf7", "#c6ecdc", "#8fd8bb", "#4fc196", "#1baf7a", "#138a5f", "#0c5f41")
col_evidence <- c(full = "#1b56a0", partial = "#5596e0", none = "#cfe0f6")
col_omark <- c(Single = "#1b56a0", Duplicated = "#5596e0", Missing = "#cfe0f6",
               Consistent = "#1b56a0", Inconsistent = "#eda100", Contamination = "#e34948", Unknown = col_other)

# font of all figure text (PNG and PDF alike): a family fontconfig finds (`fc-list : family`), or one whose
# files are in fig_font_dir (registered for the session; needs Regular, Italic, Bold/SemiBold files).
# Arimo has the metrics of Helvetica/Arial. Example: fig_font <- "PP Mori"; fig_font_dir <- "~/fonts/PP_Mori"
fig_font     <- "PP Mori"
fig_font_dir <- "/home/colinl/Proj/ID-MetaGvClimate/PP_Mori"

## ================================================================ functions

# --- font. Registers fig_font from fig_font_dir (if given) for ragg/svglite, and points the generic "sans"
# at it, so text without an explicit family follows too
font_file <- function(dir, pattern) {
  f <- list.files(path.expand(dir), pattern = pattern, full.names = TRUE, ignore.case = TRUE)
  if (length(f) == 0) stop("no font file matching ", pattern, " in ", dir)
  f[1]
}
fig_font_files <- if (is.null(fig_font_dir)) {
  lapply(list(plain = list(), italic = list(italic = TRUE), bold = list(weight = "bold"),
              bolditalic = list(italic = TRUE, weight = "bold")),
         function(a) do.call(systemfonts::match_fonts, c(list(fig_font), a))$path)
} else {
  list(plain = font_file(fig_font_dir, "-Regular\\.(otf|ttf)$"),
       italic = font_file(fig_font_dir, "-(Regular)?Italic\\.(otf|ttf)$"),
       bold = font_file(fig_font_dir, "-(Semi)?Bold\\.(otf|ttf)$"),
       bolditalic = font_file(fig_font_dir, "-(Semi)?BoldItalic\\.(otf|ttf)$"))
}
if (!fig_font %in% systemfonts::system_fonts()$family &&
    !fig_font %in% systemfonts::registry_fonts()$family) do.call(systemfonts::register_font, c(list(fig_font), fig_font_files))
if (!"sans" %in% systemfonts::registry_fonts()$family) do.call(systemfonts::register_font, c(list("sans"), fig_font_files))

# genome ID -> plain-text label: "Porites_lobata" -> "P. lobata", the new genome -> new_genome_label
pretty_genome <- function(x) {
  w <- str_split_fixed(x, "_", 2)
  ifelse(x == new_genome, new_genome_label, paste0(substr(w[, 1], 1, 1), ". ", w[, 2]))
}
genome_lab <- setNames(pretty_genome(all_genomes), all_genomes)

# --- italics. Labels are plain text; these turn them into plotmath (figures) or Markdown (tables),
# with the italic_terms in italics and everything else upright.
italic_regex <- paste0("(", paste(italic_terms, collapse = "|"), ")")
pm_quote <- function(x) paste0('"', gsub('(["\\\\])', "\\\\\\1", x), '"')
plotmath_line <- function(s) {
  m <- str_locate_all(s, italic_regex)[[1]]
  if (nrow(m) == 0) return(pm_quote(s))
  parts <- character(); pos <- 1
  for (i in seq_len(nrow(m))) {
    if (m[i, 1] > pos) parts <- c(parts, pm_quote(substr(s, pos, m[i, 1] - 1)))
    parts <- c(parts, paste0("italic(", pm_quote(substr(s, m[i, 1], m[i, 2])), ")"))
    pos <- m[i, 2] + 1
  }
  if (pos <= nchar(s)) parts <- c(parts, pm_quote(substr(s, pos, nchar(s))))
  paste(parts, collapse = "*")
}
# one label -> plotmath string; a "\n" gives two lines (atop)
plotmath_label <- function(s) {
  l <- str_split(s, "\n")[[1]]
  if (length(l) == 1) plotmath_line(l) else paste0("atop(", plotmath_line(l[1]), ",",
                                                    plotmath_line(paste(l[-1], collapse = " ")), ")")
}
# labels -> expression vector, for scale `labels =` and legend titles
label_expr <- function(s) parse(text = vapply(as.character(s), plotmath_label, ""))
# facet labeller with the same italics
labeller_italic <- function(labels) {
  labels[] <- lapply(labels, function(v) vapply(as.character(v), plotmath_label, ""))
  label_parsed(labels)
}
# Markdown: italic terms wrapped in *...*
md_label <- function(s) str_replace_all(s, italic_regex, "*\\1*")

# genome axes: italics from label_expr, the new assembly first (top on a y axis)
scale_x_genome <- function(...) scale_x_discrete(labels = function(x) label_expr(genome_lab[x]), ...)
scale_y_genome <- function(genomes = NULL, ...) {
  if (is.null(genomes)) scale_y_discrete(labels = function(x) label_expr(genome_lab[x]), ...)
  else scale_y_discrete(limits = rev(genomes), labels = function(x) label_expr(genome_lab[x]), ...)
}

theme_fig <- function(base_size = 9) {
  theme_minimal(base_size = base_size, base_family = fig_font) +
    theme(text = element_text(colour = ink), axis.text = element_text(colour = ink_muted),
          panel.grid.major = element_line(colour = grid_col, linewidth = 0.3),
          panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold", hjust = 0, colour = ink),
          plot.background = element_rect(fill = "white", colour = NA),
          legend.title = element_text(size = 8), legend.text = element_text(size = 7))
}

# every figure is written as PDF and as 300 dpi PNG, both laid out with the same font metrics (systemfonts):
# PNG with ragg; PDF as SVG (svglite) converted by rsvg-convert, fonts embedded. Not cairo_pdf: through Pango
# it rounds glyph positions to whole points, which spaces letters unevenly at small sizes.
save_fig <- function(fig, prefix, width, height) {
  ggsave(paste0(prefix, ".png"), fig, width = width, height = height, dpi = 300, bg = "white",
         device = ragg::agg_png, limitsize = FALSE)
  svg <- tempfile(fileext = ".svg"); on.exit(unlink(svg))
  ggsave(svg, fig, width = width, height = height, device = svglite::svglite, bg = "white",
         system_fonts = list(sans = fig_font), limitsize = FALSE)
  s <- readLines(svg)
  s <- sub("<svg ", "<svg xml:space='preserve' ", s, fixed = TRUE)   # keep spaces between plotmath pieces
  writeLines(s, svg)
  # rsvg finds fonts through fontconfig: add fig_font_dir to the configuration for this call
  env <- character()
  if (!is.null(fig_font_dir)) {
    base <- c(Sys.getenv("FONTCONFIG_FILE"), file.path(Sys.getenv("CONDA_PREFIX"), "etc", "fonts", "fonts.conf"),
              "/etc/fonts/fonts.conf")
    base <- base[nzchar(base) & file.exists(base)][1]
    conf <- tempfile(fileext = ".conf"); on.exit(unlink(conf), add = TRUE)
    writeLines(c('<?xml version="1.0"?>', '<!DOCTYPE fontconfig SYSTEM "fonts.dtd">', "<fontconfig>",
                 sprintf('  <include ignore_missing="yes">%s</include>', base),
                 sprintf("  <dir>%s</dir>", normalizePath(fig_font_dir)), "</fontconfig>"), conf)
    env <- paste0("FONTCONFIG_FILE=", conf)
  }
  rsvg <- Sys.which("rsvg-convert")
  if (!nzchar(rsvg)) stop("rsvg-convert not found (conda env r_env_DESeq2 has it)")
  st <- system2(rsvg, c("-f", "pdf", "-o", shQuote(paste0(prefix, ".pdf")), shQuote(svg)), env = env,
                stderr = FALSE)
  if (st != 0) stop("rsvg-convert failed for ", prefix)
}

# --- species trees (ape). Tip IDs "Plobata1_GCA_973946955.1" (genus letter, species, replicate number,
# accession); the new assembly is new_genome.
# tip ID -> plain-text label: "P. lobata (GCA_973946955.1)", the new assembly -> new_genome_label
tree_tip_label <- function(x) {
  m <- str_match(x, "^([A-Z])([a-z]+)[0-9]*_(GC[AF]_[0-9.]+)$")
  ifelse(x == new_genome, new_genome_label,
         ifelse(is.na(m[, 1]), x, sprintf("%s. %s (%s)", m[, 2], m[, 3], m[, 4])))
}
# clade of every internal node (node Ntip + i), as sorted tip IDs; with `outgroup`, a clade containing it
# is replaced by its complement, so that keys match between differently rooted trees
tree_clades <- function(t, outgroup = NULL) {
  vapply(ape::prop.part(t), function(i) {
    tips <- t$tip.label[i]
    if (!is.null(outgroup) && outgroup %in% tips) tips <- setdiff(t$tip.label, tips)
    paste(sort(tips), collapse = ",")
  }, "")
}
# read a tree, root it on the outgroup and ladderize it; node labels (supports) are read off the tree as
# written and put back by clade, because rerooting shifts them
read_tree_rooted <- function(file, outgroup) {
  t <- ape::read.tree(file)
  sup <- setNames(c(t$node.label, rep("", t$Nnode))[seq_len(t$Nnode)], tree_clades(t, outgroup))
  t <- ape::ladderize(ape::root(t, outgroup, resolve.root = TRUE), right = FALSE)
  t$node.label <- unname(sup[tree_clades(t, outgroup)])
  t$node.label[is.na(t$node.label)] <- ""
  t
}
# printed width (in) of the widest of the plain-text labels, italics included, at `size_pt`; used to size
# figures so that tip labels get exactly the space they need
label_width_in <- function(s, size_pt) {
  f <- tempfile(fileext = ".png")
  ragg::agg_png(f, width = 1, height = 1, units = "in", res = 72)   # measured with the figure font's metrics
  on.exit({ grDevices::dev.off(); unlink(f) })
  max(vapply(label_expr(s), function(e) grid::convertWidth(grid::grobWidth(
    grid::textGrob(e, gp = grid::gpar(fontsize = size_pt, fontfamily = fig_font))), "in", valueOnly = TRUE), 0))
}
# a value out of an ASTRAL node label "[...;localPP=1;q1=0.43;...]"
astral_value <- function(label, key) {
  as.numeric(str_match(label, paste0("\\b", key, "=([0-9.eE+-]+)"))[, 2])
}

# node coordinates of a tree (ggtree layout: tips at y = 1..n, x = depth), nodes kept in the order of `tr`
tree_layout <- function(tr, ...) {
  as_tibble(as.data.frame(ggtree::ggtree(tr, ladderize = FALSE, ...)$data))
}
# rectangular edges of a layout as segments: horizontal into each node, vertical at its parent;
# `colour` = one colour per node (the edge leading to it)
tree_edges <- function(d, colour = rep(ink, nrow(d))) {
  e <- d |> mutate(colour = colour) |> filter(node != parent)
  p <- match(e$parent, d$node)
  h <- tibble(x = d$x[p], xend = e$x, y = e$y, yend = e$y, colour = e$colour, node = e$node, part = "h")
  # one continuous vertical per parent (two halves meeting at the node leave a seam), in ink; the half of a
  # coloured branch is drawn over it
  v <- e |> mutate(px = d$x[p], py = d$y[p]) |> group_by(parent) |>
    summarise(x = first(px), y0 = min(y, first(py)), y1 = max(y, first(py)), .groups = "drop") |>
    transmute(x, xend = x, y = y0, yend = y1, colour = ink, node = parent, part = "v")
  vc <- tibble(x = d$x[p], xend = d$x[p], y = d$y[p], yend = e$y, colour = e$colour, node = e$node,
               part = "vc") |> filter(colour != ink)
  bind_rows(v, vc, h)
}
# layers drawing tree_edges(): verticals with butt ends, horizontals with square ends (closed corners)
geom_tree_edges <- function(e, linewidth = 0.45) {
  list(geom_segment(data = filter(e, part != "h"), aes(x, y, xend = xend, yend = yend, colour = colour),
                    linewidth = linewidth, lineend = "butt"),
       geom_segment(data = filter(e, part == "h"), aes(x, y, xend = xend, yend = yend, colour = colour),
                    linewidth = linewidth, lineend = "square"))
}
tree_height <- function(n_tip, extra = 0.3) n_tip * tree_tip_h + extra
# tree with tip labels: edges, labels in italics where due, the new assembly in its accent colour.
# The tree takes `tree_w` in, the labels exactly their width; the figure width is attr(, "fig_w").
tree_plot <- function(d, tree_w, edge_colour = rep(ink, nrow(d)), linewidth = 0.45) {
  tips <- d |> filter(isTip) |>
    mutate(lab = tree_tip_label(label), expr = vapply(lab, plotmath_label, ""),
           col = ifelse(label == new_genome, col_new, ink))
  lw <- label_width_in(tips$lab, tip_pt); off <- 0.06
  x0 <- min(d$x); upi <- (max(d$x) - x0) / tree_w            # data units per inch
  p <- ggplot() +
    geom_tree_edges(tree_edges(d, edge_colour), linewidth) +
    geom_text(data = tips, aes(x + off * upi, y, label = expr, colour = col), parse = TRUE, hjust = 0,
              size = tip_pt / .pt) +
    scale_colour_identity() +
    scale_x_continuous(limits = c(x0 - 0.03 * upi, max(d$x) + (off + lw + 0.02) * upi), expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(add = 0.6)) +
    theme_void(base_size = tip_pt, base_family = fig_font) +
    theme(plot.background = element_rect(fill = "white", colour = NA),
          plot.margin = margin(4, 4, 4, 4), text = element_text(colour = ink))
  attr(p, "fig_w") <- 0.03 + tree_w + off + lw + 0.02 + 8 / 72   # 0.03: the root line is not clipped
  attr(p, "upi") <- upi
  p
}
# rooted tree with branch lengths: rooting puts the root on the outgroup branch with a zero-length ingroup
# stem, and the outgroup is far away; the whole path between the ingroup crown and the outgroup is redrawn
# shortened (a short stem plus the outgroup branch, to be marked with break_mark()). Real and drawn lengths
# are kept as attributes.
shorten_outgroup <- function(tr, outgroup, stem = 0.08, og = 0.3) {
  og_tip <- which(tr$tip.label == outgroup)
  og_edge <- which(tr$edge[, 2] == og_tip)
  stem_edge <- which(tr$edge[, 1] == ape::Ntip(tr) + 1 & tr$edge[, 2] != og_tip)
  depth <- max(ape::node.depth.edgelength(ape::drop.tip(tr, outgroup)))
  real <- sum(tr$edge.length[c(og_edge, stem_edge)])
  tr$edge.length[stem_edge] <- stem * depth
  tr$edge.length[og_edge] <- og * depth
  attr(tr, "crown_to_outgroup") <- real; attr(tr, "drawn") <- (stem + og) * depth
  attr(tr, "ingroup_depth") <- depth
  tr
}
# break mark (gap and two slanted strokes) on a horizontal branch at (x, y); upi/yi: data units per inch
break_mark <- function(x, y, upi, yi) {
  gap <- 0.03 * upi; dx <- 0.025 * upi; dy <- 0.06 * yi
  list(annotate("rect", xmin = x - gap / 2, xmax = x + gap / 2, ymin = y - 0.4 * dy, ymax = y + 0.4 * dy,
                fill = "white", colour = NA),
       annotate("segment", x = x + c(-gap, gap) / 2 - dx / 2, xend = x + c(-gap, gap) / 2 + dx / 2,
                y = y - dy / 2, yend = y + dy / 2, colour = ink, linewidth = 0.45))
}
# pie of shares `q`, physical radius `r_mm`, centred in its viewport (annotation_custom)
pie_grob <- function(q, r_mm, cols = col_quartet) {
  a <- 2 * pi * c(0, cumsum(q / sum(q)))
  wedges <- lapply(seq_along(q), function(i) {
    t <- pi / 2 - seq(a[i], a[i + 1], length.out = max(2, ceiling(60 * q[i])))   # clockwise from 12 o'clock
    grid::polygonGrob(grid::unit(0.5, "npc") + grid::unit(c(0, r_mm * cos(t)), "mm"),
                      grid::unit(0.5, "npc") + grid::unit(c(0, r_mm * sin(t)), "mm"),
                      gp = grid::gpar(fill = cols[i], col = "white", lwd = 0.4))
  })
  ring <- grid::circleGrob(r = grid::unit(r_mm, "mm"), gp = grid::gpar(fill = NA, col = ink, lwd = 0.5))
  do.call(grid::grobTree, c(wedges, list(ring)))
}

# funannotate compare tables: first column = row names, one column per genome
read_fc_table <- function(path) {
  d <- read_csv(path, show_col_types = FALSE, name_repair = "minimal")
  names(d)[1] <- "id"
  d
}
