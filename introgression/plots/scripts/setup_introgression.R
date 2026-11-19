# setup_introgression.R: what the introgression plots need on top of the shared setup.R.
# Sources shared/plots/setup.R (names, colours, italics, theme, save_fig) and adds only
# the introgression pipeline's sample IDs and the roles of the shared colours. Nothing in setup.R is
# redefined here; change names and colours there.
# env: conda env r_env_DESeq2
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(script_dir, "..", "..", "..", "shared", "plots", "setup.R"))

## ================================================================ samples of the introgression pipeline

# the target's sample ID in the pipeline (params name): new_genome of setup.R
target <- new_genome
p1 <- "Plobata1"; p2 <- "Pharrisoni"   # parent1 / parent2 in references.tsv

# sample ID -> label. The species name alone where a species has one genome; "1"/"2" where it has two.
sample_lab <- c(setNames(new_genome_label, target),
                Plobata1 = "P. lobata 1", Plobata2 = "P. lobata 2", Pharrisoni = "P. harrisoni",
                Pcylindrica1 = "P. cylindrica", Paustraliensis = "P. australiensis",
                Pdivaricata2 = "P. divaricata", Plutea2 = "P. lutea")
# where only one P. lobata genome is involved (ancestry painting, D, Twisst: parent1)
sample_lab_short <- replace(sample_lab, "Plobata1", "P. lobata")
# sample order in figures: target, parents, the other P. lobata, outgroups, P. lutea
sample_order <- c(target, p1, p2, "Plobata2", "Pcylindrica1", "Paustraliensis", "Pdivaricata2", "Plutea2")
# outgroups of D-statistics and Twisst, in references.tsv order
outgroups <- c("Pcylindrica1", "Pdivaricata2", "Paustraliensis")

## ================================================================ colour roles (values from setup.R)

# the two parents are the two poles of setup.R's diverging pair; the neutral between them (shared
# by both: heterozygous opposite, species-tree topology) is the grey of the other genomes
col_p1   <- col_under    # P. lobata
col_p2   <- col_over     # P. harrisoni
col_both <- col_other
# ADMIXTURE components: dataviz categorical slots 2-6 (orange, aqua, yellow, magenta, green): not the blue
# and red of the parents (an ADMIXTURE component is not a parent), and not violet
# (slot 7) is left out so purple stays the accent of the new genome
col_k <- c("#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300")

## ================================================================ shared helpers

# scaffold set (00_inputs/scaffolds.bed) -> offsets on one concatenated genome axis, with a gap
# between scaffolds; `mid` is where the scaffold number is printed
genome_axis <- function(bed, gap = 3e6) {
  read_tsv(bed, col_names = c("chrom", "start", "end"), show_col_types = FALSE) |>
    mutate(len = end - start, offset = cumsum(lag(len + gap, default = 0)),
           mid = offset + len / 2, label = str_remove(chrom, "^scaffold_"))
}
scale_x_genome_axis <- function(ax, ...)
  scale_x_continuous(breaks = ax$mid, labels = ax$label, expand = expansion(mult = 0.005), ...)

# sample labels (plain) -> plotmath, for axes and geom_text(parse = TRUE)
sample_expr <- function(ids, labs = sample_lab) label_expr(labs[ids])

# key (legend) drawn as a grob, for keys whose entries need line breaks: plotmath can only centre two lines
# (atop), so each line is its own left-aligned plotmath text (italics as everywhere else), the symbol centred
# on the entry. `sections`: list of list(title, labels (plain text, "\n" = line break), colours, shapes (pch));
# sizes in pt.
key_grob <- function(sections, text_pt = 6.5, title_pt = 7, key_mm = 3, gap_mm = 2.5) {
  pt_mm <- 0.3528
  line_mm <- text_pt * pt_mm * 1.3
  g <- gtable::gtable(widths = grid::unit(c(key_mm + 1, 1), c("mm", "null")))
  txt <- function(s, pt) grid::textGrob(label_expr(s)[[1]], x = 0, hjust = 0, vjust = 0.5,
                                        gp = grid::gpar(fontfamily = fig_font, fontsize = pt, col = ink))
  add_row <- function(g, h_mm) gtable::gtable_add_rows(g, grid::unit(h_mm, "mm"))
  for (k in seq_along(sections)) {
    s <- sections[[k]]
    if (k > 1) g <- add_row(g, gap_mm)
    g <- add_row(g, title_pt * pt_mm * 1.6)
    g <- gtable::gtable_add_grob(g, txt(s$title, title_pt), t = nrow(g), l = 1, r = 2)
    for (i in seq_along(s$labels)) {
      lines <- strsplit(s$labels[i], "\n")[[1]]
      if (i > 1) g <- add_row(g, 1)                       # space between entries
      top <- nrow(g) + 1
      for (l in lines) {
        g <- add_row(g, line_mm)
        g <- gtable::gtable_add_grob(g, txt(l, text_pt), t = nrow(g), l = 2)
      }
      g <- gtable::gtable_add_grob(g, grid::pointsGrob(0.5, 0.5, pch = s$shapes[i], size = grid::unit(1.6, "mm"),
                                                       gp = grid::gpar(col = s$colours[i], fill = s$colours[i])),
                                   t = top, b = nrow(g), l = 1)
    }
  }
  g
}
