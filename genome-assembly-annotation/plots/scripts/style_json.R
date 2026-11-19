# style_json.R: write the customisation of setup.R (names, colours, font, italic terms) to a JSON file for the
# Python figure scripts (fig_style.py reads it), so that nothing is defined twice.
# usage: Rscript style_json.R out.json
# env: conda env r_env_DESeq2
args <- commandArgs(trailingOnly = TRUE)
here <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))))
source(file.path(here, "..", "..", "..", "shared", "plots", "setup.R"))

style <- list(
  new_genome_label = new_genome_label, italic_terms = italic_terms,
  font = list(family = fig_font, files = lapply(fig_font_files, normalizePath)),
  col = list(new = col_new, other = col_other, ink = ink, ink_muted = ink_muted, grid = grid_col, band = band,
             seq_ramp = seq_ramp, mito = as.list(col_mito), genome_rna = as.list(col_genome_rna),
             repeat_ = as.list(col_repeat), gc = col_gc, link = as.list(col_link))
)
names(style$col)[names(style$col) == "repeat_"] <- "repeat"
jsonlite::write_json(style, args[1], auto_unbox = TRUE, pretty = TRUE)
