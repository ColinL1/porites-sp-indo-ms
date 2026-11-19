#!/usr/bin/env Rscript
# Assembly and annotation statistics from funannotate compare (stats/genome.stats.summary.csv),
# Porites genomes only, the new assembly first. A table, not a figure: few numbers in different units.
#   <prefix>_table.tsv   raw values (metric x genome), read by composite.R
#   <prefix>_table.md    display version for the manuscript: formatted numbers, *italic* species names, group
#                        rows "Assembly" / "Gene annotation" (| **Section** | | ...) with their metrics
#                        indented below them ("&nbsp;&nbsp;" prefix)
#   <prefix>_table_transposed.md   the layout of the composite table panel: genomes as rows, five metrics
#                        as columns (an alternative to choose from; both go into the Word file)
#   titles.tsv, notes.tsv (next to the .md)   title and footnotes, read by tables_docx.R
# "Proteins without an ortholog" is counted here from the ProteinOrtho table (protortho/funannotate.poff.tsv):
# proteins in rows holding a single genome. funannotate's "Unique Proteins" counts those rows (groups), not
# proteins, so it is not used. With + without = the protein total (incl. isoforms).
# tRNA counts are left out: three of the input annotations carry no tRNA genes.
# usage: genome_stats.R <funannotate_compare dir> <out prefix>
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(script_dir, "..", "..", "..", "shared", "plots", "setup.R"))   # names, colours, genome order, shared functions

args <- commandArgs(trailingOnly = TRUE)
fc_dir <- args[1]; prefix <- args[2]

raw <- read_fc_table(file.path(fc_dir, "stats", "genome.stats.summary.csv"))
num <- function(x) as.numeric(str_remove_all(x, "[^0-9.]"))

metrics <- tribble(
  ~id,                          ~label,                          ~scale, ~digits, ~section,
  "Assembly Size",              "Assembly size (Mb)",            1e-6,   0,       "Assembly",
  "Scaffold N50",               "Scaffold N50 (Mb)",             1e-6,   1,       "Assembly",
  "Largest Scaffold",           "Largest scaffold (Mb)",         1e-6,   1,       "Assembly",
  "Num Scaffolds",              "Scaffolds",                     1,      0,       "Assembly",
  "Percent GC",                 "GC (%)",                        1,      1,       "Assembly",
  "pc_genes",                   "Protein-coding genes",          1,      0,       "Gene annotation",
  "Num Proteins",               "Proteins (incl. isoforms)",     1,      0,       "Gene annotation",
  "Prots atleast 1 ortholog",   "Proteins with an ortholog",     1,      0,       "Gene annotation",
  "no_ortholog",                "Proteins without an ortholog",  1,      0,       "Gene annotation"
)

# proteins without an ortholog: proteins in ProteinOrtho rows that contain only one genome
poff <- read_tsv(file.path(fc_dir, "protortho", "funannotate.poff.tsv"), show_col_types = FALSE,
                 col_types = cols(.default = "c")) |> rename_with(~ str_remove(.x, "\\.faa$"))
n_prot <- function(x) ifelse(x == "*", 0, str_count(x, ",") + 1)
no_orth <- poff |> filter(`# Species` == "1") |> select(-(1:3)) |>
  summarise(across(everything(), ~ sum(n_prot(.x)))) |> mutate(id = "no_ortholog", .before = 1) |>
  mutate(across(-id, as.character))
# protein-coding genes: distinct gene IDs in annotations/<genome>.all.annotations.tsv (funannotate's
# "Num Genes" also counts tRNA/rRNA genes where annotated)
pc <- vapply(names(no_orth)[-1], function(g)
  as.character(n_distinct(read_tsv(file.path(fc_dir, "annotations", paste0(g, ".all.annotations.tsv")),
                                   col_types = cols(.default = "c"), col_select = 1)[[1]])), "")
raw <- bind_rows(raw, no_orth, as_tibble(as.list(c(id = "pc_genes", pc))))

stats <- raw |> filter(id %in% metrics$id) |>
  pivot_longer(-id, names_to = "genome", values_to = "raw") |>
  filter(genome %in% porites) |>
  left_join(metrics, by = "id") |>
  mutate(value = num(raw) * scale,
         genome = factor(genome, levels = porites), label = factor(label, levels = metrics$label)) |>
  arrange(label, genome)

## raw values (metric x genome)
stats |> select(Metric = label, genome, value) |>
  mutate(value = round(value, 6)) |>
  pivot_wider(names_from = genome, values_from = value) |>
  write_tsv(paste0(prefix, "_table.tsv"))

## display table (Markdown): formatted numbers, italic species names, a group row per section
fmt <- function(v, d) formatC(v, format = "f", digits = d, big.mark = ",")
disp <- stats |> mutate(txt = mapply(fmt, value, digits)) |>
  select(section, Metric = label, genome, txt) |>
  pivot_wider(names_from = genome, values_from = txt)
row_md <- function(cells) paste("|", paste(cells, collapse = " | "), "|")
md <- c(row_md(c("Metric", md_label(genome_lab[porites]))),
        row_md(c("---", rep("---:", length(porites)))))
for (sec in unique(metrics$section)) {
  md <- c(md, row_md(c(paste0("**", sec, "**"), rep("", length(porites)))))
  d <- disp |> filter(section == sec) |> mutate(Metric = paste0("&nbsp;&nbsp;", Metric))   # indented under it
  md <- c(md, apply(d |> select(-section), 1, row_md))
}
writeLines(md, paste0(prefix, "_table.md"))

## transposed version, as in the composite table panel (composite.R): genomes as rows, five metrics as columns
tcols <- c("Assembly size (Mb)" = "Assembly (Mb)", "Scaffold N50 (Mb)" = "N50 (Mb)", "Protein-coding genes" = "Genes",
           "Proteins with an ortholog" = "Proteins w/ ortholog",
           "Proteins without an ortholog" = "Proteins w/o ortholog")
tr <- stats |> filter(label %in% names(tcols)) |> mutate(txt = mapply(fmt, value, digits)) |>
  mutate(col = factor(tcols[as.character(label)], levels = tcols)) |>
  select(genome, col, txt) |> arrange(col) |>
  pivot_wider(names_from = col, values_from = txt) |> arrange(genome)
writeLines(c(row_md(c("Genome", tcols)), row_md(c("---", rep("---:", length(tcols)))),
             apply(tr |> mutate(genome = md_label(genome_lab[as.character(genome)])), 1, row_md)),
           paste0(prefix, "_table_transposed.md"))

## titles and footnotes for tables_docx.R (names from setup.R), one entry per table
name   <- basename(paste0(prefix, "_table"))
name_t <- paste0(name, "_transposed")
title  <- md_label(paste("Assembly and annotation statistics of", new_genome_label, "and four other Porites genomes"))
note_src <- md_label(paste("Statistics from funannotate compare. Gene models are those of each source",
                           "annotation (see Methods)."))
note_orth <- md_label(paste("Proteins with / without an ortholog: ProteinOrtho (funannotate compare) across",
                            "the five Porites genomes and A. millepora; together they make up the protein total",
                            "(incl. isoforms)."))
write_tsv(tibble(table = c(name, name_t), title = title), file.path(dirname(prefix), "titles.tsv"))
write_tsv(tibble(table = c(name, name, name_t, name_t),
                 note = c(note_src, note_orth, paste(note_src, "N50: scaffold N50."), note_orth)),
          file.path(dirname(prefix), "notes.tsv"))
