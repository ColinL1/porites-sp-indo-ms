# gene_support.R: standalone figure on the release gene models, as 100% bars: protein-coding genes by extrinsic
# evidence of their model (BRAKER gene_support: full, partial or no hint support; release GFF tag evidence=);
# OMArk completeness of the conserved HOGs (single, duplicated, missing); OMArk consistency of the proteome
# (consistent, inconsistent, contamination, unknown).
# usage: Rscript gene_support.R release.gff3 omark.sum out_prefix
# env: conda env r_env_DESeq2
args <- commandArgs(trailingOnly = TRUE)
here <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))))
source(file.path(here, "..", "..", "..", "shared", "plots", "setup.R"))
gff <- args[1]; omark_sum <- args[2]; out <- args[3]

g <- read_tsv(gff, comment = "#", col_names = FALSE, col_types = cols(.default = "c"), progress = FALSE) |>
  filter(X3 == "gene", str_detect(X9, "evidence=")) |>
  transmute(evidence = str_match(X9, "evidence=([a-z]+)")[, 2])
ev <- g |> count(evidence) |> mutate(pct = 100 * n / sum(n))
s <- readLines(omark_sum)
num <- function(line, key) as.numeric(str_match(line, paste0("\\b", key, ":([0-9.]+)%"))[, 2])
comp_line <- grep("^S:[0-9.]+%", s, value = TRUE)
cons_line <- grep("^A:[0-9.]+%", s, value = TRUE)
bars <- bind_rows(
  tibble(bar = "Gene models,\nextrinsic evidence", cat = c("full", "partial", "none"),
         lab = c("Full", "Partial", "None"), pct = ev$pct[match(c("full", "partial", "none"), ev$evidence)]),
  tibble(bar = "OMArk,\nconserved HOGs", cat = c("Single", "Duplicated", "Missing"),
         lab = c("Single", "Duplicated", "Missing"), pct = c(num(comp_line, "S"), num(comp_line, "D"), num(comp_line, "M"))),
  tibble(bar = "OMArk,\nproteome", cat = c("Consistent", "Inconsistent", "Contamination", "Unknown"),
         lab = c("Consistent", "Inconsistent", "Contamination", "Unknown"),
         pct = c(num(cons_line, "A"), num(cons_line, "I"), num(cons_line, "C"), num(cons_line, "U")))
) |> mutate(bar = factor(bar, rev(unique(bar)))) |>
  group_by(bar) |> mutate(x1 = cumsum(pct), x0 = x1 - pct) |> ungroup()
cols <- c(col_evidence, col_omark)
bars <- bars |> mutate(dark = cat %in% c("full", "Single", "Consistent", "Duplicated", "partial", "Contamination"),
                       inside = pct >= 20)
# labels: inside wide segments; narrow ones (other than 0%) named above the bar, at the segment
# (a second narrow segment in the same bar is named below it)
small <- bars |> filter(!inside, pct > 0) |> group_by(bar) |>
  mutate(above = row_number() == 1, y = as.integer(bar) + ifelse(above, 0.38, -0.38)) |> ungroup()
p_sup <- ggplot(bars) +
  geom_rect(aes(xmin = x0, xmax = x1, ymin = as.integer(bar) - 0.34, ymax = as.integer(bar) + 0.34, fill = cat),
            colour = "white", linewidth = 0.4) +
  geom_text(data = filter(bars, inside), aes((x0 + x1) / 2, as.integer(bar), label = sprintf("%s\n%.0f%%", lab, pct),
                                              colour = dark), size = 7 / .pt, lineheight = 0.9, family = fig_font) +
  geom_text(data = small, aes((x0 + x1) / 2, y, label = sprintf("%s %.1f%%", lab, pct), vjust = ifelse(above, 0, 1)),
            hjust = 0.8, size = 6.5 / .pt, colour = ink_muted, family = fig_font) +
  scale_fill_manual(values = cols, guide = "none") +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = ink), guide = "none") +
  scale_x_continuous("Share (%)", breaks = seq(0, 100, 25), expand = expansion(mult = c(0, 0.03))) +
  scale_y_continuous(NULL, breaks = seq_along(levels(bars$bar)), labels = levels(bars$bar),
                     expand = expansion(add = c(0.6, 0.6))) +
  coord_cartesian(clip = "off") +
  theme_fig() + theme(panel.grid.major.y = element_blank(), axis.text.y = element_text(colour = ink, hjust = 1))
save_fig(p_sup, out, 4.8, 2.0)
write_tsv(bars |> select(bar, cat, pct), paste0(out, ".tsv"))

