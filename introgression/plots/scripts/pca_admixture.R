#!/usr/bin/env Rscript
# PCA (PC1 against PC2) and ADMIXTURE (K = 2 to 5) of the 8 genomes, from 09_pca_admixture.
# ADMIXTURE components are coloured in the order of the first genome (sample_order) where each is
# largest, so a component keeps its colour across K as long as the same genome carries it.
# usage: pca_admixture.R <results dir> <out prefix>
# writes <prefix>.pdf/.png, <prefix>.admixture.tsv (plotted Q values), <prefix>.panels.rds
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(script_dir, "setup_introgression.R"))
suppressPackageStartupMessages(library(ggrepel))

args <- commandArgs(trailingOnly = TRUE)
res <- args[1]; prefix <- args[2]
dir <- file.path(res, "09_pca_admixture")

## ---- PCA
pca <- read_tsv(file.path(dir, "pca.tsv"), show_col_types = FALSE)
ve <- unlist(filter(pca, sample == "variance_explained")[, c("PC1", "PC2")])
pca <- filter(pca, sample != "variance_explained") |>
  mutate(is_new = sample == target,
         crowded = !is_new & abs(PC1 - median(PC1)) < 0.1 & abs(PC2 - median(PC2)) < 0.1,
         lab = vapply(sample_lab[sample], plotmath_label, ""))
pPCA <- ggplot(pca, aes(PC1, PC2)) +
  geom_point(data = filter(pca, !is_new), colour = col_other, size = 1.6) +
  geom_point(data = filter(pca, is_new), colour = col_new, size = 2.2) +
  # genomes stacked on one spot: labels in a column to their lower left (the y axis is extended to make room)
  geom_text_repel(data = filter(pca, crowded), aes(label = lab), parse = TRUE, size = 2.3, seed = 1, family = fig_font,
                  colour = ink, nudge_x = -0.3, direction = "y", hjust = 1, ylim = c(NA, -0.05),
                  min.segment.length = 0, segment.colour = ink_muted, segment.size = 0.2,
                  box.padding = 0.1, max.overlaps = Inf) +
  # the other genomes: new genome above its point, P. divaricata (top) left of it, P. lutea above it
  geom_text(data = filter(pca, is_new), aes(x = PC1 - 0.05, y = PC2 + 0.09, label = lab),
            parse = TRUE, size = 2.3, hjust = 0, vjust = 0, colour = col_new) +
  geom_text(data = filter(pca, !crowded, !is_new, PC2 == max(PC2)), aes(x = PC1 - 0.07, label = lab),
            parse = TRUE, size = 2.3, hjust = 1, vjust = 0.5, colour = ink) +
  geom_text(data = filter(pca, !crowded, !is_new, PC2 != max(PC2)), aes(y = PC2 + 0.07, label = lab),
            parse = TRUE, size = 2.3, hjust = 0.5, vjust = 0, colour = ink) +
  scale_x_continuous(expand = expansion(mult = 0.12)) +
  expand_limits(y = min(pca$PC2) - 0.4) +
  scale_y_continuous(expand = expansion(mult = 0.06)) +
  labs(x = sprintf("PC1 (%.1f%%)", 100 * ve[1]), y = sprintf("PC2 (%.1f%%)", 100 * ve[2])) +
  theme_fig()

## ---- ADMIXTURE
samples <- readLines(file.path(dir, "admixture_samples.txt"))
Ks <- sort(as.integer(str_match(list.files(file.path(dir, "admixture"), "^ad\\.\\d+\\.Q$"), "\\.(\\d+)\\.")[, 2]))
q <- bind_rows(lapply(Ks, function(k) {
  m <- as.matrix(read.table(file.path(dir, "admixture", sprintf("ad.%d.Q", k))))
  rownames(m) <- samples
  m <- m[sample_order, , drop = FALSE]
  ord <- order(apply(m, 2, which.max))          # component order: first genome where it peaks
  m <- m[, ord, drop = FALSE]; colnames(m) <- paste0("C", seq_len(k))
  as_tibble(m, rownames = "sample") |> mutate(K = k)
})) |>
  pivot_longer(starts_with("C"), names_to = "component", values_to = "q") |>
  filter(!is.na(q))
write_tsv(q, paste0(prefix, ".admixture.tsv"))

q <- mutate(q, sample = factor(sample, levels = rev(sample_order)), K = factor(paste0("K = ", K)),
            component = factor(component, levels = rev(paste0("C", 1:5))))
pADM <- ggplot(q, aes(x = q, y = sample, fill = component)) +
  geom_col(width = 0.85) +
  facet_wrap(~ K, nrow = 1) +
  scale_fill_manual(values = setNames(col_k, paste0("C", 1:5)), guide = "none") +
  scale_y_discrete(labels = function(x) sample_expr(x)) +
  scale_x_continuous(breaks = NULL, expand = c(0, 0)) +   # each bar spans 0 to 1
  labs(x = "Ancestry proportion", y = NULL) +
  theme_fig() +
  theme(panel.grid = element_blank(), panel.spacing.x = unit(1.2, "mm"),
        strip.text = element_text(face = "plain", hjust = 0.5, size = 7.5))

saveRDS(list(pca = pPCA, admixture = pADM), paste0(prefix, ".panels.rds"))
fig <- (pPCA | pADM) + plot_layout(widths = c(1, 1.15)) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 11, face = "bold"))
save_fig(fig, prefix, width = 7.5, height = 2.9)
