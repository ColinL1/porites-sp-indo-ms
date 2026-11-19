#!/usr/bin/env Rscript
# GO enrichment from funannotate compare (goatools find_enrichment.py: each genome's genes vs the pooled
# genes of all genomes, Fisher's exact test, resampling FDR, alpha 0.001). Only significant terms are in
# the goatools output, so an empty cell means "not significant", not "no change".
#   A. number of over- (e) and under-represented (p, "purified") terms per genome and GO namespace
#   B. dot matrix genome x term: colour = log2 fold enrichment (ratio_in_study / ratio_in_pop),
#      size = genes annotated with the term; terms shallower than min depth are left out (generic parents),
#      and redundant parents are collapsed: a term is dropped when an is_a descendant (the relation goatools
#      propagates over) is also shown and is significant in the same direction in every genome where the
#      parent is. Parents that add information (significant in a genome where no shown child is) are kept.
# usage: go_enrichment.R <funannotate_compare dir> <out prefix> <go.obo> [min depth, default 3]
#        go.obo: the release funannotate used (data/go_2024-11-03.obo, releases/2024-11-03)
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(script_dir, "..", "..", "..", "shared", "plots", "setup.R"))   # names, colours, genome order, shared functions

args <- commandArgs(trailingOnly = TRUE)
fc_dir <- args[1]; prefix <- args[2]; obo_file <- args[3]
min_depth <- if (length(args) >= 4) as.integer(args[4]) else 3L

# genomes plotted, new assembly first. Acropora millepora is left out of the figure (not informative for
# a Porites genome description) but its genes are still part of the goatools background pool.
genomes <- porites
n_pool <- length(list.files(file.path(fc_dir, "go_enrichment"), pattern = "\\.go\\.enrichment\\.txt$"))
ns_names <- c(BP = "Biological process", MF = "Molecular function", CC = "Cellular component")

files <- file.path(fc_dir, "go_enrichment", paste0(genomes, ".go.enrichment.txt"))
stopifnot(all(file.exists(files)))
frac <- function(r) { p <- str_split_fixed(r, "/", 2); as.numeric(p[, 1]) / as.numeric(p[, 2]) }

go <- bind_rows(lapply(seq_along(genomes), function(i) {
  read_tsv(files[i], col_types = cols(.default = "c"), comment = "") |>
    rename(GO = `# GO`) |>
    mutate(genome = genomes[i])
})) |>
  mutate(depth = as.integer(depth), study_count = as.integer(study_count),
         p_uncorrected = as.numeric(p_uncorrected), p_fdr = as.numeric(p_fdr),
         fold = frac(ratio_in_study) / frac(ratio_in_pop), log2fe = log2(fold),
         direction = recode(enrichment, e = "over", p = "under"),
         NS = factor(NS, levels = names(ns_names)),
         genome = factor(genome, levels = genomes))

write_tsv(go |> select(genome, GO, NS, name, depth, direction, ratio_in_study, ratio_in_pop,
                       fold, study_count, p_uncorrected, p_fdr) |> arrange(genome, NS, desc(abs(log2(fold)))),
          paste0(prefix, ".terms.tsv"))


## A. counts per genome / namespace / direction (under drawn to the left)
counts <- go |> count(genome, NS, direction) |>
  complete(genome, NS, direction = c("over", "under"), fill = list(n = 0)) |>
  mutate(signed = ifelse(direction == "under", -n, n))
lim_a <- max(abs(counts$signed))
pA <- ggplot(counts, aes(y = fct_rev(genome), x = signed, fill = direction)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 0, colour = ink_muted, linewidth = 0.3) +
  geom_text(data = filter(counts, n > 0), aes(label = n, hjust = ifelse(signed < 0, 1.2, -0.2)),
            size = 2.4, colour = ink_muted) +
  facet_wrap(~NS, nrow = 1, labeller = as_labeller(setNames(str_replace(ns_names, " ", "\n"), names(ns_names)))) +
  scale_fill_manual(values = c(over = col_over, under = col_under),
                    labels = c(over = "over-represented", under = "under-represented"), name = NULL) +
  scale_x_continuous(limits = c(-lim_a, lim_a) * 1.25, labels = abs) +
  scale_y_genome() +
  labs(x = "significant GO terms (FDR < 0.001)", y = NULL) +
  theme_fig() + theme(legend.position = "top", legend.justification = "left",
                      panel.grid.major.y = element_blank())

## B. dot matrix of the more specific terms
mat <- go |> filter(depth >= min_depth)

# is_a parents from the obo, then all is_a ancestors of each shown term
obo <- readLines(obo_file)
term_start <- which(obo == "[Term]")
blk <- findInterval(seq_along(obo), c(term_start, which(obo == "[Typedef]")[1]))
id_line <- grepl("^id: GO:", obo); isa_line <- grepl("^is_a: GO:", obo)
blk_id <- setNames(sub("^id: ", "", obo[id_line]), blk[id_line])
is_a <- tibble(child = blk_id[as.character(blk[isa_line])], parent = substr(obo[isa_line], 7, 16))
parents_of <- split(is_a$parent, is_a$child)
ancestors <- function(id) {
  out <- character(); todo <- parents_of[[id]]
  while (length(todo)) { out <- union(out, todo); todo <- setdiff(unlist(parents_of[todo]), out) }
  out
}
shown <- unique(mat$GO)
anc <- lapply(setNames(shown, shown), ancestors)
# signature of a term: the genomes and directions where it is significant
sig <- mat |> group_by(GO) |> summarise(sig = list(paste(genome, direction)), .groups = "drop")
sig <- setNames(sig$sig, sig$GO)
redundant <- vapply(shown, function(p) {
  kids <- shown[vapply(shown, function(k) p %in% anc[[k]], logical(1))]
  length(kids) > 0 && all(sig[[p]] %in% unlist(sig[kids]))
}, logical(1))
message(sprintf("collapsed %d redundant parent terms: %s", sum(redundant),
                paste(unique(mat$name[mat$GO %in% shown[redundant]]), collapse = "; ")))
mat <- mat |> filter(!GO %in% shown[redundant])
# order rows within each namespace by clustering their log2FE profiles (not significant = 0)
row_order <- mat |> group_by(NS) |> group_modify(function(d, key) {
  w <- d |> select(name, genome, log2fe) |>
    pivot_wider(names_from = genome, values_from = log2fe, values_fill = 0)
  m <- as.matrix(w[, -1]); rownames(m) <- w$name
  ord <- if (nrow(m) > 2) hclust(dist(m))$order else seq_len(nrow(m))
  tibble(name = rownames(m)[ord], row = seq_along(ord))
}) |> ungroup()
mat <- mat |> left_join(row_order, by = c("NS", "name")) |>
  mutate(term = fct_reorder(str_wrap(name, 48), row))
lim_b <- max(abs(mat$log2fe))
pB <- ggplot(mat, aes(x = genome, y = term)) +
  annotate("rect", xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf, fill = band) +
  geom_point(aes(size = study_count, fill = log2fe), shape = 21, colour = ink_muted, stroke = 0.25) +
  facet_grid(NS ~ ., scales = "free_y", space = "free_y",
             labeller = as_labeller(setNames(str_replace(ns_names, " ", "\n"), names(ns_names)))) +
  scale_fill_gradient2(low = col_under, mid = col_mid, high = col_over, midpoint = 0,
                       limits = c(-lim_b, lim_b), name = "log2 fold\nenrichment") +
  scale_size_area(max_size = 5, breaks = c(10, 100, 1000), trans = "sqrt", name = "genes") +
  scale_x_genome(drop = FALSE) +
  labs(x = NULL, y = NULL) +
  theme_fig() + theme(axis.text.x = element_text(angle = 35, hjust = 1),
                    axis.text.y = element_text(size = 7, lineheight = 0.85),
                    strip.text.y = element_text(angle = -90, hjust = 0.5),
                    panel.spacing.y = unit(6, "pt"))

n_rows <- nrow(row_order)
row_h <- 0.2    # inches per term row (room for two-line labels)
fig <- pA / pB + plot_layout(heights = c(2.2, n_rows * row_h)) +
  plot_annotation(tag_levels = "A")

h <- 2.6 + n_rows * row_h + 1.4
save_fig(fig, prefix, width = 7.5, height = h)
saveRDS(list(A = pA, B = pB, n_rows_B = n_rows), paste0(prefix, ".panels.rds"))   # for composite.R
message(sprintf("%d significant terms (%d unique); %d shown in B (depth >= %d)",
                nrow(go), n_distinct(go$GO), n_rows, min_depth))
