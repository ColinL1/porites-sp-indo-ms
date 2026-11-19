#!/usr/bin/env Rscript
# Functional category composition from funannotate compare, Porites genomes only, as one heatmap:
#   COG categories (% of COG-assigned proteins), CAZy classes, MEROPS peptidase classes, secreted proteins
#   (SignalP) and transcription factor families (InterPro; families with a max count >= --min-tf).
# Cells are coloured by log2(value / median of the Porites for that row), so rows with different units sit
# side by side; the printed number is the value itself.
# usage: functional_categories.R <funannotate_compare dir> <out prefix> [min TF count, default 10]
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(script_dir, "..", "..", "..", "shared", "plots", "setup.R"))   # names, colours, genome order, shared functions

args <- commandArgs(trailingOnly = TRUE)
fc_dir <- args[1]; prefix <- args[2]
min_tf <- if (length(args) >= 3) as.numeric(args[3]) else 10
fc <- function(...) file.path(fc_dir, ...)

to_long <- function(d, source) d |> select(id, all_of(porites)) |>
  pivot_longer(all_of(porites), names_to = "genome", values_to = "value") |>
  mutate(value = as.numeric(value), source = source)

cog <- read_fc_table(fc("cogs", "COGS.all.results.csv")) |> to_long("COG (% of proteins with a COG)") |>
  group_by(genome) |> mutate(value = 100 * value / sum(value)) |> ungroup() |>
  mutate(id = str_replace(id, "^\\((\\w)\\) ", "\\1  "))

cazy_names <- c(AA = "Auxiliary activities", CBM = "Carbohydrate-binding modules", CE = "Carbohydrate esterases",
                GH = "Glycoside hydrolases", GT = "Glycosyltransferases", PL = "Polysaccharide lyases")
cazy <- read_fc_table(fc("cazy", "CAZyme.summary.results.csv")) |> to_long("CAZy classes") |>
  mutate(id = paste0(id, "  ", cazy_names[id]))

merops_names <- c(A = "Aspartic", C = "Cysteine", G = "Glutamic", M = "Metallo", P = "Mixed",
                  S = "Serine", T = "Threonine", I = "Inhibitors")
merops <- read_fc_table(fc("merops", "MEROPS.summary.results.csv")) |> to_long("MEROPS peptidases") |>
  mutate(id = paste0(id, "  ", merops_names[id]))

signalp <- read_fc_table(fc("signalp", "signalp.csv")) |> rename(genome = id, value = 2) |>
  filter(genome %in% porites) |>
  mutate(id = "Proteins with a signal peptide", source = "SignalP", value = as.numeric(value))

# the funannotate TF table lists some InterPro families twice under different names (e.g. IPR004827);
# keep the first row per InterPro ID
tf_raw <- read_csv(fc("tfs", "transcription_factor_counts.csv"), show_col_types = FALSE,
                   name_repair = "minimal") |>
  select(-1) |> distinct(InterPro, .keep_all = TRUE) |> mutate(id = paste0(Description, "  ", InterPro))
tf <- tf_raw |> to_long("Transcription factor families") |>
  group_by(id) |> filter(max(value) >= min_tf) |> ungroup()

all <- bind_rows(cog, cazy, merops, signalp, tf) |>
  group_by(source, id) |>
  mutate(rel = log2((value + 0.5) / (median(value) + 0.5))) |> ungroup() |>
  mutate(genome = factor(genome, levels = porites),
         source = factor(source, levels = unique(c(cog$source, cazy$source, merops$source,
                                                   signalp$source, tf$source))),
         txt = ifelse(str_starts(source, "COG"), sprintf("%.1f", value),
                      formatC(value, format = "d", big.mark = ",")))
# rows top to bottom: COG in letter order, the other blocks by decreasing median value
row_levels <- all |> group_by(source, id) |> summarise(size = median(value), .groups = "drop") |>
  mutate(key = ifelse(str_starts(source, "COG"), rank(id), -size)) |>
  arrange(source, key) |> pull(id) |> as.character()
all <- all |> mutate(id = factor(id, levels = row_levels))
write_tsv(all |> select(source, category = id, genome, value, log2_vs_median = rel),
          paste0(prefix, ".tsv"))

lim <- 1.5
fig <- ggplot(all, aes(x = genome, y = fct_rev(id), fill = pmax(pmin(rel, lim), -lim))) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = txt), size = 2.2, colour = ink) +
  facet_grid(source ~ ., scales = "free_y", space = "free_y", labeller = label_wrap_gen(14)) +
  scale_fill_gradient2(low = col_under, mid = col_mid, high = col_over, midpoint = 0,
                       limits = c(-lim, lim), breaks = c(-1, 0, 1), labels = c("half", "median", "double"),
                       name = label_expr("Relative to\nPorites median")) +
  scale_x_genome() +
  labs(x = NULL, y = NULL) +
  theme_fig() + theme(axis.text.x = element_text(angle = 30, hjust = 1),
                      axis.text.y = element_text(size = 7), panel.grid = element_blank(),
                      strip.text.y = element_text(angle = 0, size = 7.5))

n_rows <- n_distinct(paste(all$source, all$id))
save_fig(fig, prefix, width = 7.5, height = 1.6 + n_rows * 0.155)
message(sprintf("%d rows", n_rows))
