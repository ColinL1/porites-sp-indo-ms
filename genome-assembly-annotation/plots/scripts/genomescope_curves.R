# genomescope_curves.R: run GenomeScope2 exactly as the GenomeMaker GENOMESCOPE process does, and also write
# the curves of its linear plot (observed, full model, unique sequence, errors, k-mer peaks) as a table, so that
# kmer_spectrum.R can redraw them. The curves are captured from GenomeScope's own drawing calls (trace on
# lines/points/abline), so they are exactly what GenomeScope plots; nothing of the model is re-implemented.
# Its statistics lines (the plot title) are written per plot as <plot>.title.txt, from its own title() call.
# Optional plot_px ("W" or "WxH"): GenomeScope's plots are drawn on a W x H pixel device (300 dpi) instead of
# its fixed 2000 x 2000, so that its own plot fits the panel shape of a figure; nothing else changes.
# Run inside the pipeline's container (quay.io/biocontainers/genomescope2:2.1.0--py313r44hdfd78af_0):
# usage: Rscript genomescope_curves.R kmer_db.meryl.hist out_dir k ploidy kcov_init [plot_px]
args <- commandArgs(trailingOnly = TRUE)
hist <- args[1]; out <- args[2]; k <- args[3]; p <- args[4]; l <- args[5]
plot_px <- if (length(args) >= 6 && nzchar(args[6])) as.numeric(strsplit(args[6], "x")[[1]]) else NA
if (length(plot_px) == 1) plot_px <- rep(plot_px, 2)
dir.create(out, showWarnings = FALSE, recursive = TRUE)

# as the pipeline: keep the first two columns of the meryl histogram
clean <- file.path(out, "kmer_db.meryl.hist.clean")
h <- read.table(hist)[, 1:2]
write.table(h, clean, sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)

# record every drawing call with the plot file it goes to (GenomeScope opens all its PNGs first and then
# switches between them, so the file is looked up from the active device)
rec <- new.env(); rec$calls <- list(); rec$dev <- list()
rec$px <- plot_px
trace("png", tracer = quote(if (!anyNA(rec$px)) { width <- rec$px[1]; height <- rec$px[2] }),
      exit = quote(rec$dev[[as.character(grDevices::dev.cur())]] <- basename(filename)), print = FALSE)
dev_file <- function() { f <- rec$dev[[as.character(grDevices::dev.cur())]]; if (is.null(f)) "" else f }
# lines(x, ...) / points(x, ...): y is the first unnamed argument in ..., col a named one
grab <- function(type) bquote({
  a <- list(...); un <- a[names(a) %in% c("", NA) | is.null(names(a))]
  xy <- grDevices::xy.coords(x, if ("y" %in% names(a)) a$y else if (length(un)) un[[1]] else NULL)
  rec$calls[[length(rec$calls) + 1]] <- data.frame(file = dev_file(), call = .(type),
    col = if (is.null(a$col)) NA_character_ else as.character(a$col)[1], x = xy$x, y = xy$y)
})
trace("lines", grab("lines"), print = FALSE)
trace("points", grab("points"), print = FALSE)
# title(): GenomeScope's statistics lines (len, uniq, aa, ab, kcov, err, dup, k, p), per plot file
rec$titles <- list()
trace("title", quote(if (is.character(main) && grepl("len:", main[1])) rec$titles[[dev_file()]] <- main[1]),
      print = FALSE)
trace("abline", quote(if (!is.null(v)) rec$calls[[length(rec$calls) + 1]] <-
  data.frame(file = dev_file(), call = "abline", col = NA_character_, x = v, y = NA)), print = FALSE)

# the CLI script, with its command line replaced by the pipeline's arguments
cli <- Sys.which("genomescope2")
src <- readLines(cli)
call <- "arguments <- parser$parse_args()"
if (sum(src == call) != 1) stop("unexpected genomescope2 script: no single '", call, "' line")
src[src == call] <- sprintf("arguments <- parser$parse_args(c('-i', '%s', '-o', '%s', '-k', '%s', '-p', '%s', '-l', '%s'))",
                            clean, out, k, p, l)
eval(parse(text = src), envir = globalenv())

# one text file per plot with its statistics lines, as GenomeScope prints them (leading blank lines dropped)
for (f in names(rec$titles))
  writeLines(sub("^\\n+", "", rec$titles[[f]]), file.path(out, sub("\\.png$", ".title.txt", f)))

d <- do.call(rbind, rec$calls)
d <- d[d$file == "linear_plot.png", ]
role <- c("#56B4E9" = "observed", "black" = "model", "#F0E442" = "unique", "#D55E00" = "errors")
d$curve <- ifelse(d$call == "abline", "peak", unname(role[d$col]))
d <- d[!is.na(d$curve), c("curve", "x", "y")]
write.table(d, file.path(out, "linear_plot.curves.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("curves:", paste(names(table(d$curve)), table(d$curve), collapse = ", "), "\n")
