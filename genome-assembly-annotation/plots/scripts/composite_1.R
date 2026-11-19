# composite_1.R: genome assembly figure from the tools' own plots, A | B, both the same height.
#   A: GenomeScope2 k-mer spectrum (its transformed linear plot, coverage x frequency, drawn by GenomeScope on a
#      device of the panel's shape). Two versions:
#      with <plot>.title.txt (composite_1a): GenomeScope's title area is cropped and its statistics lines (len,
#        uniq, aa, ab, kcov, err, dup, k, p; the exact strings of its title() call) are set above the plot in
#        the figure font, centred on the plot frame as GenomeScope places them;
#      with "none" (composite_1b): GenomeScope's title area cropped, no header (A is the plot alone; drawn by
#        GenomeScope on a taller device so that A matches B's height without a header).
#   B: BlobToolKit snail plot of the release (viewer, large fonts)
# The figure is 7.49 in wide; the widths are solved so that A (header + plot) and B (square) are equally tall.
# usage: Rscript composite_1.R genomescope_plot.png genomescope_plot.title.txt|none snail.png out_prefix
# env: conda env r_env_DESeq2
args <- commandArgs(trailingOnly = TRUE)
here <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))))
source(file.path(here, "..", "..", "..", "shared", "plots", "setup.R"))
gs_png <- args[1]; gs_title <- args[2]; snail_png <- args[3]; out <- args[4]

no_header <- gs_title == "none"             # plot alone, without the statistics lines
gs_crop_top <- 330                           # px: GenomeScope's title area (the frame starts below it)
fig_w <- 7.49; gap <- 0.25                   # in: figure width, room for tags and margins
hdr_pt <- 8                                  # pt: statistics lines
hdr_lines <- if (no_header) character(0) else trimws(readLines(gs_title))
hdr_in <- if (no_header) 0 else (length(hdr_lines) * 1.25 + 0.4) * hdr_pt / 72   # in: header block above the plot

gs <- png::readPNG(gs_png)
gs <- gs[(gs_crop_top + 1):dim(gs)[1], , , drop = FALSE]
r <- dim(gs)[1] / dim(gs)[2]                 # height / width of the cropped plot
# plot frame (columns of the box drawn by GenomeScope) to centre the header on it, as GenomeScope does
dark <- which(colMeans(gs[, , 1] < 0.5) > 0.5)
frame_mid <- (min(dark) + max(dark)) / 2 / dim(gs)[2]

# A width a_w, B width b_w: a_w + b_w = fig_w - gap and a_w * r + hdr_in = b_w
b_w <- ((fig_w - gap) * r + hdr_in) / (1 + r)
a_w <- fig_w - gap - b_w
hdr_u <- hdr_in / a_w                         # header height in A's data units (x spans 0 to 1)

tag_theme <- theme(plot.tag = element_text(family = fig_font, face = "bold", size = 11, colour = ink),
                   plot.tag.position = "topleft", plot.margin = margin(2, 4, 2, 2))
pa <- ggplot() + annotation_raster(gs, 0, 1, 0, r, interpolate = TRUE)
if (!no_header)
  pa <- pa + annotate("text", x = frame_mid, y = r + hdr_u / 2, label = paste(hdr_lines, collapse = "\n"),
                      family = fig_font, size = hdr_pt / .pt, lineheight = 1.05, colour = ink)
pa <- pa +
  coord_fixed(ratio = 1, xlim = c(0, 1), ylim = c(0, r + hdr_u), expand = FALSE) +
  theme_void() + labs(tag = "A") + tag_theme
snail <- png::readPNG(snail_png)
pb <- ggplot() + annotation_raster(snail, 0, 1, 0, 1, interpolate = TRUE) +
  coord_fixed(ratio = dim(snail)[1] / dim(snail)[2], xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  theme_void() + labs(tag = "B") + tag_theme

fig <- pa + pb + plot_layout(widths = c(a_w, b_w))
fig_h <- b_w + gap
save_fig(fig, out, fig_w, fig_h)
cat(sprintf("composite_1 %.2f x %.2f in (A %.2f in wide, plot %.2f + header %.2f in tall; B %.2f in)\n",
            fig_w, fig_h, a_w, a_w * r, hdr_in, b_w))
