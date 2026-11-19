#!/usr/bin/env Rscript
# Main composite of the introgression tests, sized for a portrait US Letter page (7.5 in wide,
# <= 10 in tall; heights follow the content):
#   A. ancestry painting along the genome, full width
#   below, side by side (1/3 : 2/3):
#   B. Twisst T1 - T2 weight difference along the genome (P. cylindrica outgroup), turned 90 degrees
#      (scaffolds down the y axis)
#   C. Twisst mean topology weights as stacked bars: per outgroup at default settings, above the
#      robustness grid (P. cylindrica) by representation and window size
#   keys below their panels
# The panels are the plot objects saved by ancestry_painting.R and twisst.R (.rds); C is joined as in
# twisst.R.
# usage: composite.R <out prefix> <ancestry .rds> <twisst .panels.rds>
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(script_dir, "setup_introgression.R"))

args <- commandArgs(trailingOnly = TRUE)
prefix <- args[1]
pA <- readRDS(args[2]); tw <- readRDS(args[3])

page_w <- 7.5
base   <- 7     # pt; text size across the composite
small <- theme(text = element_text(size = base), axis.text = element_text(size = base - 1),
               legend.text = element_text(size = base - 0.5), legend.title = element_text(size = base),
               strip.text = element_text(size = base),
               plot.tag = element_text(size = base + 3, face = "bold"))

# C: stacked bars as in twisst.R (outgroups above the robustness grid), topology key centred below them (title on
# top, entries in one row); text size set before wrapping (a wrapped panel does not receive the
# composite-wide theme)
key_below <- theme(legend.position = "bottom", legend.justification = "center", legend.title.position = "top",
                   legend.key.spacing.x = unit(2.5, "mm"), legend.text = element_text(margin = margin(l = 1.2, unit = "mm")),
                   legend.key.size = unit(3, "mm"), legend.margin = margin(0, 0, 0, 0),
                   legend.box.spacing = unit(1.5, "mm"))
pC <- (tw$outgroup_bars / tw$grid_bars) + plot_layout(heights = tw$bars_heights, guides = "collect") &
  small & key_below & theme(strip.text = element_text(size = base - 0.5), legend.text = element_text(size = base - 0.5))
pC <- wrap_elements(full = pC)

# heights (in) per row, from the content: A = legend + panel + axis; the B | C row = the bars of C
# with its keys (B's genome axis is continuous, so it takes the height of C without stretching anything)
h_A <- 1.55; h_BC <- 3.3
H <- h_A + h_BC
stopifnot(H + 0.2 <= 10)

row_BC <- (tw$along_vertical | pC) + plot_layout(widths = c(1, 2))
fig <- (pA / row_BC) + plot_layout(heights = c(h_A, h_BC)) +
  plot_annotation(tag_levels = "A") & small
save_fig(fig, prefix, width = page_w, height = H + 0.2)
