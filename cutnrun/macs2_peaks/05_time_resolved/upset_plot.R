#!/usr/bin/env Rscript
# Usage: Rscript upset_plot.R <multiinter.txt> <output.pdf> <title> [nintersects]
suppressPackageStartupMessages(library(UpSetR))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript upset_plot.R <multiinter.txt> <output.pdf> <title> [nintersects]")
}
input_file  <- args[1]
output_file <- args[2]
title       <- args[3]
nintersects <- if (length(args) >= 4) as.integer(args[4]) else 40

d <- read.table(input_file, header = TRUE, sep = "\t", comment.char = "")
colnames(d)[1] <- sub("^#", "", colnames(d)[1])

set_cols <- setdiff(colnames(d), c("chrom", "start", "end", "num", "list"))
sets_df <- d[, set_cols, drop = FALSE]

pdf(output_file, width = 10, height = 6, onefile = FALSE)
print(upset(
  sets_df,
  sets = rev(set_cols),
  nsets = length(set_cols),
  keep.order = TRUE,
  nintersects = nintersects,
  order.by = "freq",
  mainbar.y.label = paste0(title, ": Peak Intersections"),
  sets.x.label = "Peaks per Condition",
  text.scale = 1.3
))
invisible(dev.off())

cat("Wrote", output_file, "\n")
