# Load libraries
library(Gviz)
library(dplyr)
library(tidyr)
library(rtracklayer)
library(GenomicFeatures)
library(txdbmaker)

# Set location
setwd("~/AlbaoRunx3Manuscript/tracks/qc")

# Load Pre-Computed Data first -- load() overwrites any object in the current
# environment that shares a name with something saved in the RData, so it must
# run BEFORE sourcing the functions file, or a stale function baked into the
# RData (from whatever was in memory when preprocessTracks.R last called
# save.image()) would silently clobber the current, correct function definitions.
# Skipped if bigwig_data is already in memory (this file is large and slow to
# reload) -- rm(bigwig_data) first, or restart R, to force a fresh reload after
# re-running preprocessTracks.R.
if (!exists("bigwig_data")) {
    load("rds/plotTracks_essential.RData")
}

# Read significantly accessible peaks
daps <- read.csv("../source_data/daps.csv")

# Load Functions (must come after load() above, not before)
source("~/AlbaoRunx3Manuscript/tracks/scripts/plotTracks_functions.R")

# Define genome
genome <- "mm10"

# Define Autoscale Groups
# Each named vector is an autoscale group
# The names of the vector are the sample names
# The values of the vector are the line colors to be used for each sample
cnr_k4me3.samps   <- c( "K4me3.shCon"   = "#4DBBD5",
                        "K4me3.shR3"    = "#F95335"
)
cnr_k36me3.samps <- c(  "K36me3.shCon"  = "#F39B7F",
                        "K36me3.shR3"   = "#F95335"
)
cnr_d5_runx3.samps  <- c(   "Runx3.shCon"   = "#00A087",
                            "Runx3.shR3"    = "#F95335"
)
cnr_d5_runx1.samps  <- c(   "Runx1.shCon"   = "#8491B4",
                            "Runx1.shR3"    = "#F95335"
)
cnr_d8_runx3.samps  <- c(   "Runx3.mem"     = "#00A087",
                            "Runx3.early"   = "#00A087",
                            "Runx3.late"    = "#00A087",
                            "Runx3.terminal"= "#00A087"
)
cnr_d8_runx1.samps  <- c(   "Runx1.mem"     = "#8491B4",
                            "Runx1.early"   = "#8491B4",
                            "Runx1.late"    = "#8491B4",
                            "Runx1.terminal"= "#8491B4"
)
atc_d5.samps <- c(  "atac.d5.shCon" = "#1f0c74",
                    "atac.d5.shR3"  = "#F95335"
)
atc_d8.samps <- c(  "atac.d8.gfpPos.shCon" = "#1f0c74",
                    "atac.d8.gfpPos.shR3"  = "#F95335",
                    "atac.d8.gfpNeg.shCon" = "#1f0c74",
                    "atac.d8.gfpNeg.shR3"  = "#F95335",
                    "atac.d8.cx3cr1Pos.shCon" = "#1f0c74",
                    "atac.d8.cx3cr1Pos.shR3"  = "#F95335"
)

# Combine Autoscale Groups into a list
# This is used to define the plotting order of the autoscale groups
all <- list(cnr_d5_runx3  = cnr_d5_runx3.samps,
            cnr_d5_runx1  = cnr_d5_runx1.samps,
            atc_d8        = atc_d8.samps
        )

# Define panel colors for the title panes of each autoscale group
group_color <- c(   cnr_d5_runx3  = "#00A087",
                    cnr_d5_runx1  = "#8491B4",
                    atc_d8        = "white"
                )

# Define per-group Y-axis scaling factors (1 = no extra scaling)
norm_factor <- c(   cnr_d5_runx3  = 1,
                    cnr_d5_runx1  = 1,
                    cnr_d8_runx3  = 1,
                    cnr_d8_runx1  = 1,
                    atc_d8        = 1
                )

# Define the DataTrack plot shape per autoscale group (passed to Gviz's "type")
group_type <- c(    cnr_d5_runx3  = "polygon",
                    cnr_d5_runx1  = "polygon",
                    atc_d8        = "polygon"
                )

# Define whether to fill zero-coverage gaps per autoscale group
group_densify <- c( cnr_d5_runx3  = FALSE,
                    cnr_d5_runx1  = FALSE,
                    atc_d8        = FALSE
                    )


# Degfine highlight colors for DAPs
dap_highlights <- c(    Runx3opened_GFPpos  = "#faedcb",     # yellow
                        Runx3opened_all     = "#faedcb",     # yellow
                        Runx3closed_GFPpos  = "#c9e4df",     # green
                        Runx3closed_all     = "#c5def2"      # blue
)


# Change display of gene_track and axis_track
# Here I change some display parameters that my helper scripts already preset
# This allowed me to change the look of the final product
# Without having to change the helper scripts
displayPars(gene_track)$background.title <- "transparent"       # Makes gene_track title transparent
# displayPars(gene_track)$collapseTranscripts <- "meta"           # Collapses transcripts of the same gene in gene_track
# displayPars(gene_track)$transcriptAnnotation <- TRUE            # Use ENSEMBL gene IDs
displayPars(gene_track)$showTitle <- FALSE                      # Do not show a title for gene_track
displayPars(gene_track)$size <- 2                               # Expand the size of the gene_track
displayPars(gene_track)$fontcolor.group <- "black"              # Gene label text color
displayPars(gene_track)$fontface.group <- 4                     # Gene label text style (4 = bold italic)
displayPars(gene_track)$col <- "black"                          # Gene model line/exon color
displayPars(gene_track)$collapseTranscripts <- "meta"                          # Gene model line/exon color
displayPars(gene_track)$col.line <- "black"                     # Gene model connecting line color
displayPars(axis_track)$labelPos <- "alternating"               # Makes the DNA position labels alternate in axis_track
displayPars(axis_track)$col <- "black"                          # Axis line/tick color
displayPars(axis_track)$fontcolor <- "black"                    # Axis position label color



############### PEAK TRACK ###############
# Prepare BED file to make the peak annotation track
bed_file1 <- "../source_data/cluster1.clean.noCluster.bed.zero"  # Replace with the path to your BED file
peak_data1 <- import(bed_file1, format = "BED")
bed_file2 <- "../source_data/cluster2.clean.noCluster.bed.zero"  # Replace with the path to your BED file
peak_data2 <- import(bed_file2, format = "BED")
bed_atac <- "../source_data/atac.bed"  # Replace with the path to your BED file
peak_data_atac <- import(bed_atac, format = "BED")

# Set strand to "*"
strand(peak_data1) <- "*"
strand(peak_data2) <- "*"
strand(peak_data_atac) <- "*"
# Prepend "chr" to the chromosome names
seqlevels(peak_data1) <- paste0("chr", seqlevels(peak_data1))
seqlevels(peak_data2) <- paste0("chr", seqlevels(peak_data2))
seqlevels(peak_data_atac) <- paste0("chr", seqlevels(peak_data_atac))
peak_data1 <- keepSeqlevels(peak_data1, paste0("chr", c(1:19)), pruning.mode = "coarse")
peak_data2 <- keepSeqlevels(peak_data2, paste0("chr", c(1:19)), pruning.mode = "coarse")
peak_data_atac <- keepSeqlevels(peak_data_atac, paste0("chr", c(1:19)), pruning.mode = "coarse")

# Create an AnnotationTrack
# Actually create the peak annotation track
peak_track1 <- AnnotationTrack(
    range = peak_data1,
    genome = genome,
    name = "Type 1",
    stacking = "dense",  # Squish peaks
    fill = "darkred",
    col = "darkred",         # Border color
    showFeatureNames = FALSE,  # Hide feature names
    arrowHeadWidth = 0,         # Remove arrows
    showTitle = FALSE
)
# Create an AnnotationTrack
# Actually create the peak annotation track
peak_track2 <- AnnotationTrack(
    range = peak_data2,
    genome = genome,
    name = "Type 2",
    stacking = "dense",  # Squish peaks
    fill = "darkgreen",
    col = "darkgreen",         # Border color
    showFeatureNames = FALSE,  # Hide feature names
    arrowHeadWidth = 0,         # Remove arrows
    showTitle = FALSE
)
peak_track_atac <- AnnotationTrack(
    range = peak_data_atac,
    genome = genome,
    name = "ATAC",
    stacking = "dense",  # Squish peaks
    fill = "black",
    col = "black",         # Border color
    showFeatureNames = FALSE,  # Hide feature names
    arrowHeadWidth = 0,         # Remove arrows
    showTitle = FALSE
)


# Compute peak width as score
# Keep peak_data to inform peak highlighting later on
mcols(peak_data1)$width <- width(peak_data1)
mcols(peak_data2)$width <- width(peak_data2)
mcols(peak_data_atac)$width <- width(peak_data_atac)


# Title options
# this changes how the title of the tracks looks
rotation.title <- 0
cex.title <- 0.8

# Create model tracks
# These are common tracks in all my plots
model_tracks <- list(gene_track, axis_track)


plot_genes <- function(
    gene,
    chromosome,
    start_pos,
    end_pos,
    custom_norm_factor = NULL,
    custom_group_type = NULL,
    custom_group_densify = NULL
){
    if (!is.null(custom_norm_factor)) {
        norm_factor <- custom_norm_factor
    }
    if (!is.null(custom_group_type)) {
        group_type <- custom_group_type
    }
    if (!is.null(custom_group_densify)) {
        group_densify <- custom_group_densify
    }

    result.all <- plot_gene_auto(
        gene = gene,                              # Gene of interest
        genome = genome,                          # Genome build
        bw_data = bigwig_data,                    # List containing bigwig GRanges
        norm_factor = norm_factor,                # Per-group Y-axis scaling factors
        normalization_groups = all,               # BigWig Tracks to plot grouped by normalization groups
        group_panel_colors = group_color,         # Dictionary of colors for the title panes of each normalization group
        pretty_names = pretty_names,              # Pretty names of the samples
        daps = daps,
        gene_data = gene_data,                    # df containing gene locations
        highlight_colors = dap_highlights,
        model_tracks = model_tracks,              # axis_track and gene_track
        peak_data = peak_data_atac,
        type = group_type,                        # DataTrack plot shape per autoscale group
        densify = group_densify                   # Densify zeroes?
    )

    pdf(paste0(gene,".all.pdf"), width = 8.5, height = 11)
    # Plot the tracks
    plotTracks( c(result.all$tracks, peak_track1, peak_track2, peak_track_atac),
                from = result.all$start_pos,
                to = result.all$end_pos,
                chromosome = result.all$chromosome,
                rotation.title = rotation.title,
                cex.title = cex.title)
    dev.off()

}

plot_genes <- function(
    gene,
    chromosome,
    start_pos,
    end_pos,
    custom_norm_factor = NULL,
    custom_group_type = NULL,
    custom_group_densify = NULL
){
    if (!is.null(custom_norm_factor)) {
        norm_factor <- custom_norm_factor
    }
    if (!is.null(custom_group_type)) {
        group_type <- custom_group_type
    }
    if (!is.null(custom_group_densify)) {
        group_densify <- custom_group_densify
    }

    result.all <- plot_gene_manual(
        gene = gene,                              # Gene of interest
        genome = genome,                          # Genome build
        bw_data = bigwig_data,                    # List containing bigwig GRanges
        chromosome = chromosome,                  # Chromosome to plot
        start_pos = start_pos,                    # Start position of region to plot
        end_pos = end_pos,                        # End position of region to plot
        norm_factor = norm_factor,                # Per-group Y-axis scaling factors
        normalization_groups = all,               # BigWig Tracks to plot grouped by normalization groups
        group_panel_colors = group_color,         # Dictionary of colors for the title panes of each normalization group
        pretty_names = pretty_names,              # Pretty names of the samples
        daps = daps,
        gene_data = gene_data,                    # df containing gene locations
        highlight_colors = dap_highlights,
        model_tracks = model_tracks,              # axis_track and gene_track
        peak_data = peak_data_atac,
        type = group_type,                        # DataTrack plot shape per autoscale group
        densify = group_densify                   # Densify zeroes?
    )


    pdf(paste0(gene,".manual.pdf"), width = 5, height = 7)
    # Plot the tracks
    plotTracks( c(result.all$tracks, peak_track1, peak_track2, peak_track_atac),
                from = result.all$start_pos,
                to = result.all$end_pos,
                chromosome = result.all$chromosome,
                rotation.title = rotation.title,
                cex.title = cex.title)
    dev.off()

}


# plot_genes(gene = "Cxcr6", chromosome = "chr9", start_pos = 123790000, end_pos = 123854000)

plot_genes(gene = "Tcf7", chromosome = "chr11", start_pos = 52210000, end_pos = 52340000)
plot_genes(gene = "Rora", chromosome = "chr9", start_pos = 68550000, end_pos = 69410000)
