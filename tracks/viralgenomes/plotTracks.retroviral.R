# Load libraries
library(Gviz)
library(rtracklayer)
library(GenomicFeatures)
library(txdbmaker)

# Set location
setwd("~/AlbaoRunx3Manuscript/tracks/viralgenomes")

# Define gem
gem <- "gem1"

# Define retrovirus
rv <- "shRunx3"

# Pretty viral chromosome names
viralchrom_names <- 
    c(  "pMIA_Empty"    = "pMIA_Empty",
        "pMIA_Runx3"    = "pMIA_Runx3",
        "shCd19"        = "LMPD_shCd19",
        "shRunx3"       = "LMPD_shRunx3"
    ) 


# Load Pre-Computed Data first -- load() overwrites any object in the current
# environment that shares a name with something saved in the RData, so it must
# run BEFORE sourcing the functions file, or a stale function baked into the
# RData (from whatever was in memory when preprocessTracks.R last called
# save.image()) would silently clobber the current, correct function definitions.
# Skipped if bigwig_data is already in memory (this file is large and slow to
# reload) -- rm(bigwig_data) first, or restart R, to force a fresh reload after
# re-running preprocessTracks.R.
if (!exists("bigwig_data")) {
    load(paste0("rds/", gem, "_", rv, "_plotTracks_essential.RData"))
}

# Source helper functions
source("/home/dalbao/AlbaoRunx3Manuscript/tracks/scripts/plotTracks_functions.R") # Load custom functions for plotting tracks

# Define genome
genome <- "mm10"

options(ucscChromosomeNames = FALSE)


# Define Autoscale Groups
# Each named vector is an autoscale group
# The names of the vector are the sample names
# The values of the vector are the line colors to be used for each sample
scrna_seq.samps   <- c( "naive"         = "#4A4D4C",
                        "day5"          = "#7E6148",
                        "day8"          = "#8491B4",
                        "mock"          = "#808180",
                        "shRunx3"       = "#DC0000",
                        "OE"            = "#00A087",
                        "dAD"           = "#3C5488",
                        "dID"           = "#F39B7F",
                        "dVWRPY"        = "#4DBBD5"
)



# Combine Autoscale Groups into a list
# This is used to define the plotting order of the autoscale groups
all <- list(scrna_seq  = scrna_seq.samps)

# Define panel colors for the title panes of each autoscale group
group_color <- c(   scrna_seq  = "white")

# Define per-group Y-axis scaling factors (1 = no extra scaling)
norm_factor <- c(   scrna_seq  = 1)

# Define the DataTrack plot shape per autoscale group (passed to Gviz's "type")
group_type <- c(scrna_seq  = "polygon")

# Define whether to fill zero-coverage gaps per autoscale group
group_densify <- c( scrna_seq  = FALSE)

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
displayPars(gene_track)$col.line <- "black"                     # Gene model connecting line color
displayPars(axis_track)$labelPos <- "revAlternating"               # Makes the DNA position labels alternate in axis_track
displayPars(axis_track)$col <- "black"                          # Axis line/tick color
displayPars(axis_track)$fontcolor <- "black"                    # Axis position label color



############### PEAK TRACK ###############
# Prepare BED file to make the peak annotation track
bed_file <- "../source_data/all.clean.noCluster.bed.zero"  # Replace with the path to your BED file
peak_data <- import(bed_file, format = "BED")

# Set strand to "*"
strand(peak_data) <- "*"
# Prepend "chr" to the chromosome names
seqlevels(peak_data) <- paste0("chr", seqlevels(peak_data))
peak_data <- keepSeqlevels(peak_data, paste0("chr", c(1:19)), pruning.mode = "coarse")

# Create an AnnotationTrack
# Actually create the peak annotation track
peak_track <- AnnotationTrack(
    range = peak_data,
    genome = genome,
    name = "Peaks",
    stacking = "dense",  # Squish peaks
    fill = "black",
    col = "black",         # Border color
    showFeatureNames = FALSE,  # Hide feature names
    arrowHeadWidth = 0,         # Remove arrows
    showTitle = FALSE
)

# Compute peak width as score
# Keep peak_data to inform peak highlighting later on
mcols(peak_data)$width <- width(peak_data)



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

    result.all <- plot_gene_manual(
        gene = gene,                              # Gene of interest
        genome = genome,                          # Genome build
        chromosome = chromosome,                  # Chromosome to plot
        start_pos = start_pos,                    # Start position of region to plot
        end_pos = end_pos,                        # End position of region to plot
        bw_data = bigwig_data,                    # List containing bigwig GRanges
        norm_factor = norm_factor,                # Per-group Y-axis scaling factors
        normalization_groups = all,               # BigWig Tracks to plot grouped by normalization groups
        group_panel_colors = group_color,         # Dictionary of colors for the title panes of each normalization group
        pretty_names = pretty_names,              # Pretty names of the samples
        gene_data = gene_data,                    # df containing gene locations
        model_tracks = model_tracks,              # axis_track and gene_track
        type = group_type,                        # DataTrack plot shape per autoscale group
        densify = group_densify,                  # Densify zeroes?
        includeIdeogram = FALSE
    )

    result.tracks <- plot_gene_manual(
        gene = gene,                              # Gene of interest
        genome = genome,                          # Genome build
        chromosome = chromosome,                  # Chromosome to plot
        start_pos = start_pos,                    # Start position of region to plot
        end_pos = end_pos,                        # End position of region to plot
        bw_data = bigwig_data,                    # List containing bigwig GRanges
        norm_factor = norm_factor,                # Per-group Y-axis scaling factors
        normalization_groups = all,               # BigWig Tracks to plot grouped by normalization groups
        group_panel_colors = group_color,         # Dictionary of colors for the title panes of each normalization group
        pretty_names = pretty_names,              # Pretty names of the samples
        gene_data = gene_data,                    # df containing gene locations
        model_tracks = list(),                    # axis_track and gene_track
        type = group_type,                        # DataTrack plot shape per autoscale group
        densify = group_densify,                  # Densify zeroes?
        includeIdeogram = FALSE
    )

    pdf(paste0(gene,".all.pdf"), width = 7, height = 7)
    # Plot the tracks
    plotTracks( c(result.all$tracks, peak_track),
                from = result.all$start_pos,
                to = result.all$end_pos,
                chromosome = result.all$chromosome,
                rotation.title = rotation.title,
                cex.title = cex.title)
    dev.off()

    pdf(paste0(gene,".tracks.pdf"), width = 7, height = 6)
    # Plot the tracks
    plotTracks( c(axis_track, result.tracks$tracks, peak_track),
                from = result.tracks$start_pos,
                to = result.tracks$end_pos,
                chromosome = result.tracks$chromosome,
                rotation.title = rotation.title,
                cex.title = cex.title)
    dev.off()


}

plot_genes(gene = "mir30loop", chromosome = "LMPD_shRunx3", start_pos = -500, end_pos = 3600)
plot_genes(gene = "Runx3", chromosome = "chr4", start_pos = 135106309, end_pos = 135192327)
plot_genes(gene = "Runx1", chromosome = "chr16", start_pos = 92545295, end_pos = 92882320)


