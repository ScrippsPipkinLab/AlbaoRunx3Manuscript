# preprocessTracks.R
# R script to prepare and pre-process genomic data for genome track plotting

# Load libraries
# Below packages help with genome visualization and annotation
library(Gviz)         # For genome-based plots
library(rtracklayer)  # For importing/exporting track files
library(GenomicFeatures) # For managing transcript-related data
library(txdbmaker)    # For creating TxDb from GTF/GFF

# Load Functions (must come after load() above, not before)
source("~/AlbaoRunx3Manuscript/tracks/scripts/plotTracks_functions.R")

# Genome
genome <- "mm10"

############### GENOME AXIS TRACK ###############
axis_track <- GenomeAxisTrack(labelPos = "above") # Genome-wide axis with tick marks

############### GENE TRACK ###############
# Path to the GTF file
# This GTF is from ENSEMBL.
# Chromosomes in ENSEMBL GTFs do not have the "chr" prefix.
# Thus it is not compatible with ideogram_track which is UCSC-annotated, and uses the chr prefix.
# I therefore parsed the GTF to add the "chr" prefix to the chromosome names before loading here!
# Moreover, this version has a modified gene_id attribute to contain gene symbol rather than ENSEMBL gene ID.
# This is to facilitate plotting of gene symbols in the gene track.
gtf_file <- "~/AlbaoRunx3Manuscript/tracks/source_data/genome/genesUCSC-Symbol.gtf"  # Path to GTF file containing gene annotations


# Import the full GTF, then also make a coding-only subset that drops non-coding
# transcripts (lncRNA, miRNA, pseudogenes, etc.), keeping only protein_coding ones
gtf_data <- import(gtf_file, format = "gtf")
gtf_data_coding <- gtf_data[!is.na(mcols(gtf_data)$transcript_biotype) & mcols(gtf_data)$transcript_biotype == "protein_coding"]

# Convert both to TxDb objects
txdb <- makeTxDbFromGRanges(gtf_data_coding) # Coding-only TxDb
print(txdb)  # Inspect the TxDb object

txdb_full <- makeTxDbFromGRanges(gtf_data) # Full TxDb (all transcripts, including non-coding)
print(txdb_full)  # Inspect the full TxDb object

# Create gene_data object
# Keep gene_data to inform gene highlighting later on
gene_data <- genes(txdb) # Extract gene-level information
mcols(gene_data)$start <- start(gene_data) # Store start positions in metadata
mcols(gene_data)$end <- end(gene_data)
mcols(gene_data)$width <- width(gene_data)

# Full version of gene_data (includes genes made up entirely of non-coding transcripts)
gene_data_full <- genes(txdb_full)
mcols(gene_data_full)$start <- start(gene_data_full)
mcols(gene_data_full)$end <- end(gene_data_full)
mcols(gene_data_full)$width <- width(gene_data_full)

# Create a GeneRegionTrack from TxDb
gene_track <- GeneRegionTrack( # Create a track displaying gene regions
    txdb,
    genome = genome,
    transcriptAnnotation = "gene",
    stacking = "squish",
)

# Full version of gene_track (includes non-coding transcripts)
gene_track_full <- GeneRegionTrack(
    txdb_full,
    genome = genome,
    transcriptAnnotation = "gene",
    stacking = "squish",
)

gc(full = T)

############### MULTIPLE BIGWIG FILES ###############
# List of BigWig files
# This defined code-friendly names for the BigWig files to be used in plotting.
common_path <- "~/AlbaoRunx3Manuscript/tracks/source_data/"
bigwig_files <- c(  "atac.naive" = paste0(common_path, "atac/", "naive.mRp.clN.bigWig"),
                    "atac.d5.shCon" = paste0(common_path, "atac/", "D5_shCD19.mRp.clN.bigWig"),
                    "atac.d5.shR3"  = paste0(common_path, "atac/", "D5_shRunx3.mRp.clN.bigWig"),
                    "atac.d8.gfpPos.shCon" = paste0(common_path, "atac/", "GFPpos_shCD19.mRp.clN.bigWig"),
                    "atac.d8.gfpPos.shR3"  = paste0(common_path, "atac/", "GFPpos_shRunx3.mRp.clN.bigWig"),
                    "atac.d8.gfpNeg.shCon" = paste0(common_path, "atac/", "GFPneg_shCD19.mRp.clN.bigWig"),
                    "atac.d8.gfpNeg.shR3"  = paste0(common_path, "atac/", "GFPneg_shRunx3.mRp.clN.bigWig"),
                    "atac.d8.cx3cr1Pos.shCon" = paste0(common_path, "atac/", "CX3CR1pos_shCD19.mRp.clN.bigWig"),
                    "atac.d8.cx3cr1Pos.shR3"  = paste0(common_path, "atac/", "CX3CR1pos_shRunx3.mRp.clN.bigWig"),
                    "K4me3.shCon"   = paste0(common_path, "cutnrun1/", "shCd19_H3K4me3_R1.bigWig"),
                    "K4me3.shR3"    = paste0(common_path, "cutnrun1/", "shRunx3_H3K4me3_R1.bigWig"),
                    "K36me3.shCon"  = paste0(common_path, "cutnrun1/", "shCd19_H3K36me3_R1.bigWig"),
                    "K36me3.shR3"   = paste0(common_path, "cutnrun1/", "shRunx3_H3K36me3_R1.bigWig"),
                    "Runx3.shCon"   = paste0(common_path, "cutnrun/", "shCd19_Runx3_log2.bigWig"),
                    "Runx3.shR3"    = paste0(common_path, "cutnrun/", "shRunx3_Runx3_log2.bigWig"),
                    "Runx1.shCon"   = paste0(common_path, "cutnrun/", "shCd19_Runx1_log2.bigWig"),
                    "Runx1.shR3"    = paste0(common_path, "cutnrun/", "shRunx3_Runx1_log2.bigWig"),
                    "Runx3.mem"     = paste0(common_path, "cutnrun/", "memory_Runx3_log2.bigWig"),
                    "Runx3.early"   = paste0(common_path, "cutnrun/", "early_Runx3_log2.bigWig"),
                    "Runx3.late"    = paste0(common_path, "cutnrun/", "late_Runx3_log2.bigWig"),
                    "Runx3.terminal"= paste0(common_path, "cutnrun/", "terminal_Runx3_log2.bigWig"),
                    "Runx1.mem"     = paste0(common_path, "cutnrun/", "memory_Runx1_log2.bigWig"),
                    "Runx1.early"   = paste0(common_path, "cutnrun/", "early_Runx1_log2.bigWig"),
                    "Runx1.late"    = paste0(common_path, "cutnrun/", "late_Runx1_log2.bigWig"),
                    "Runx1.terminal"= paste0(common_path, "cutnrun/", "terminal_Runx1_log2.bigWig")
                )

# Make a Pretty Names Vector (this is for plotting!)
pretty_names <- c(  "naive",
                    "D5 shCd19",
                    "D5 shRunx3",
                    "GFPpos shCd19",
                    "GFPpos shRunx3",
                    "GFPneg shCd19",
                    "GFPneg shRunx3",
                    "CX3CR1pos shCd19",
                    "CX3CR1pos shRunx3",
                    "H3K4me3 shCd19",
                    "H3K4me3 shRunx3",
                    "H3K36me3 shCd19",
                    "H3K36me3 shRunx3",
                    "Runx3 shCd19",
                    "Runx3 shRunx3",
                    "Runx1 shCd19",
                    "Runx1 shRunx3",
                    "Runx3 memory",
                    "Runx3 early",
                    "Runx3 late",
                    "Runx3 terminal",
                    "Runx1 memory",
                    "Runx1 early",
                    "Runx1 late",
                    "Runx1 terminal"
                )
names(pretty_names) <- names(bigwig_files)


# Create a list to store the BigWig tracks
bigwig_data <- list()

# Import BigWig files
# Use bigwig_files vector as a key to the list of BigWig tracks (bigwig_data)
for(i in seq_along(bigwig_files)){
    # Path to the BigWig file
    bigwig_id <- names(bigwig_files[i])
    # Import the BigWig file
    bigwig_data[[bigwig_id]] <- import(bigwig_files[i], format = "BigWig")
    # Prepend "chr" to the chromosome names

    # Prepending chr is important because ideogram tracks are UCSC-annotated and chromsomes have the chr prefix
    if(! bigwig_id %in% c("mock", "shRunx3", "OE")){ # These are scrna-seq BigWigs, which are UCSC-annotated and already have the chr prefix
        seqlevels(bigwig_data[[bigwig_id]]) <- paste0("chr", seqlevels(bigwig_data[[bigwig_id]]))
    }  # These are scrna-seq BigWigs, which are UCSC-annotated and already have the chr prefix

    
    # Remove X and Y chromosomes
    bigwig_data[[bigwig_id]] <- keepSeqlevels(bigwig_data[[bigwig_id]], paste0("chr", c(1:19)), pruning.mode = "coarse")

    # Print first few items of the imported BigWig data for verification
    print(paste("Imported BigWig:", bigwig_id))
    print(head(bigwig_data[[bigwig_id]]))

    gc(full = T)
}

# if rds directory does not exist, create it
if (!dir.exists("rds")) {
    dir.create("rds")
}

# Save data
save.image("rds/plotTracks_essential.RData")
saveRDS(bigwig_data, "rds/bigwig_data.rds")
saveRDS(gene_data, "rds/gene_data.rds")
saveRDS(gene_track, "rds/gene_track.rds")
saveDb(txdb, "rds/txdb.rds")
saveRDS(gene_data_full, "rds/gene_data_full.rds")
saveRDS(gene_track_full, "rds/gene_track_full.rds")
saveDb(txdb_full, "rds/txdb_full.rds")