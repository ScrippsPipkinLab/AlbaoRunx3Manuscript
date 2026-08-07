# preprocessTracks.R
# R script to prepare and pre-process genomic data for genome track plotting

# Load libraries
# Below packages help with genome visualization and annotation
library(Gviz)         # For genome-based plots
library(rtracklayer)  # For importing/exporting track files
library(GenomicFeatures) # For managing transcript-related data
library(txdbmaker)    # For creating TxDb from GTF/GFF

# Source helper functions
source("/home/dalbao/AlbaoRunx3Manuscript/tracks/scripts/plotTracks_functions.R") # Load custom functions for plotting tracks

# The viral contig (e.g. "pMIA_Runx3") is not a real UCSC chromosome, and Gviz
# validates chromosome names against UCSC style whenever genome is a recognized
# UCSC assembly (mm10 here). Disable that check so the raw viral contig name
# can be used as-is.
options(ucscChromosomeNames = FALSE)

# Define gem
gem <- "gem1"

# Define retrovirus
rv <- "pMIA_Runx3"

# Pretty viral chromosome names
viralchrom_names <- 
    c(  "pMIA_Empty"    = "pMIA_Empty",
        "pMIA_Runx3"    = "pMIA_Runx3",
        "shCd19"        = "LMPD_shCd19",
        "shRunx3"       = "LMPD_shRunx3"
    ) 

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

# Additional Retroviral GTFs
viral_gtf <- paste0("~/AlbaoRunx3Manuscript/tracks/source_data/viralgenomes/", rv, ".gtf")

# Import the full GTF and fold in the viral construct (as chr<rv>, tagged mm10),
# then also make a coding-only subset that drops non-coding
# transcripts (lncRNA, miRNA, pseudogenes, etc.), keeping only protein_coding ones
gtf_data <- import(gtf_file, format = "gtf")

# Combine host and viral GTFs
gtf_data <- combine_host_viral_gtf(gtf_data, viral_gtf, genome_tag = "mm10", add_chr_prefix = FALSE)

# Select only protein-coding transcripts for the coding-only TxDb
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
common_path <- "~/AlbaoRunx3Manuscript/tracks/source_data/viralgenomes"
bigwig_files <- c(  "mock"          = paste(common_path, paste0(gem, "_", rv), "Base.bw", sep = "/"),
                    "shRunx3"       = paste(common_path, paste0(gem, "_", rv), "Null.bw", sep = "/"),
                    "OE"            = paste(common_path, paste0(gem, "_", rv), "WT.bw", sep = "/"),
                    "dAD"           = paste(common_path, paste0(gem, "_", rv), "dAD.bw", sep = "/"),
                    "dID"           = paste(common_path, paste0(gem, "_", rv), "dID.bw", sep = "/"),
                    "dVWRPY"        = paste(common_path, paste0(gem, "_", rv), "dVWRPY.bw", sep = "/"),
                    "day5"          = paste(common_path, paste0(gem, "_", rv), "d5.bw", sep = "/"),
                    "day8"          = paste(common_path, paste0(gem, "_", rv), "d8.bw", sep = "/"),
                    "naive"         = paste(common_path, paste0(gem, "_", rv), "Naive.bw", sep = "/")
                )

# Make a Pretty Names Vector (this is for plotting!)
pretty_names <- c(  "mock",
                    "shRunx3",
                    "RUNX3-OE",
                    "dAD",
                    "dID",
                    "dVWRPY",
                    "day 5",
                    "day 8",
                    "naive"
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
    # if(! bigwig_id %in% c("mock", "shRunx3", "OE")){ # These are scrna-seq BigWigs, which are UCSC-annotated and already have the chr prefix
    #     seqlevels(bigwig_data[[bigwig_id]]) <- paste0("chr", seqlevels(bigwig_data[[bigwig_id]]))
    # }  # These are scrna-seq BigWigs, which are UCSC-annotated and already have the chr prefix

    
    # Remove X and Y chromosomes
    # Host chromosomes carry the "chr" prefix already; the viral contig does not (see viralchrom_names).
    bigwig_data[[bigwig_id]] <- keepSeqlevels(bigwig_data[[bigwig_id]], c(paste0("chr", 1:19), viralchrom_names[rv]), pruning.mode = "coarse")

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
save.image(paste0("rds/", gem, "_", rv, "_plotTracks_essential.RData"))
saveRDS(bigwig_data, paste0("rds/", gem, "_", rv, "_bigwig_data.rds"))
saveRDS(gene_data, paste0("rds/", gem, "_", rv, "_gene_data.rds"))
saveRDS(gene_track, paste0("rds/", gem, "_", rv, "_gene_track.rds"))
saveDb(txdb, paste0("rds/", gem, "_", rv, "_txdb.rds"))
saveRDS(gene_data_full, paste0("rds/", gem, "_", rv, "_gene_data_full.rds"))
saveRDS(gene_track_full, paste0("rds/", gem, "_", rv, "_gene_track_full.rds"))
saveDb(txdb_full, paste0("rds/", gem, "_", rv, "_txdb_full.rds"))