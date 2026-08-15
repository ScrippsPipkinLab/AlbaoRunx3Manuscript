# plotTracks_scripts.R

# This script contains functions used to prepare ATAC-seq data for plotting

# Bigwig/bedGraph-derived GRanges often omit zero-coverage regions entirely --
# they're absent, not present as zero-valued rows. type="polygon" (and "l",
# "smooth", etc.) draws straight segments between consecutive intervals
# regardless of genomic distance, so a gap gets linearly interpolated instead of
# dropping to baseline. Anchoring each gap's edges with a narrow (1bp) score=0
# point -- rather than one wide score=0 interval spanning the whole gap -- makes
# the interpolation between the two zero edges trivially flat, since Gviz would
# otherwise treat one wide interval as a single point and draw a smooth V/triangle
# dipping down to it and back up, instead of a flat baseline with sharp edges.
#
# The same collapse happens for wide intervals that are already IN the data
# (not gaps at all): Gviz's position() for a DataTrack is rowMeans(start, end)
# per range -- one x per range, regardless of width -- so any run the bigWig
# encodes as a single wide interval (its native RLE compression does this for
# every flat run, zero or not) also collapses to one midpoint and produces the
# same V/triangle when interpolated against neighboring points. Anchoring both
# edges of every wide interval at 1bp (below), not just synthetic gap edges,
# keeps those plateaus flat too.
densifyZeros <- function(subset_data, region) {
    mcols(subset_data)$score[mcols(subset_data)$score < 0] <- 0

    holes <- GenomicRanges::setdiff(region, subset_data)
    if (length(holes) > 0) {
        edges <- c(resize(holes, width = 1, fix = "start"),
                   resize(holes, width = 1, fix = "end"))
        mcols(edges)$score <- 0
        subset_data <- c(subset_data, edges)
    }

    wide <- width(subset_data) > 1
    if (any(wide)) {
        flat <- subset_data[wide]
        subset_data <- c(subset_data[!wide],
                          resize(flat, width = 1, fix = "start"),
                          resize(flat, width = 1, fix = "end"))
    }

    sort(subset_data)
}

subsetNormalizeBigwig <- function(  bigwig_data,
                                    region, samples,
                                    pretty_names,
                                    normalization_color,
                                    normalization_factor = 1,
                                    type = "polygon",
                                    densify = TRUE){
    # Bigwig data is already a list of BigWig files in GRanges format,
    # Identified by name in samples
    # region is the region of interest in which to subset the big_wig data
    ## subsetting by region is necessary to keep RAM usage low!!!
    # pretty_names is a named vector of pretty names for the samples

    # Initialize list to hold DataTracks
    data_tracks <- list()
    colors  <- samples
    samples <- names(samples)

    range_values <- c(0, 0)

    # Loop over samples to calculate normalization of y ranges
    for(i in samples){
        # Subset to the region of interest, filling zero-coverage gaps if requested
        subset_data <- subsetByOverlaps(bigwig_data[[i]], region)
        if (densify) {
            subset_data <- densifyZeros(subset_data, region)
        }

        # Extract scores within the regions
        # The identify a Y value range
        score_values <- mcols(subset_data)$score
        range_values_sample <- range(score_values, na.rm = TRUE)

        # Set minimum Y value
        if(range_values[1] > range_values_sample[1]){
            range_values[1] <- range_values_sample[1]
        }

        # Set maximum Y value
        if(range_values[2] < range_values_sample[2]){
            range_values[2] <- range_values_sample[2]
        }

        # Set Y floor to 0 if the minimum Y value is less than 0
        if(range_values[1] < 0){
            range_values[1] <- 0
        }

    }

    # Normalize the Y-axis limits
    range_values[2] <- range_values[2] * normalization_factor

    # Only mark the top value on the axis, labeled in scientific notation
    top_tick <- range_values[2]
    names(top_tick) <- formatC(top_tick, format = "e", digits = 1)

    # Loop over BigWig files
    for (i in samples) {
        # Subset to the region of interest, filling zero-coverage gaps if requested
        subset_data <- subsetByOverlaps(bigwig_data[[i]], region)
        if (densify) {
            subset_data <- densifyZeros(subset_data, region)
        }

        # Extract scores within the regions
        # The identify a Y value range
        score_values <- mcols(subset_data)$score

        # Determine title color based on normalization color
        title_color <- get_contrasting_text_color(normalization_color)

        # Create a DataTrack for each file
        data_tracks[[i]] <- DataTrack(
            range = subset_data,
            data = score_values,
            type = type,
            name = pretty_names[i],  # Assign a name based on the sample
            col = colors[i],   # Line color
            fill.mountain = c(colors[i], colors[i]),  # Fill color
            fill = colors[i],  # Fill color
            ylim = range_values,
            background.title = normalization_color,
            col.title = title_color,  # Title color
            col.axis = title_color,  # Axis label color, matched to the background
            showAxis = TRUE,
            yTicksAt = top_tick,  # Only the top value, labeled in scientific notation
            fontface.title = 1  # Gviz bolds titles (fontface.title=2) by default; keep it plain

        )
    }

    # Normalize the Y-axis limits witin the same normalization group
    ylims <- range(unlist(lapply(data_tracks, function(track) displayPars(track)$ylim)))
    for (i in samples) {
        displayPars(data_tracks[[i]])$ylim <- ylims
    }

    return(data_tracks)

}

peakHighlighter <- function(    data_tracks,
                                gene,
                                gene_data,
                                peaks,
                                peak_data,
                                region){

    # Separate highlight color from peak identifier
    colors <- peaks
    peaks <- names(peaks)


    # Make Granges object of slected peaks
    peakers <- peak_data[mcols(peak_data)$name %in% peaks]
    mcols(peakers)$color <- colors[mcols(peakers)$name]

    mcols(peakers) <- mcols(peakers) %>% 
        as.data.frame() %>% 
        dplyr::select(name, color)

    # Make Granges object of gene
    chromosome <- as.character(seqnames(gene_data[mcols(gene_data)$gene_id == gene]))
    mcols(gene_data)$name <- mcols(gene_data)$gene_id
    genes_gr <- gene_data[mcols(gene_data)$name == gene]
    mcols(genes_gr)$color <- "#f5f5f5"
    mcols(genes_gr) <- mcols(genes_gr) %>% 
        as.data.frame() %>% 
        dplyr::select(name, color)

    # Combine gene and peaks
    combined_gr <- c(genes_gr, peakers)
    # Subset by region
    combined_gr <- subsetByOverlaps(combined_gr, region)


    # Add start and end positions to the combined data
    mcols(combined_gr)$start <- start(combined_gr)
    mcols(combined_gr)$end <- end(combined_gr)

    # Extact mcols
    mcols_combined <- as.data.frame(mcols(combined_gr))
    mcols_combined

    # Check mcols_combined length
    # If it is greater than 0, proceed
    if (nrow(mcols_combined) > 0) {

        # Extract start and end of region
        start_pos <- start(region)
        end_pos <- end(region)

        # Check start and end position of mcols_combined
        # If start is less than start_pos, set it to start_pos
        if (min(mcols_combined$start) < start_pos) {
            mcols_combined$start[mcols_combined$start < start_pos] <- start_pos
        }
        # If end is greater than end_pos, set it to end_pos
        if (max(mcols_combined$end) > end_pos) {
            mcols_combined$end[mcols_combined$end > end_pos] <- end_pos
        }

        # Compute widths
        mcols_combined$width <- mcols_combined$end - mcols_combined$start +1

        color <- mcols_combined$color

        # # Make highlights
        ht  <- HighlightTrack(  trackList = data_tracks,
                                start = mcols_combined$start, width = mcols_combined$width,
                                chromosome = chromosome,
                                col = color,
                                fill = color)

        return(ht)
    } else {
        return(data_tracks)
    }

}

modelHighlighter <- function(   model_tracks,
                                gene,
                                gene_data, 
                                region){

    # Make Granges object of gene
    chromosome <- as.character(seqnames(gene_data[mcols(gene_data)$gene_id == gene]))
    mcols(gene_data)$name <- mcols(gene_data)$gene_id
    genes_gr <- gene_data[mcols(gene_data)$name == gene]
    mcols(genes_gr)$color <- "#f5f5f5"
    mcols(genes_gr) <- mcols(genes_gr) %>% 
        as.data.frame() %>% 
        dplyr::select(name, color)

    combined_gr <- subsetByOverlaps(genes_gr, region)

    # Add start and end positions to the combined data
    mcols(combined_gr)$start <- start(combined_gr)
    mcols(combined_gr)$end <- end(combined_gr)

    # Extact mcols
    mcols_combined <- as.data.frame(mcols(combined_gr))
    mcols_combined

    # Check mcols_combined length
    # If it is greater than 0, proceed
    if (nrow(mcols_combined) > 0) {

        # Extract start and end of region
        start_pos <- start(region)
        end_pos <- end(region)

        # Check start and end position of mcols_combined
        # If start is less than start_pos, set it to start_pos
        if (min(mcols_combined$start) < start_pos) {
            mcols_combined$start[mcols_combined$start < start_pos] <- start_pos
        }
        # If end is greater than end_pos, set it to end_pos
        if (max(mcols_combined$end) > end_pos) {
            mcols_combined$end[mcols_combined$end > end_pos] <- end_pos
        }

        # Compute widths
        mcols_combined$width <- mcols_combined$end - mcols_combined$start +1

        color <- mcols_combined$color

        # Make highlights
        ht  <- HighlightTrack(  trackList = model_tracks,
                                start = mcols_combined$start, width = mcols_combined$width,
                                chromosome = chromosome,
                                col = color,
                                fill = color)

        return(ht)
    } else {
        return(model_tracks)
    }

}


assemble_data_tracks <- function(   bigwig_data,
                                    normalization_groups,
                                    group_title_panel_colors = NULL,
                                    pretty_names,
                                    region,
                                    normalization_factor,
                                    type = NULL,
                                    densify = NULL
){
    # Get normalization group names
    groups <- names(normalization_groups)

    # Initialize data_tracks
    data_tracks <- list()

    # Compute subset and normalize along normalization groups
    for(i in seq_along(normalization_groups)){

        # Define group name
        group_name <- groups[i]

        if (is.null(group_title_panel_colors)) {
            # Recolor title panel
            norm_color <- "lightgray"

        } else {
            # Recolor title panel
            norm_color <- group_title_panel_colors[group_name]
        }

        if (is.null(type)) {
            # Default plot shape
            group_type <- "polygon"

        } else {
            # Plot shape for this group
            group_type <- type[group_name]
        }

        if (is.null(densify)) {
            # Default: fill zero-coverage gaps
            group_densify <- TRUE

        } else {
            # Whether to fill zero-coverage gaps for this group
            group_densify <- densify[group_name]
        }

        data_tracks[[group_name]] <- subsetNormalizeBigwig( bigwig_data = bigwig_data,
                                                            region = region,
                                                            samples = normalization_groups[[group_name]],
                                                            pretty_names = pretty_names,
                                                            normalization_color = norm_color,
                                                            normalization_factor = normalization_factor[group_name],
                                                            type = group_type,
                                                            densify = group_densify)

    }

    data_tracks <- unlist(data_tracks)

    return(data_tracks)
}

get_contrasting_text_color <- function(bg_color) {
    rgb <- col2rgb(bg_color) / 255  # Normalize
    linear_rgb <- ifelse(rgb <= 0.03928, rgb / 12.92, ((rgb + 0.055) / 1.055)^2.4)
    luminance <- 0.2126 * linear_rgb[1] + 0.7152 * linear_rgb[2] + 0.0722 * linear_rgb[3]
    if (luminance > 0.5) "black" else "white"
}

# Function to identify relevant daps
identifyDaps <- function(   daps,
                            gene,
                            gene_data,
                            highlight_colors,
                            naive_relative = FALSE
                            ) {

    # Determine the chromosome of the gene, and it's start and end positions
    chromosome <- as.character(seqnames(gene_data[mcols(gene_data)$gene_id == gene]))
    print(chromosome)

    # Make naive relative filters:
    if(naive_relative){
        # First isolate shCD19 vs naive comparisons
        control_comparisons <- daps %>% 
            filter(comparison %in% c("GFPpos_shCD19", "GFPneg_shCD19", "CX3CR1pos_shCD19"))

        # Then isolate regions that open relative to naive (union of all)
        # Convert to vector of names
        open_vs_naive <- control_comparisons %>%
            filter(log2FC > 0) %>%
            dplyr::select(name) %>% unique()
        open_vs_naive <- open_vs_naive$name

        # Then isolate regions that open relative to naive (union of all)
        # Convert to vector of names
        close_vs_naive <- control_comparisons %>%
            filter(log2FC < 0) %>%
            dplyr::select(name) %>% unique()
        close_vs_naive <- close_vs_naive$name
    } else {
        open_vs_naive <- close_vs_naive <- unique(daps$name)  # If not naive relative, just use all names
    }

    # Make a gene filter
    filter <- paste0("^", gene, "$|;", gene, ";|^", gene, ";|;", gene, "$")

    select <- daps %>%
        filter( grepl(filter, flank_geneSymbols),                               # filter for gene
                grepl("shRunx3_vs_shCD19", comparison)) %>%                     # filter for comparisons
        dplyr::select(comparison, name, log2FC, location) %>%                          # location is of the form chr:start-end, split into three columns
        separate(location, into = c("chr", "start", "end"), sep = "[:-]") %>%   # append "chr" to prefix of chr
        mutate(chr = paste0("chr", chr)) %>%                                    # convert start and end to numeric
        mutate(start = as.numeric(start), end = as.numeric(end)) %>%            # remove comparison with CX3CR1pos or GFPneg
        filter(!grepl("CX3CR1pos|GFPneg", comparison)) %>%                      # select only daps in the same chromosome as the gene
        filter(chr == chromosome)

    # In select, make a new column called type
    # If name is in open_vs_naive, comparison is "GFPpos_shRunx3_vs_shCD19" and log2FC < 0, then type is "opened by Runx3 in GFPpos"
    select <- select %>%
        mutate(type = case_when(
            name %in% open_vs_naive & comparison == "GFPpos_shRunx3_vs_shCD19" & log2FC < 0 ~ "Runx3opened_GFPpos",
            name %in% open_vs_naive & comparison == "shRunx3_vs_shCD19" & log2FC < 0 ~ "Runx3opened_all",
            name %in% close_vs_naive & comparison == "GFPpos_shRunx3_vs_shCD19" & log2FC > 0 ~ "Runx3closed_GFPpos",
            name %in% close_vs_naive & comparison == "shRunx3_vs_shCD19" & log2FC > 0 ~ "Runx3closed_all",
            TRUE ~ NA_character_  # Default case if none of the conditions are met
        )) %>%
        filter(!is.na(type))  # Remove rows where type is NA

    # Order select by type, in the following order:
    # "Runx3 opened, GFPpos", "Runx3 opened, all", "Runx3 closed, GFPpos", "Runx3 closed, all"
    type_order <- c("Runx3opened_GFPpos", "Runx3opened_all", "Runx3closed_GFPpos", "Runx3closed_all")
    select$type <- factor(select$type, levels = type_order)
    select <- select %>%
        arrange(type, chr, start) %>%  # Make a column to see if name is duplicated
        mutate(duplicated_name = duplicated(name)) %>% # Select not duplicated names
        filter(!duplicated_name) %>%  # Remove duplicated names
        dplyr::select(-duplicated_name)  # Remove the duplicated_name column

    # Make a new column called highlight, this is dependent on the type and named vector highlight_colors
    select <- select %>%
        mutate(highlight = case_when(
            type == "Runx3opened_GFPpos" ~ highlight_colors["Runx3opened_GFPpos"],
            type == "Runx3opened_all" ~ highlight_colors["Runx3opened_all"],
            type == "Runx3closed_GFPpos" ~ highlight_colors["Runx3closed_GFPpos"],
            type == "Runx3closed_all" ~ highlight_colors["Runx3closed_all"],
            TRUE ~ NA_character_  # Default case if none of the conditions are met
        ))

    return(select)
}

# A function to determine plot bounds based on gene location, and optionally DAPs
# When identifiedDaps is NULL (or empty), bounds fall back to just the gene's own
# coordinates padded by `margin` on each side.
determinePlotBoundsAuto <- function(    gene,
                                        gene_data,
                                        identifiedDaps = NULL,
                                        margin = 5000){
    # Subset GRanges by gene_id
    genes_df <- gene_data[mcols(gene_data)$gene_id %in% gene, ]

    # Convert to data frame and select only the desired columns
    genes_df <- as.data.frame(genes_df)[, c("gene_id", "seqnames", "start", "end")] %>%
        rename(chr = seqnames, name = gene_id)

    # If genes_df is empty or more than one gene is selected, die
    if (nrow(genes_df) == 0 || nrow(genes_df) > 1) {
        stop("Error: No gene found or more than one gene selected.")
    }

    # If DAPs were identified for this gene, fold their bounds in too
    if (!is.null(identifiedDaps) && nrow(identifiedDaps) > 0) {
        # Select only name, chr, start, and end columns of identifiedDaps
        identifiedDaps <- identifiedDaps %>%
            dplyr::select(name, chr, start, end)

        # Row bind the two data frames
        genes_df <- rbind(genes_df, identifiedDaps)
    }

    # Isolate chomosome
    chromosome <- unique(genes_df$chr)
    # Isolate start and end positions, padded by margin on each side
    start_pos <- min(genes_df$start) - margin
    end_pos <- max(genes_df$end) + margin

    # Make Granges object only of largest region
    region <- GRanges(chromosome, IRanges(start_pos, end_pos))

    return(region)
}

plot_gene_auto <- function( gene,
                            genome,
                            bw_data,
                            norm_factor,
                            normalization_groups,
                            group_panel_colors,
                            pretty_names,
                            gene_data,
                            model_tracks,
                            daps = NULL,
                            highlight_colors = NULL,
                            peak_data = NULL,
                            margin = 5000,
                            naive_relative = FALSE,
                            type = NULL,
                            densify = NULL,
                            includeIdeogram = TRUE
)
{

    identifiedDaps <- NULL

    # Only identify DAPs if DAP/peak data was actually supplied
    if (!is.null(daps) && !is.null(peak_data)) {
        identifiedDaps <- identifyDaps( daps = daps,
                                        gene = gene,
                                        gene_data = gene_data,
                                        highlight_colors = highlight_colors,
                                        naive_relative = naive_relative)
    }

    # Automatically determine the region of interest.
    # Falls back to just the gene's own coordinates +/- margin when no DAPs are supplied.
    region <- determinePlotBoundsAuto(  gene = gene,
                                        gene_data = gene_data,
                                        identifiedDaps = identifiedDaps,
                                        margin = margin)

    # Determine regions
    chromosome <- as.character(seqnames(region))
    start_pos <- start(region)
    end_pos <- end(region)

    #### RUN PLOTTING ####
    # See plotTracks_scripts.R for more details
    # This returns a list of tracks to plot
    data_tracks <- assemble_data_tracks(
        bigwig_data = bw_data,                          # List containing bigwig GRanges
        normalization_groups = normalization_groups,    # BigWig Tracks to plor grouped by normalization groups
        group_title_panel_colors = group_panel_colors,  # Dictionary of colors for the title panes of each normalization group
        pretty_names = pretty_names,                    # Pretty names of the samples
        region = region,                                # Region to plot
        normalization_factor = norm_factor,             # Normalization factor for each normalization group
        type = type,                                    # DataTrack plot type (e.g. "polygon", "l", "histogram")
        densify = densify                                # Whether to fill zero-coverage gaps, per normalization group
    )

    # Un-bold gene labels in the model tracks (Gviz's GeneRegionTrack defaults
    # fontface.group to 2/bold, and preprocessTracks.R/plotTracks.R may set it
    # even bolder) -- override here so no track text renders bold. Must run
    # before modelHighlighter() below, since that call can replace model_tracks
    # with a single HighlightTrack (not a subsettable list) when DAPs overlap
    # the plotted region.
    for (i in seq_along(model_tracks)) {
        if (is(model_tracks[[i]], "GeneRegionTrack")) {
            displayPars(model_tracks[[i]])$fontface.group <- 3  # italic, not bold
        }
    }

    # Only highlight if DAPs were actually identified
    if (!is.null(identifiedDaps)) {

        # Peaks of interest
        # Define colors and peaks from identifiedDaps
        poi <- identifiedDaps$highlight
        names(poi) <- identifiedDaps$name

        # Highlight Peaks in Data
        data_tracks <- peakHighlighter( data_tracks = data_tracks,     # List of data tracks to highlight
                                    gene = gene,                            # Gene to highlight
                                    gene_data = gene_data,                  # Gene data
                                    peaks = poi,                            # Peaks to highlight
                                    peak_data = peak_data,
                                    region = region
        )

        # Highlight Gene in Model
        model_tracks <- modelHighlighter(   model_tracks = model_tracks,
                                        gene = gene,
                                        gene_data = gene_data,
                                        region = region
        )
    }

    # Create Ideogram Track
    # Skipped for non-UCSC chromosomes (e.g. viral contigs) where an ideogram
    # cannot be computed -- IdeogramTrack() would otherwise fail or hang trying
    # to fetch cytoband data for a chromosome UCSC doesn't know about.
    if (includeIdeogram) {
        ideo_track <- IdeogramTrack(
            genome = genome,
            showBandId = TRUE,
            chromosome = chromosome,
            fontcolor = "black",
            fontface.title = 1  # Gviz bolds titles (fontface.title=2) by default; keep it plain
        )
    } else {
        ideo_track <- NULL
    }

    # Combine lists of tracks
    # This is the final list of tracks to plot
    plot <- c(  ideo_track,
                model_tracks,
                data_tracks
    )


    return(list(tracks = plot, chromosome = chromosome, start_pos = start_pos, end_pos = end_pos))
}


plot_gene_manual <- function(   gene,
                                genome,
                                chromosome,
                                start_pos,
                                end_pos,
                                bw_data,
                                norm_factor,
                                normalization_groups,
                                group_panel_colors,
                                pretty_names,
                                gene_data,
                                model_tracks,
                                daps = NULL,
                                highlight_colors = NULL,
                                peak_data = NULL,
                                naive_relative = FALSE,
                                type = NULL,
                                densify = NULL,
                                includeIdeogram = TRUE
)
{

    # Region is taken directly from the manually supplied coordinates
    region <- GRanges(chromosome, IRanges(start_pos, end_pos))

    #### RUN PLOTTING ####
    # See plotTracks_scripts.R for more details
    # This returns a list of tracks to plot
    data_tracks <- assemble_data_tracks(
        bigwig_data = bw_data,                          # List containing bigwig GRanges
        normalization_groups = normalization_groups,    # BigWig Tracks to plor grouped by normalization groups
        group_title_panel_colors = group_panel_colors,  # Dictionary of colors for the title panes of each normalization group
        pretty_names = pretty_names,                    # Pretty names of the samples
        region = region,                                # Region to plot
        normalization_factor = norm_factor,             # Normalization factor for each normalization group
        type = type,                                    # DataTrack plot type (e.g. "polygon", "l", "histogram")
        densify = densify                                # Whether to fill zero-coverage gaps, per normalization group
    )

    # Un-bold gene labels in the model tracks (Gviz's GeneRegionTrack defaults
    # fontface.group to 2/bold, and preprocessTracks.R/plotTracks.R may set it
    # even bolder) -- override here so no track text renders bold. Must run
    # before modelHighlighter() below, since that call can replace model_tracks
    # with a single HighlightTrack (not a subsettable list) when DAPs overlap
    # the plotted region.
    for (i in seq_along(model_tracks)) {
        if (is(model_tracks[[i]], "GeneRegionTrack")) {
            displayPars(model_tracks[[i]])$fontface.group <- 3  # italic, not bold
        }
    }

    # Only highlight if DAP/peak information was actually supplied
    if (!is.null(daps) && !is.null(peak_data)) {

        # Determine cis-DAPs for the gene
        identifiedDaps <- identifyDaps( daps = daps,
                                        gene = gene,
                                        gene_data = gene_data,
                                        highlight_colors = highlight_colors,
                                        naive_relative = naive_relative)

        # Peaks of interest
        # Define colors and peaks from identifiedDaps
        poi <- identifiedDaps$highlight
        names(poi) <- identifiedDaps$name

        # Highlight Peaks in Data
        data_tracks <- peakHighlighter( data_tracks = data_tracks,     # List of data tracks to highlight
                                    gene = gene,                            # Gene to highlight
                                    gene_data = gene_data,                  # Gene data
                                    peaks = poi,                            # Peaks to highlight
                                    peak_data = peak_data,
                                    region = region
        )

        # Highlight Gene in Model
        model_tracks <- modelHighlighter(   model_tracks = model_tracks,
                                        gene = gene,
                                        gene_data = gene_data,
                                        region = region
        )
    }

    # Create Ideogram Track
    # Skipped for non-UCSC chromosomes (e.g. viral contigs) where an ideogram
    # cannot be computed -- IdeogramTrack() would otherwise fail or hang trying
    # to fetch cytoband data for a chromosome UCSC doesn't know about.
    if (includeIdeogram) {
        ideo_track <- IdeogramTrack(
            genome = genome,
            showBandId = TRUE,
            chromosome = chromosome,
            fontcolor = "black",
            fontface.title = 1  # Gviz bolds titles (fontface.title=2) by default; keep it plain
        )
    } else {
        ideo_track <- NULL
    }

    # Combine lists of tracks
    # This is the final list of tracks to plot
    plot <- c(  ideo_track,
                model_tracks,
                data_tracks
    )


    return(list(tracks = plot, chromosome = chromosome, start_pos = start_pos, end_pos = end_pos))
}


#' Merge a viral genome GTF into the host gene annotation
#'
#' The viral GTFs (tracks/source_data/viralgenomes/<rv>.gtf) carry only `exon`
#' lines and a minimal attribute set, while the host GTF is a full ENSEMBL
#' annotation. This helper reconciles the two so that c() succeeds and so that
#' makeTxDbFromGRanges() treats the viral construct like any other contig:
#'   1. prefixes the viral seqlevel with "chr" (pMIA_Runx3 -> chrpMIA_Runx3)
#'   2. fills the attributes the viral GTF omits (transcript_biotype, exon_number, ...)
#'   3. synthesises the transcript- and gene-level rows the viral GTF lacks
#'   4. pads both mcols sets to a common set of columns and flattens factors
#'   5. stamps a single genome tag on the combined object
#'
#' @param host_data   GRanges from import()-ing the host GTF
#' @param viral_gtf   path to the viral GTF
#' @param genome_tag  genome label applied to ALL seqlevels of the result
#' @param add_chr_prefix  prepend "chr" to viral seqlevels that lack it
#' @return a single GRanges holding host + viral annotation
combine_host_viral_gtf <- function( host_data,
                                    viral_gtf,
                                    genome_tag = "mm10",
                                    add_chr_prefix = TRUE) {

    viral_data <- import(viral_gtf, format = "gtf")

    ## ---- 1. chr prefix on the viral contig -------------------------------
    # The ideogram/UCSC-style tracks and the host GTF are chr-prefixed, so the
    # viral contig has to follow suit or it will not match anything downstream.
    if (add_chr_prefix) {
        lv <- seqlevels(viral_data)
        needs <- !grepl("^chr", lv)
        lv[needs] <- paste0("chr", lv[needs])
        seqlevels(viral_data) <- lv
    }

    ## ---- factor columns ---------------------------------------------------
    # rtracklayer returns `source` and `type` as factors. Binding DataFrames
    # whose factor columns have disjoint level sets is the usual source of
    # silently-introduced NAs, so flatten them to character on both objects.
    flatten_factors <- function(gr) {
        for (j in which(vapply(mcols(gr), is.factor, logical(1)))) {
            mcols(gr)[[j]] <- as.character(mcols(gr)[[j]])
        }
        gr
    }
    host_data  <- flatten_factors(host_data)
    viral_data <- flatten_factors(viral_data)

    ## ---- 2. fill the attributes the viral GTF does not carry --------------
    m <- mcols(viral_data)
    if (is.null(m$gene_name))       m$gene_name       <- m$gene_id
    if (is.null(m$gene_biotype))    m$gene_biotype    <- NA_character_
    if (is.null(m$transcript_name)) m$transcript_name <- m$transcript_id
    if (is.null(m$exon_number))     m$exon_number     <- NA_character_
    if (is.null(m$transcript_biotype)) m$transcript_biotype <- NA_character_

    m$gene_biotype[is.na(m$gene_biotype)] <- "protein_coding"
    # The coding-only subset filters on transcript_biotype, which these GTFs omit;
    # inherit it from gene_biotype so viral features survive that filter.
    na_tb <- is.na(m$transcript_biotype)
    m$transcript_biotype[na_tb] <- m$gene_biotype[na_tb]
    mcols(viral_data) <- m

    ## ---- exon ranks -------------------------------------------------------
    # makeTxDbFromGRanges() takes exon rank from exon_number when that column
    # exists, and an all-NA column will trip it up. Rank exons along the
    # direction of transcription within each transcript.
    is_exon <- mcols(viral_data)$type == "exon"
    if (any(is_exon)) {
        ex   <- viral_data[is_exon]
        key  <- ifelse(as.character(strand(ex)) == "-", -start(ex), start(ex))
        ord  <- order(mcols(ex)$transcript_id, key)
        rk   <- integer(length(ex))
        rk[ord] <- ave(seq_along(ord), mcols(ex)$transcript_id[ord], FUN = seq_along)
        mcols(viral_data)$exon_number[is_exon] <- as.character(rk)
    }

    ## ---- 3. synthesise transcript- and gene-level rows --------------------
    # Collapse the exons of each transcript / gene into a single parent range and
    # inherit that feature's attributes from its first exon.
    make_parents <- function(exons, id_col, type_label) {
        ids     <- mcols(exons)[[id_col]]
        parents <- unlist(range(split(exons, ids)), use.names = TRUE)
        if (length(parents) != length(unique(ids))) {
            warning("Some ", id_col, " values span multiple contigs or strands; ",
                    "their parent features were split.")
        }
        mcols(parents) <- mcols(exons)[match(names(parents), ids), , drop = FALSE]
        mcols(parents)$type        <- type_label
        mcols(parents)$exon_number <- NA_character_
        if (type_label == "gene") {
            # Mirrors ENSEMBL: gene rows carry no transcript-level attributes,
            # which is also why they drop out of the coding-only subset.
            mcols(parents)$transcript_id      <- NA_character_
            mcols(parents)$transcript_name    <- NA_character_
            mcols(parents)$transcript_biotype <- NA_character_
        }
        names(parents) <- NULL
        parents
    }

    if (any(is_exon)) {
        ex <- viral_data[is_exon]
        if (!any(mcols(viral_data)$type %in% c("transcript", "mRNA"))) {
            viral_data <- c(viral_data, make_parents(ex, "transcript_id", "transcript"))
        }
        if (!any(mcols(viral_data)$type == "gene")) {
            viral_data <- c(viral_data, make_parents(ex, "gene_id", "gene"))
        }
    }

    ## ---- 4. reconcile mcols before binding --------------------------------
    # import() only creates columns for attributes a file actually contains, and
    # c() on GRanges requires identical mcols. Pad both to the union, host order first.
    common_cols <- union(colnames(mcols(host_data)), colnames(mcols(viral_data)))
    pad_mcols <- function(gr, cols) {
        for (mi in setdiff(cols, colnames(mcols(gr)))) mcols(gr)[[mi]] <- NA_character_
        mcols(gr) <- mcols(gr)[, cols, drop = FALSE]
        gr
    }
    host_data  <- pad_mcols(host_data,  common_cols)
    viral_data <- pad_mcols(viral_data, common_cols)

    ## ---- 5. one genome tag across the board -------------------------------
    # c() refuses to bind GRanges whose genome() tags disagree, and everything
    # downstream (Gviz, the API) is told the genome is mm10 - so the viral contig
    # is carried as an extra "mm10" sequence.
    genome(host_data)  <- genome_tag
    genome(viral_data) <- genome_tag

    sort(c(host_data, viral_data), ignore.strand = TRUE)
}

