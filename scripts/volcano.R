library(ggplot2)
library(ggrepel)
library(ggsci)
library(ggpubr)
library(ggprism)
library(ggpp)

# processDegs is a function that removes non-significant genes and randomly selects genes to remove
# This is because SVG files for volcano plots are too complicated to be edited in Illustrator
# This function will help you to remove non-significant genes and randomly select genes to remove
# It will also mark genes to keep: genes in the 'genes' vector and top 10% p-value ranked genes

# degs is the data frame with the DEGs
# genes is a vector with the genes you want to keep (for naming purposes)

processDegs <- function(degs, genes, seed = 42) {
    # Make maximum P value
    degs$pvals_adj[degs$pvals_adj == 0] <- .Machine$double.xmin
    # Marking non-significant and significant genes
    degs$Significant <- degs$pvals_adj < 0.05
    
    # Mark genes to keep: genes in the 'genes' vector and top 10% p-value ranked genes
    top_10_percentile <- quantile(degs$pvals_adj[degs$Significant], 0.1)
    degs$Keep <- with(degs, names %in% genes | (Significant & pvals_adj <= top_10_percentile))
    
    # For non-significant points, randomly select 95% to remove
    set.seed(seed) # For reproducibility
    degs$Remove <- with(degs, ifelse(!Significant & !Keep, runif(nrow(degs)) < 95/100, FALSE))
    
    # Calculate IQR for logFC among significant genes
    iqr_logFC <- IQR(degs$logfoldchanges[degs$Significant & !degs$Remove])
    
    # Identify points within the IQR for logFC
    lower_bound <- quantile(degs$logfoldchanges[degs$Significant & !degs$Remove], 0.25)
    upper_bound <- quantile(degs$logfoldchanges[degs$Significant & !degs$Remove], 0.75)
    degs$WithinIQR <- with(degs, logfoldchanges >= lower_bound & logfoldchanges <= upper_bound & Significant & !Remove)
    
    # Randomly remove 3/4 of genes within the IQR for logFC
    degs$Remove <- with(degs, ifelse(WithinIQR & !Keep, runif(nrow(degs)) < 3/4, Remove))

    # Make a new column, label
    # If names is in genes AND Significant == TRUE then label = names,
    # Otherwise label = NA
    degs <- degs %>%
        mutate(label = ifelse(Significant & names %in% genes, names, NA))

    return(degs)
}

processDegsStrict <- function(degs, genes, seed = 42) {
    # Make maximum P value
    degs$pvals_adj[degs$pvals_adj == 0] <- .Machine$double.xmin
    # Marking non-significant and significant genes
    degs$Significant <- degs$pvals_adj < 0.05
    
    # Mark genes to keep: genes in the 'genes' vector and top 10% p-value ranked genes
    top_10_percentile <- quantile(degs$pvals_adj[degs$Significant], 0.1)
    degs$Keep <- with(degs, names %in% genes | (Significant & pvals_adj <= top_10_percentile))
    
    # For non-significant points, randomly select 99% to remove
    set.seed(seed) # For reproducibility
    degs$Remove <- with(degs, ifelse(!Significant & !Keep, runif(nrow(degs)) < 99/100, FALSE))
    
    # Calculate IQR for logFC among significant genes
    iqr_logFC <- IQR(degs$logfoldchanges[degs$Significant & !degs$Remove])
    
    # Identify points within the IQR for logFC
    lower_bound <- quantile(degs$logfoldchanges[degs$Significant & !degs$Remove], 0.25)
    upper_bound <- quantile(degs$logfoldchanges[degs$Significant & !degs$Remove], 0.75)
    degs$WithinIQR <- with(degs, logfoldchanges >= lower_bound & logfoldchanges <= upper_bound & Significant & !Remove)
    
    # Randomly remove 9/10 of genes within the IQR for logFC
    degs$Remove <- with(degs, ifelse(WithinIQR & !Keep, runif(nrow(degs)) < 9/10, Remove))

    # Make a new column, label
    # If names is in genes AND Significant == TRUE then label = names,
    # Otherwise label = NA
    degs <- degs %>%
        mutate(label = ifelse(Significant & names %in% genes, names, NA))

    return(degs)
}


# Assuming your data frame is named degs
# Function to plot volcanoes
volcanoPlot <- function(degs,
                        genes,
                        theme,
                        pval = 0.05,
                        logfc = 0.2,
                        up_color = "#3C5488",
                        down_color = "#DC0000",
                        xlim = c(-8, 8),
                        ylim = c(-2, 400),
                        label_x_clearance = c(-4, 4),
                        overlaps = 30,
                        title = "Volcano Plot of p-values vs log fold change",
                        seed = 42) {
    degs$Significant <- with(degs, pvals_adj < pval & abs(logfoldchanges) > logfc)
    degs$Color <- with(degs, ifelse(logfoldchanges > 0 & Significant, up_color, ifelse(logfoldchanges < 0 & Significant, down_color, "grey50")))

    # Remove observations if Remove == TRUE
    degs <- degs[!degs$Remove, ]

    # Randomize obeservation order of the data frame
    set.seed(seed) # For reproducibility
    degs <- degs[sample(nrow(degs)), ]

    # Creating the plot
    volcano_plot <- ggplot(degs, aes(x = logfoldchanges, y = -log10(pvals_adj), color = Color, label = label)) +
        geom_point(alpha = 1) +
        scale_color_identity() +
        # Downregulated labels
        geom_text_repel(data = degs %>% filter(!is.na(label), logfoldchanges < 0),
                        aes(label = label, color = Color),
                        fontface = "italic",
                        position = position_nudge_keep(x = -0.5, y = 50),
                        box.padding = 0.5,
                        point.padding = 0.5,
                        segment.color = 'black', # Color of the line
                        segment.size = 1 /.pt, # Thickness of the line
                        min.segment.length = 0,
                        segment.alpha = 1,
                        max.overlaps = overlaps, # Maximum number of overlaps
                        size = 12 /.pt,
                        force = 10,
                        force_pull = 0.2,
                        max.iter = 50000,
                        max.time = 60,
                        seed = seed,
                        xlim = c(min(xlim), min(label_x_clearance)),
                        ylim = c(min(ylim), max(ylim))
        ) +
        # Upregulated labels
        geom_text_repel(data = degs %>% filter(!is.na(label), logfoldchanges > 0),
                        aes(label = label, color = Color),
                        fontface = "italic",
                        position = position_nudge_keep(x = 0.5, y = 50),
                        box.padding = 0.5,
                        point.padding = 0.5,
                        segment.color = 'black', # Color of the line
                        segment.size = 1 / .pt, # Thickness of the line
                        min.segment.length = 0,
                        segment.alpha = 1,
                        max.overlaps = overlaps, # Maximum number of overlaps
                        size = 12 /.pt,
                        force = 10,
                        force_pull = 0.2,
                        max.iter = 50000,
                        max.time = 60,
                        seed = seed,
                        xlim = c(max(label_x_clearance), max(xlim)),
                        ylim = c(min(ylim), max(ylim))
        ) +
        theme +
        labs(title = title,
            x = "Log Fold Change",
            y = "-log(adjusted p-value)") +
        geom_vline(xintercept = c(-logfc, logfc), linetype = "dashed") +
        geom_hline(yintercept = -log10(pval), linetype = "dashed") +
        coord_cartesian(xlim = xlim, ylim = ylim)

    return(volcano_plot)
}

# Use ggprism, base 12, not bold, and remove legend
theme <- theme_bw(      base_size = 14,
                        base_family = "sans"
                        ) +
    theme(  axis.text = element_text(size = 12, color = "black"),
            axis.title = element_text(size = 12),
            strip.text = element_text(size = 12),
            legend.position = "none",
            plot.title = element_text(size = 12),  # Set main title size
            panel.grid = element_blank())  # Remove gridlines


# Assuming your data frame is named degs
# Function to plot volcanoes
volcanoPlotBak <- function(degs,
                        genes,
                        theme,
                        pval = 0.05,
                        logfc = 0.2,
                        up_color = "#3C5488",
                        down_color = "#DC0000",
                        xlim = c(-8, 8),
                        ylim = c(-2, 400),
                        overlaps = 30,
                        title = "Volcano Plot of p-values vs log fold change",
                        seed = 42) {
    degs$Significant <- with(degs, pvals_adj < pval & abs(logfoldchanges) > logfc)
    degs$Color <- with(degs, ifelse(logfoldchanges > 0 & Significant, up_color, ifelse(logfoldchanges < 0 & Significant, down_color, "grey50")))

    # Randomize obeservation order of the data frame
    set.seed(seed) # For reproducibility
    degs <- degs[sample(nrow(degs)), ]

    # Creating the plot
    volcano_plot <- ggplot(degs, aes(x = logfoldchanges, y = -log10(pvals_adj), color = Color, label = label)) +
        geom_point(alpha = 1) +
        scale_color_identity() +
        geom_text_repel(box.padding = 0.5,
                        point.padding = 0.5,
                        segment.color = 'black', # Color of the line
                        segment.size = 1, # Thickness of the line
                        min.segment.length = 0,
                        segment.alpha = 1,
                        max.overlaps = overlaps, # Maximum number of overlaps
                        size = 6,
                        force = 50,
                        force_pull = 0.2,
                        max.iter = 50000,
                        max.time = 60,
                        seed = seed) + # Transparency of the line
        theme +
        labs(title = title,
            x = "Log Fold Change",
            y = "-log(adj. p-value)") +
        geom_vline(xintercept = c(-logfc, logfc), linetype = "dashed") +
        geom_hline(yintercept = -log10(pval), linetype = "dashed") +
        
        coord_cartesian(xlim = xlim, ylim = ylim)

    return(volcano_plot)
}

# Use ggprism, base 12, not bold, and remove legend
theme <- theme_bw(      base_size = 14,
                        base_family = "sans"
                        ) +
    theme(  axis.text = element_text(size = 12, color = "black"),
            axis.title = element_text(size = 12),
            strip.text = element_text(size = 12),
            legend.position = "none",
            plot.title = element_text(size = 12),  # Set main title size
            panel.grid = element_blank())  # Remove gridlines