# plotSEA.R
# R scripts for plotting GSEA/PSEA results

plotSEAbubble <- function(  sea_results,
                            comparison,
                            enrichment_type,
                            cluster_labels = NULL,
                            color_scale = viridis::cividis(3)) {
    # Define color scale
    #color_scale <- viridis::cividis(3)  # Generate a color scale with 3 colors

    # Ensure limits are symmetric around 0
    NES_range <- range(sea_results$NES, na.rm = TRUE)
    max_abs_NES <- max(abs(NES_range))
    limits <- c(-max_abs_NES, max_abs_NES)

    # Create bubble plot
    p <- ggplot(sea_results, aes(x = cluster, y = ID)) +
        geom_point(aes(color = NES, size = log_qvalue), stroke = 1, shape = 21, fill = NA, colour = "black") +
        geom_point(aes(color = NES, size = log_qvalue)) +
        scale_color_gradientn(
            colors = color_scale,
            limits = limits,
            guide = guide_colorbar(title = "NES")
        ) +
        scale_size_continuous(
            breaks = seq(0, max(sea_results$log_qvalue), by = 0.5),
            labels = seq(0, max(sea_results$log_qvalue), by = 0.5),
            guide = guide_legend(title = "-log10(Q-value)")
        ) +
        labs(
            title = paste(enrichment_type, "Set Enrichment", comparison),  # Set the plot title
            x = "Cluster",  # Label for the X-axis
            y = "ID"  # Label for the Y-axis
        ) +
        theme_minimal() +  # Set the plot theme to minimal
        theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate and align X-axis labels

    # Apply user-provided cluster labels if given
    if (!is.null(cluster_labels)) {
        p <- p + scale_x_discrete(labels = cluster_labels)
    }

    return(p)
}

plotSEAbubble2 <- function(sea_results, comparison, enrichment_type, cluster_labels) {
    # Define color scale
    color_scale <- viridis::cividis(3)  # Generate a color scale with 3 colors

    # Ensure limits are symmetric around 0
    NES_range <- range(sea_results$NES, na.rm = TRUE)
    max_abs_NES <- max(abs(NES_range))
    limits <- c(-max_abs_NES, max_abs_NES)

    # Create bubble plot
    p <- ggplot(sea_results, aes(x = cluster, y = ID)) +
        geom_point(aes(color = NES, size = log_qvalue)) +
        scale_color_gradientn(
            colors = color_scale,
            limits = limits,
            guide = guide_colorbar(title = "NES")
        ) +
        scale_size_continuous(
            breaks = seq(0, max(sea_results$log_qvalue), by = 0.5),
            labels = seq(0, max(sea_results$log_qvalue), by = 0.5),
            guide = guide_legend(title = "-log10(Q-value)")
        ) +
        scale_x_discrete(labels = cluster_labels) +  # Apply user-provided cluster labels
        labs(
            title = paste(enrichment_type, "Set Enrichment", comparison),  # Set the plot title
            x = NULL,  # Label for the X-axis
            y = NULL  # Label for the Y-axis
        ) +
        theme_minimal() +  # Set the plot theme to minimal
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate and align X-axis labels
            legend.position = "bottom",  # Move legends to bottom
            legend.box = "horizontal"   # Arrange legends side by side
        )

    return(p)
}

plotGREATbubble <- function(great_results, title = "GREAT Enrichment", fill_scale = rev(RColorBrewer::brewer.pal(11, "Spectral"))) {
    # great_results columns expected: group (x), id (y), fold_enrichment (fill), p_adjust (size)
    # p_adjust can round to exactly 0 at floating-point precision; floor it at the smallest
    # nonzero p_adjust in the data so -log10() doesn't produce an outlier that swamps the
    # size scale (flooring at .Machine$double.xmin instead makes that point ~20x any other)
    nonzero_padj <- great_results$p_adjust[great_results$p_adjust > 0]
    padj_floor <- if (length(nonzero_padj) > 0) min(nonzero_padj) else .Machine$double.xmin
    great_results$neg_log10_padj <- -log10(pmax(great_results$p_adjust, padj_floor))

    # Create bubble plot
    p <- ggplot(great_results, aes(x = group, y = id)) +
        geom_point(aes(size = neg_log10_padj), stroke = 1, shape = 21, fill = NA, colour = "black") +
        geom_point(aes(fill = fold_enrichment, size = neg_log10_padj), shape = 21, colour = "black") +
        scale_fill_gradientn(
            colors = fill_scale,
            guide = guide_colorbar(title = "Fold Enrichment")
        ) +
        scale_size_continuous(
            guide = guide_legend(title = "-log10(p adjust)")
        ) +
        labs(
            title = title,  # Set the plot title
            x = NULL,  # Label for the X-axis
            y = NULL  # Label for the Y-axis
        ) +
        theme_minimal() +  # Set the plot theme to minimal
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate and align X-axis labels
            legend.position = "bottom",  # Move legends to bottom
            legend.box = "horizontal"   # Arrange legends side by side
        )

    return(p)
}

