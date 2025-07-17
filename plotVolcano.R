library(ggplot2)
library(ggrepel)
library(gridExtra)
library(RColorBrewer)
library(ggrastr)

plotVolcano <- function(
  data,                      # Input data with all genes
  gene.col = "Gene",         # Column name for gene names
  logFC.col = "logFC",       # Column name for log fold changes
  pval.col = "P.Value",      # Column name for p-values
  selected_genes,            # Data frame with selected genes and their pathways
  plot_mode = "advanced",    # Plotting mode ("classic" or "advanced")
  logFCcut = 1.5,            # Log2 fold-change threshold
  logFCcut2 = 2.5,           # Secondary log2 fold-change threshold (advanced mode)
  logFCcut3 = 5,             # Tertiary log2 fold-change threshold (advanced mode)
  pvalCut = 0.05,            # P-value threshold
  pvalCut2 = 0.0001,         # Secondary p-value threshold (advanced mode)
  pvalCut3 = 0.00001,        # Tertiary p-value threshold (advanced mode)
  pathway_colors = c("blueviolet","#223D6C","darkgreen","chocolate4","#D20A13","#088247","#58CDD9","#7A142C","#5D90BA","#431A3D"), # Colors for pathways
  threshold_colors = brewer.pal(4, "PiYG"), # Threshold-based colors
  title = "Volcano Plot",    # Plot title
  x_range = NULL,            # Range for x-axis, e.g., c(-5, 5)
  y_range = NULL,            # Range for y-axis, e.g., c(0, 10)
  label = TRUE               # Whether to label genes and pathways
) {
  # Ensure required columns are present in data
  data$Gene <- data[[gene.col]]
  data$logFC <- data[[logFC.col]]
  data$P.Value <- data[[pval.col]]

  # Merge selected genes into the main data
  merged_data <- merge(data, selected_genes, by.x = "Gene", by.y = "gsym", all.x = TRUE)
  merged_data$Pathway <- ifelse(is.na(merged_data$Pathway), "Other", merged_data$Pathway)
  merged_data$label <- ifelse(label & !is.na(merged_data$Pathway) & merged_data$Pathway != "Other", merged_data$Gene, "")

  # Apply threshold-based color mapping
  cols <- rep("grey", nrow(merged_data))
  cols[merged_data$P.Value < pvalCut & merged_data$logFC > logFCcut] <- threshold_colors[2]
  cols[merged_data$P.Value < pvalCut2 & merged_data$logFC > logFCcut2] <- threshold_colors[1]
  cols[merged_data$P.Value < pvalCut & merged_data$logFC < -logFCcut] <- threshold_colors[3]
  cols[merged_data$P.Value < pvalCut2 & merged_data$logFC < -logFCcut2] <- threshold_colors[4]
  merged_data$color_transparent <- adjustcolor(cols, alpha.f = 0.75)

  # Assign sizes based on thresholds
  merged_data$size <- rep(1, nrow(merged_data))
  merged_data$size[merged_data$P.Value < pvalCut & abs(merged_data$logFC) > logFCcut] <- 2
  merged_data$size[merged_data$P.Value < pvalCut2 & abs(merged_data$logFC) > logFCcut2] <- 4
  merged_data$size[merged_data$P.Value < pvalCut3 & abs(merged_data$logFC) > logFCcut3] <- 6

  # Define axis limits
  if (is.null(x_range)) x_range <- c(min(merged_data$logFC) - 1, max(merged_data$logFC) + 1)
  if (is.null(y_range)) y_range <- c(0, max(-log10(merged_data$P.Value)) * 1.1)

  # Construct the plot
  p <- ggplot(data = merged_data, aes(x = logFC, y = -log10(P.Value), color = Pathway)) +
    geom_point_rast(aes(size = size), alpha = 0.6, colour = merged_data$color_transparent) +
    labs(x = bquote(~Log[2]~"(fold change)"), 
         y = bquote(~-Log[10]~italic("FDR")), 
         title = title) +
    scale_x_continuous(limits = x_range) +
    scale_y_continuous(limits = y_range) +
    geom_vline(xintercept = c(-logFCcut, logFCcut), color = "grey40", 
               linetype = "longdash", lwd = 0.5) +
    geom_hline(yintercept = -log10(pvalCut), color = "grey40", 
               linetype = "longdash", lwd = 0.5) +
    theme_classic(base_size = 12) +
    theme(panel.grid = element_blank()) +
    scale_color_manual(values = pathway_colors) +
    guides(color = guide_legend(title = "Pathways"))

  # Add advanced mode lines
  if (plot_mode == "advanced") {
    p <- p +
      geom_vline(xintercept = c(-logFCcut2, logFCcut2), color = "grey40", 
                 linetype = "longdash", lwd = 0.5) +
      geom_hline(yintercept = -log10(pvalCut2), color = "grey40", 
                 linetype = "longdash", lwd = 0.5)
  }

  # Optionally highlight selected genes with circles and pathway-specific colors
  if (label) {
    p <- p +
      geom_point_rast(data = subset(merged_data, Gene %in% selected_genes$gsym),
                 aes(size = size), alpha = 1, shape = 1, stroke = 1, color = "black") +
      geom_text_repel(data = subset(merged_data, Gene %in% selected_genes$gsym),
                      aes(label = label, color = Pathway), max.overlaps = 100, size = 2.5,
                      box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines"))
  }

  return(p)
}
