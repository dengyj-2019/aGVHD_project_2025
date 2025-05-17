# ============ Code to Reproduce Fig. 6b–f (Ruxolitinib Response)====================
# Environment Setup
# ================================
library(Seurat)
library(dplyr)
library(ggplot2)
library(viridis)
library(reshape2)
library(ggrepel)
library(ggpubr)
library(tidyr)
library(tidyverse)
library(Matrix)
library(cowplot)
library(EnhancedVolcano)
library(circlize)
library(patchwork)
library(stringr)

# Set working directory and define color palettes
setwd("D:/Project/GVHDscRNA/20250430")
clustcol <- c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72", "#B17BA6",
              "#FF7F00", "#FDB462", "#E7298A", "#E78AC3", "#33A02C", "#B2DF8A",
              "#55A1B1", "#8DD3C7", "#A6761D", "#E6AB02", "#7570B3", "#BEAED4",
              "#666666", "#999999", "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3",
              "#808000", "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00")
gpcol <- c("#D9534F", "#55A1B1") # Red: before (SR); Blue: after (CR)

# Load processed Seurat object of CD3+ T cells
Tcell <- readRDS("Tcell_seurat_object.rds")

# Output directory
outdir <- "./out/"
dir.create(outdir)
prefix <- "Visualization"

# ================================
# Subsetting Metadata
# ================================
# Subset only "before" (Steroid-Refractory, SR) and "after" (Complete Response, CR) samples
Tcell <- Tcell[, Tcell$former_type %in% c('before', 'after')]  # before = SR, after = CR
Idents(Tcell) <- Tcell$Identity

# Further subset CD8+ Tem cells only
celltype <- c("CD8+TEM")
TB <- subset(Tcell, idents = celltype)
levels(TB) <- celltype

# ================================
# Fig. 6b: UMAP of CD3+ T Cells (Before vs. After)
# ================================

# UMAP with contour (SR group)
beforeTcell <- Tcell[, Tcell$former_type == "before"]
plot3 <- DimPlot(beforeTcell, label = FALSE, cols = clustcol, pt.size = 0.01) +
  labs(title = "SR") +
  theme(panel.border = element_rect(color = "black", size = 1),
        axis.text = element_blank(), axis.ticks = element_blank()) +
  geom_density_2d(data = data.frame(UMAP1 = Embeddings(beforeTcell, "umap")[, 1],
                                    UMAP2 = Embeddings(beforeTcell, "umap")[, 2]),
                  aes(x = UMAP1, y = UMAP2), color = "black")
ggsave(paste0(outdir, '/', prefix, '.SRumap.legend_with_contours.pdf'), plot3, width = 7.5, height = 6)

# UMAP with contour (CR group)
afterTcell <- Tcell[, Tcell$former_type == "after"]
plot4 <- DimPlot(afterTcell, label = FALSE, cols = clustcol, pt.size = 0.01) +
  labs(title = "CR") +
  theme(panel.border = element_rect(color = "black", size = 1),
        axis.text = element_blank(), axis.ticks = element_blank()) +
  geom_density_2d(data = data.frame(UMAP1 = Embeddings(afterTcell, "umap")[, 1],
                                    UMAP2 = Embeddings(afterTcell, "umap")[, 2]),
                  aes(x = UMAP1, y = UMAP2), color = "black")
ggsave(paste0(outdir, '/', prefix, '.CRumap.legend_with_contours.pdf'), plot4, width = 7.5, height = 6)

# ================================
# Fig. 6c: Bar Plot of T Cell Subtype Fractions
# ================================
freq.table <- prop.table(table(Tcell$Identity, Tcell$former_type), margin = 2)
tr <- melt(freq.table)
colnames(tr) <- c("Identity", "former_type", "fraction")
tr$Identity <- factor(tr$Identity, levels = c("CD8+TN", "CD8+TCM", "CD8+TEM", "CD8+TEFF",
                                              "CD4+TN", "CD4+TCM", "Th17", "Th1", "Treg",
                                              "MAIT", "IFN_T", "Vδ1_T", "Vδ2_T", "NKT"))
y.max <- max(tr$fraction) * 1.05
p <- ggplot(tr, aes(x = Identity, y = fraction, fill = former_type)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", width = 0.7) +
  theme_classic() +
  scale_fill_manual(values = gpcol) +
  ylab("Fraction of Cells") +
  ylim(0, y.max) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(paste0(outdir, '/', prefix, '.CellFraction_dodg.pdf'), p, width = 4.5, height = 4.5)

# ================================
# Fig. 6d: DEGs Between SR and CR CD8+ Tem for Metascape
# ================================
CD8TEM <- subset(Tcell, idents = "CD8+TEM")
diff_CD8TEM <- FindMarkers(CD8TEM, group.by = "former_type", ident.1 = "before", ident.2 = "after",
                           logfc.threshold = 0.2, min.pct = 0.3)
diff_CD8TEM <- diff_CD8TEM[diff_CD8TEM$p_val_adj < 0.01, ]
write.csv(diff_CD8TEM, "diff_CD8TEM.beforevsafter_sig.csv")

# ================================
# Fig. 6e: Pathway Activity (AUCell-like Program Scores)
# ================================
library(xlsx)
geneset <- read.xlsx("pathways.xlsx", sheetIndex = 1)
for (i in unique(geneset$term)) {
  subset_genes <- geneset %>% filter(term == i)
  TB <- AddModuleScore(TB, features = list(subset_genes$gene), name = i)
}
n_start <- ncol(TB@meta.data) - length(unique(geneset$term)) + 1
n_end <- ncol(TB@meta.data)
colnames(TB@meta.data)[n_start:n_end] <- str_replace(colnames(TB@meta.data)[n_start:n_end], "1$", "")

# Boxplots for each module

library(ggrastr)

library(dplyr)

data<- FetchData(TB,vars = c("former_type","CYTOKINE_SIGNALING","GC_RESISTANCE",
                             "memory_like_feature","T_cell_activation"))


colors <- c("#D9534F","#55A1B1")

# Convert novel_type to factor with desired levels order

data$former_type <- factor(data$former_type, levels = c("before", "after"))

comparisons <- list(c("before", "after"))

# Define plotting function
plotAUCellProgram <- function(var, title_name) {
  p <- ggplot(data, aes_string(x = "former_type", y = var, fill = "former_type", color = "former_type")) +
    theme_bw() + RotatedAxis() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(size = 12, angle = 0, vjust = 0.5, hjust = 0.5),
          axis.text.y = element_text(size = 12),
          plot.title = element_text(hjust = 0.5),
          legend.position = 'none') +
    labs(x = NULL, y = NULL, title = title_name) +
    geom_jitter_rast(col = "#00000033", pch = 19, cex = 2, position = position_jitter(0.2)) +
    geom_boxplot(position = position_dodge(0), color = 'black',
                 outlier.colour = NA, outlier.fill = NA, outlier.shape = NA) +
    scale_fill_manual(values = colors) +
    stat_compare_means(comparisons = comparisons, method = "wilcox.test", label = "p.signif")
  
  # Export to PDF
  pdf(paste0(outdir, '/', prefix, title_name, '.pdf'), width = 3, height = 4)
  print(p)
  dev.off()
}

# Generate all four plots for Fig. 6e
plotAUCellProgram("T_cell_activation", "T_cell_activation")
plotAUCellProgram("memory_like_feature", "memory_like_feature")
plotAUCellProgram("CYTOKINE_SIGNALING", "CYTOKINE_SIGNALING")
plotAUCellProgram("GC_RESISTANCE", "GC_RESISTANCE")


# ================================
# Fig. 6f:  Boxplot + Jitter Plot 
# Description: This script generates boxplots with overlaid jitter points for each specified gene.
# Author: [Your Name]
# Date: [Date]
# ================================

# Load required libraries
library(ggplot2)
library(ggpubr)

# Load expression data from Seurat object (TB)
# Ensure that 'TB' is a Seurat object with metadata including 'former_type' and the genes of interest
data <- FetchData(TB, vars = c("former_type", "IL2RA", "IL27RA", "STAT1", "NR3C1"))

# Convert 'former_type' to a factor variable with the specified order
data$former_type <- factor(data$former_type, levels = c("before", "after"))

# Define custom colors for the two groups (consistent, no transparency change)
colors <- c("before" = "#D9534F", "after" = "#55A1B1")  # Red for 'before', Blue for 'after'

# Define a function to create boxplots with mean line and small jittered points
plot_box_mean <- function(gene_name) {
  
  # Generate the plot using ggplot2
  p <- ggplot(data, aes_string(x = "former_type", y = gene_name, fill = "former_type")) +
    geom_boxplot(outlier.shape = NA, width = 0.4, color = "black") +  # Boxplot without outliers
    geom_jitter(position = position_jitter(width = 0.15, height = 0), 
                size = 0.5, alpha = 0.6, color = "black") +  # Small, semi-transparent jittered points
    stat_summary(fun = mean, geom = "crossbar", 
                 width = 0.4, color = "black", size = 0.8) +  # Black mean line only
    scale_fill_manual(values = colors) +  # Custom fill colors (no transparency change)
    labs(x = NULL, y = "Expression Level", title = gene_name) +  # Axis labels and title
    theme_minimal() +  # Clean theme for better visualization
    theme(
      legend.position = "none",  # Remove legend
      axis.text.x = element_text(size = 12),  # Customize X-axis text size
      axis.text.y = element_text(size = 12),  # Customize Y-axis text size
      axis.ticks = element_line(),  # Show axis ticks
      axis.title = element_text(size = 12, face = "bold"),  # Bold axis titles
      plot.title = element_text(hjust = 0.5, face = "bold"),  # Centered bold title
      panel.background = element_rect(color = "black", fill = NA),  # Black border around the plot
      panel.grid = element_blank()  # Remove grid lines
    ) +
    # Add statistical significance (Wilcoxon rank-sum test)
    stat_compare_means(comparisons = list(c("before", "after")), 
                       method = "wilcox.test", 
                       label = "p.signif")
  
  return(p)
}

# Define the genes of interest for plotting
genes_of_interest <- c("IL2RA", "IL27RA", "STAT1", "NR3C1")

# Generate plots for each gene using the defined function
plots <- lapply(genes_of_interest, plot_box_mean)

# Arrange the plots in a 2x2 grid
combined_plot <- ggarrange(
  plots[[1]], plots[[2]], 
  plots[[3]], plots[[4]], 
  ncol = 2, nrow = 2, align = "hv"
)

# Display the combined plot
print(combined_plot)

# Save the combined plot as a PDF
ggsave(
  filename = "Combined_Boxplot_Mean_Only_2x2.pdf",
  plot = combined_plot,
  device = "pdf",
  width = 8,
  height = 8
)