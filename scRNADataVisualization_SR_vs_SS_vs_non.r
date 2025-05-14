library(Seurat)
library(dplyr)
#install.packages("dplyr")
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
library("patchwork")
library("stringr")
setwd("D:\\Project\\GVHDscRNA\\20250430")
#Color
clustcol <- c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
              "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
              "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D",
              "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999",
              "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000",
              "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00")
gpcol <- c("#337AB7","#F0AD4E","#D9534F")
clustcol.order <- clustcol

#Data preparation
rds <- "Tcell_seurat_object_V1.rds"
Tcell <- readRDS(rds)
Tcell<-Tcell[, Tcell$novel_type %in% c('non','SS','SR')]
outdir <- "./out/"
dir.create(outdir)
prefix <- "Visualization"

CD8TEM <- Tcell[, Idents(Tcell) %in% c('CD8+TEM')]
unique(CD8TEM$Identity)
Idents(CD8TEM)<-CD8TEM$Identity
celltype <- c("CD8+TN","CD8+TCM","CD8+TEM","CD8+TEFF","CD4+TN","CD4+TCM", "Th17",
              "Th1","Treg","MAIT", "IFN_T","Vδ1_T", "Vδ2_T", "NKT")
TB <- subset(Tcell,idents = celltype)
levels(TB) <- celltype
levels(TB)

# ============== Differential Expressed gene list for Extended Data Fig. 5d =====================

# Perform differential expression between SR vs SS CD8+ Tem cells
diff_CD8TEM <- FindMarkers(
  CD8TEM,                      # Seurat object containing CD8+ Tem cells
  min.pct = 0.2,               # Only test genes expressed in at least 20% of cells
  logfc.threshold = 0.2,       # Only return genes with log2FC > 0.2
  group.by = "novel_type",     # Group by novel_type (e.g. SR, SS, non)
  ident.1 = "SR",              # Case group: SR-aGvHD
  ident.2 = "SS"               # Control group: SS-aGvHD
)

# Filter for statistically significant genes (adjusted P < 0.01)
diff_CD8TEM <- diff_CD8TEM[diff_CD8TEM$p_val_adj < 0.01, ]

# Save result for downstream dot plot and pathway analysis (e.g., GR pathway, T cell activation genes)
write.csv(diff_CD8TEM, file = "diff_CD8TEM.SRvs_SS_novel_type.csv")


# Perform differential expression between SS vs non-aGvHD CD8+ Tem cells
diff_CD8TEM <- FindMarkers(
  CD8TEM,
  min.pct = 0.2,
  logfc.threshold = 0.2,
  group.by = "novel_type",
  ident.1 = "SR",              # Case group: SS-aGvHD
  ident.2 = "non"              # Control group: non-aGvHD
)

# Filter for significant genes
diff_CD8TEM <- diff_CD8TEM[diff_CD8TEM$p_val_adj < 0.01, ]

# Save result for Fig. 5a–c dot plots
write.csv(diff_CD8TEM, file = "diff_CD8TEM.SRvs_non_novel_type.csv")


# ===================== Marker Genes for Metascape Pathway Enrichment (Extended Data Fig. 4f) =====================

# Subset only SR-aGvHD samples
SRTcell <- Tcell[, Tcell$novel_type %in% c("SR")]
DefaultAssay(SRTcell) <- "RNA"

# Identify cluster-specific markers across all T-cell subsets in SR samples
all.markers <- FindAllMarkers(
  SRTcell,
  only.pos = TRUE,            # Only return upregulated genes per cluster
  min.pct = 0.25,             # Minimum percent of cells expressing gene
  logfc.threshold = 0.4       # Log2FC threshold for significance
)

# Retain only significantly upregulated markers (adjusted P < 0.01)
significant.markers <- all.markers[all.markers$p_val_adj < 0.01, ]

# Export full list for Metascape functional enrichment input
write.csv(significant.markers, file = "SR_significant.markers.csv")




#Cell fraction
#===================================Fig 4b===========================

freq.table <- prop.table(x=table(Tcell$Identity, Tcell$novel_type),margin=2)
head(freq.table)
tr <- melt(freq.table)
colnames(tr) <- c("Identity","novel_type", "fraction")
tr$Identity<- factor(tr$Identity,levels=c("CD8+TN","CD8+TCM","CD8+TEM","CD8+TEFF","CD4+TN","CD4+TCM", 
                                          "Th17","Th1","Treg","MAIT", "IFN_T","Vδ1_T", "Vδ2_T", "NKT"))

#Part8.1 Non-stacked bar plot of cell-type proportions
tr$Identity<- factor(tr$Identity,levels=c("CD8+TN","CD8+TCM","CD8+TEM","CD8+TEFF","CD4+TN","CD4+TCM", 
                                          "Th17","Th1","Treg","MAIT", "IFN_T","Vδ1_T", "Vδ2_T", "NKT"))
y.max <- max(tr$fraction)*1.05
p <- ggplot(tr,aes(x=Identity,y=fraction,fill=novel_type))+
  geom_bar(stat="identity",position = "dodge",width=0.7,color = "black")+
  theme_classic()+ 
  scale_fill_manual(values = gpcol) +
  xlab ('') + ylab ('Fraction of Cells') +
  ylim(c(0,y.max))+
  theme (axis.line.x = element_line(colour = "transparent"),
         axis.ticks.x = element_blank(),
         axis.text.x = element_text(angle = 90,color = "black", hjust = 1),
         aspect.ratio=1,
         panel.background=element_rect(fill='transparent', color='black')) +
  guides(fill=guide_legend(title='',ncol = 1))


pdf(paste0(outdir,'/',prefix,'.CellFraction_dodg.pdf'),width=4.5, height=4.5)
p
dev.off()

#Dimplot
#-----------------------Fig 4a left--------------------------------------------------
p.legend <- DimPlot(Tcell, cols = clustcol, pt.size = 0.01, reduction = "umap")
pdf(paste0(outdir,'/',prefix,'.umap.legend.pdf'),width=8, height=6)
p.legend
dev.off()

#-----------------------Fig 4a right-----------------------------
nonTcell<-Tcell[, Tcell$novel_type %in% c('non')]

library(ggplot2)

umap_data <- data.frame(UMAP1 = nonTcell@reductions$umap@cell.embeddings[, 1],
                        UMAP2 = nonTcell@reductions$umap@cell.embeddings[, 2])

plot3 <- DimPlot(nonTcell, label = FALSE, cols = clustcol, pt.size = 0.01) +
  labs(x = "UMAP1", y = "UMAP2", title = "non") +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

plot3_with_contours <- plot3 + 
  geom_density_2d(data = umap_data, aes(x = UMAP1, y = UMAP2), color = "black")

pdf(paste0(outdir, '/', prefix, '.nonumap.legend_with_contours.pdf'), width = 7.5, height = 6)
print(plot3_with_contours)
dev.off()


SSTcell<-Tcell[, Tcell$novel_type %in% c('SS')]

library(ggplot2)

umap_data <- data.frame(UMAP1 = SSTcell@reductions$umap@cell.embeddings[, 1],
                        UMAP2 = SSTcell@reductions$umap@cell.embeddings[, 2])

plot3 <- DimPlot(SSTcell, label = FALSE, cols = clustcol, pt.size = 0.01) +
  labs(x = "UMAP1", y = "UMAP2", title = "SS") +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

plot3_with_contours <- plot3 + 
  geom_density_2d(data = umap_data, aes(x = UMAP1, y = UMAP2), color = "black")

pdf(paste0(outdir, '/', prefix, '.SSumap.legend_with_contours.pdf'), width = 7.5, height = 6)
print(plot3_with_contours)
dev.off()


SRTcell<-Tcell[, Tcell$novel_type %in% c('SR')]

library(ggplot2)

umap_data <- data.frame(UMAP1 = SRTcell@reductions$umap@cell.embeddings[, 1],
                        UMAP2 = SRTcell@reductions$umap@cell.embeddings[, 2])

plot3 <- DimPlot(SRTcell, label = FALSE, cols = clustcol, pt.size = 0.01) +
  labs(x = "UMAP1", y = "UMAP2", title = "SR") +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

plot3_with_contours <- plot3 + 
  geom_density_2d(data = umap_data, aes(x = UMAP1, y = UMAP2), color = "black")

pdf(paste0(outdir, '/', prefix, '.SRumap.legend_with_contours.pdf'), width = 7.5, height = 6)
print(plot3_with_contours)
dev.off()

#-------------------heatmap of top 30 genes---------------------------------------------

#-------------------Fig 4c---------------------------------------------
library(Seurat)
library(dplyr)
library(dittoSeq)
SRTcell<-Tcell[, Tcell$novel_type %in% c('SR')]
DefaultAssay(SRTcell) <- "RNA"
Idents(SRTcell) <- SRTcell$Identity
celltype_markers  <- FindAllMarkers(SRTcell, only.pos = TRUE, 
                                    min.pct = 0.3, logfc.threshold = 0.3)
significant.markers  <- celltype_markers [celltype_markers$p_val_adj < 0.01, ]
top30 = significant.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)

gene_cell_exp <- AverageExpression(SRTcell,
                                   features = top30$gene,
                                   group.by = 'Identity',
                                   slot = 'data') 
gene_cell_exp <- as.data.frame(gene_cell_exp$RNA)
write.csv(gene_cell_exp,file="markergene_exp.csv")

marker_exp <- read.csv("markergene_exp.csv", header = T, row.names = 1)
top30 <- read.csv("top.csv", header = T, row.names = 1)

# Note: 'top.csv' was manually edited to:
# set gene names as rownames
# 'top.csv' and "markergene_exp.csv" were edited to:
# assign unique cluster labels (e.g., 'cluster_a_CD8+TN') to allow row splitting in ComplexHeatmap.
group <- top30[,c(6,7)]
marker_exp <- t(scale(t(marker_exp),scale = T,center = T))
rownames(group) <- rownames(marker_exp)


library(ComplexHeatmap)
library(circlize)
my_color_range <- colorRamp2(c(-2, -1,0,1, 2), c('#559643',"#C7E0A1",'white',"#F2C1DA","#CB3885"))

p <- Heatmap(
  marker_exp,                    # The expression matrix to be plotted
  cluster_rows = FALSE,          # Do not cluster rows
  cluster_columns = FALSE,       # Do not cluster columns
  show_column_names = TRUE,      # Show column names
  show_row_names = TRUE,         # Show row (gene) names
  row_split = group$cluster,     # Split rows based on cluster assignment
  row_title = NULL,              # Do not display row group titles
  column_title = NULL,           # Do not display column title
  row_gap = unit(0, "mm"),       # Set spacing between row groups to zero
  heatmap_legend_param = list(
    title = "Exp"                # Title for the heatmap color legend
  ),
  col = my_color_range,          # Color scale for heatmap (e.g., from colorRamp2)
  column_names_gp = gpar(fontsize = 8),  # Set font size for column labels
  border = "black",              # Add a black border around the entire heatmap
  border_gp = gpar(col = "black", lwd = 2)  # Border styling: color and line width
)

genes <- c("TCF7","LEF1","CD69",#CD8TN
           "GZMK","HMGB1","HMGB2",
           "HMGN2","TIMD4","MKI67","CDKN3", "MCM5",#CD8TEM
           "CCL5","CX3CR1","GZMH","FGFBP2","GZMB","PRF1","GZMA",#CD8TEFF
           "IL4I1","KLRB1","KLRG1","IL7R","GPR65","IFNGR1", #MAITT
           "S100A4","S100A10","CD28","SLAMF1",#Th1
           "FOXP3","IL2RA","IKZF2","IL32",#Treg
           "IFIT1","IFIT3","IFIT2","ISG15",#IFNT
           "S1PR5","KIR2DL3","GZMM","HLA-DRB1", "CD300A",#gdT
           "KIR2DL4","KLRC1","XCL1","XCL2","GNLY"#NKT
)
genes <- as.data.frame(genes)


p1<-p + rowAnnotation(link = anno_mark(at = which(rownames(marker_exp) %in% genes$genes), 
                                       labels = genes$genes, labels_gp = gpar(fontsize = 4)))

pdf(paste0(outdir,'/',prefix,'.heatmap_top30.pdf'),width=4.5, height=6)
p1
dev.off()



#=============================== DotPlot Extended Data Fig. 4e======================================#
markers.to.plot <- c("IL18R1", "IL10RA", "IL7R", "IL6R", "IL2RB", "IL2RA","IL6ST",
                     "TNF", "LGALS9", "FALSG", "IL12RA","IL10", "TGFB1", 
                     "IL23A","IL22","IL21", "IL17A", "IL2", "IFNG")
celltype <- c("Treg","Th1","Th17","CD4+TCM","CD4+TN")

CD4T <- subset(Tcell, idents = celltype)
levels(CD4T) <- celltype
levels(CD4T)
pd1 <- DotPlot(CD4T, features = markers.to.plot)+coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=0.5))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'), guide = guide_legend(reverse = TRUE))

pdf(paste0(outdir,'/',prefix,'.dotplot.CD4T.pdf'),width = 3.5, height = 4)
pd1
dev.off()


#-------------AddModuleScore for gene signatures --------------------------
library(xlsx)
library(dplyr)
library(stringr)
library(Seurat)

# Step 1: Load gene sets from Excel file
# The Excel file should contain two columns: 'term' (gene set name) and 'gene' (gene symbol)
geneset <- read.xlsx("pathways.xlsx", sheetIndex = 1)

# Step 2: Loop through each gene set and calculate module score
for (term in unique(geneset$term)) {
  genes <- geneset %>% filter(term == !!term) %>% pull(gene)
  TB <- AddModuleScore(object = TB, features = list(genes), name = term)
}

# Step 3: Remove Seurat's default '1' suffix (e.g., 'Pathway1' → 'Pathway') from column names
n_term <- length(unique(geneset$term))
meta_cols <- (ncol(TB@meta.data) - n_term + 1):ncol(TB@meta.data)
colnames(TB@meta.data)[meta_cols] <- str_replace(colnames(TB@meta.data)[meta_cols], "1$", "")

##============================Boxplot of GR pathway socre==============================###
##============================Fig.4g==============================###
library(ggrastr)

library(dplyr)

data<- FetchData(TB,vars = c("novel_type","T_cell_activation","GR_PATHWAY"))


colors <- c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
            "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
            "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7")

colors <- c("#337AB7", "#F0AD4E", "#D9534F")

# Convert novel_type to factor with desired levels order

data$novel_type <- factor(data$novel_type, levels = c("non", "SS", "SR"))

comparisons <- list( c("CD8+TEM", "CD8+TN"),c("CD8+TEM", "CD8+TCM"),c("CD8+TEM", "CD8+TEFF"), 
                     c("CD8+TEM", "CD4+TN"), c("CD8+TEM", "CD4+TCM"), c("CD8+TEM", "Th17"),
                     c("CD8+TEM", "Th1"), c("CD8+TEM", "Treg"), c("CD8+TEM", "MAIT"),
                     c("CD8+TEM", "IFN_T"), c("CD8+TEM", "Vδ1_T"), c("CD8+TEM", "Vδ2_T"),
                     c("CD8+TEM", "NKT")
)

comparisons <- list(c("non", "SR"), c("SS", "SR"))

p <- ggplot(data, aes(x = novel_type, y = GR_PATHWAY, fill = novel_type, color = novel_type)) +
  theme_bw() + RotatedAxis() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = 12, angle = 0, vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'none') +
  labs(x = NULL, y = NULL, title = "GR_PATHWAY") +
  geom_jitter_rast(col = "#00000033", pch = 19, cex = 2, position = position_jitter(0.2)) +
  geom_boxplot(position = position_dodge(0)) +
  scale_fill_manual(values = colors) +  # Specify fill colors
  geom_boxplot(position = position_dodge(0), color = 'black',
               outlier.colour = NA, outlier.fill = NA, outlier.shape = NA)+
  stat_compare_means(comparisons = comparisons, method = "wilcox.test", label = "p.signif") 

pdf(paste0(outdir, '/', prefix, 'GR_PATHWAY socore.pdf'), width = 4, height = 4)
print(p)
dev.off()
png(paste0(outdir, '/', prefix, '_GR_PATHWAY_score.png'), width = 4, height = 4, units = "in", res = 300)
print(p)
dev.off()



# ===========================  Correlation of Pathway Module Scores =========================== #

# Load required libraries
library(ggplot2)
library(ggpubr)
library(ggExtra)

# Fetch relevant module scores from Seurat object
data <- FetchData(TB, vars = c("novel_type", "T_cell_activation", "GR_PATHWAY", "STAT1_TARGETS"))

# ====================================================================================
# Function to generate scatter plot with Spearman correlation and marginal densities
# ====================================================================================
plot_correlation <- function(x, y, x_label, y_label, output_file) {
  
  # Perform Spearman correlation test
  cor_test <- cor.test(x, y, method = "spearman")
  cor_value <- round(cor_test$estimate, 3)
  p_value <- signif(cor_test$p.value, 3)
  
  # Create data frame for plotting
  df <- data.frame(x = x, y = y)
  
  # Generate scatter plot with regression line and correlation coefficient
  p <- ggplot(df, aes(x = x, y = y)) +
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ x, color = "blue") +
    xlab(x_label) +
    ylab(y_label) +
    theme_bw() +
    stat_cor(method = 'spearman', aes(x = x, y = y),
             label.x.npc = "left", label.y.npc = "top") +
    labs(title = paste("Correlation:", cor_value, "P-value:", p_value))
  
  # Add marginal density plots
  p_with_marginals <- ggMarginal(p, type = "density",
                                 xparams = list(fill = "orange"),
                                 yparams = list(fill = "blue"))
  
  # Save to PDF
  pdf(file = output_file, width = 4, height = 4)
  print(p_with_marginals)
  dev.off()
}

# ====================================================================================
# Panel 1: GR_PATHWAY vs. T_cell_activation  Fig. 4d
# ====================================================================================
plot_correlation(
  x = data$GR_PATHWAY,
  y = data$T_cell_activation,
  x_label = "GR_PATHWAY",
  y_label = "T_cell_activation",
  output_file = "Correlation_GR_PATHWAY_vs_T_cell_activation.pdf"
)

# ====================================================================================
# Panel 2: STAT1_TARGETS vs. T_cell_activation Extended Data Fig. 7c upper panel
# ====================================================================================
plot_correlation(
  x = data$STAT1_TARGETS,
  y = data$T_cell_activation,
  x_label = "STAT1_TARGETS",
  y_label = "T_cell_activation",
  output_file = "Correlation_STAT1_TARGETS_vs_T_cell_activation.pdf"
)

# ====================================================================================
# Panel 3: STAT1_TARGETS vs. GR_PATHWAY  Extended Data Fig. 7c lower panel
# ====================================================================================
plot_correlation(
  x = data$STAT1_TARGETS,
  y = data$GR_PATHWAY,
  x_label = "STAT1_TARGETS",
  y_label = "GR_PATHWAY",
  output_file = "Correlation_STAT1_TARGETS_vs_GR_PATHWAY.pdf"
)


##=============Extended Data Fig. 6a Heatmap of GC-related Pathway Module Scores (Z-scored) ========================##

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Extract pathway scores for selected gene sets from Seurat object
data <- FetchData(SRTB, vars = c("Identity", "RESPONSE_TO_GC", "GC_THERAPY_UP",
                                 "GC_THERAPY_DN", "GR_PATHWAY", "GC_RESISTANCE"))

# Define pathway module names
gene <- c("RESPONSE_TO_GC", "GC_THERAPY_UP", "GC_THERAPY_DN", "GR_PATHWAY", "GC_RESISTANCE")

# Step 1: Calculate the average score for each pathway within each cell subset
avg_gene_expression <- data %>%
  group_by(Identity) %>%
  summarize(across(all_of(gene), mean, na.rm = TRUE))

# Step 2: Convert wide-format dataframe to long-format for ggplot
avg_gene_expression_long <- avg_gene_expression %>%
  pivot_longer(cols = -Identity, names_to = "Gene", values_to = "Mean_Expression")

# Step 3: Apply z-score normalization within each pathway (row-wise)
avg_gene_expression_long <- avg_gene_expression_long %>%
  group_by(Gene) %>%
  mutate(Mean_Expression = (Mean_Expression - mean(Mean_Expression)) / sd(Mean_Expression))

# Step 4: Set custom order of pathways for plotting
gene_order <- c("GC_RESISTANCE", "GR_PATHWAY", "GC_THERAPY_DN", "GC_THERAPY_UP", "RESPONSE_TO_GC")
avg_gene_expression_long$Gene <- factor(avg_gene_expression_long$Gene, levels = gene_order)

# Step 5: Generate heatmap of z-scored pathway activity
p <- ggplot(avg_gene_expression_long, aes(x = Identity, y = Gene, fill = Mean_Expression)) +
  geom_tile(color = "black", size = 0.2) +
  scale_fill_gradient2(low = "#043161", mid = "#F4F6F5", high = "#660020") +
  labs(x = "Cell Subtype", y = "Pathway Module", fill = "Z-score") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black")
  )

# Step 6: Export heatmap to PDF
pdf(paste0(outdir, '/', prefix, 'GC-pathway_heatmap.pdf'), width = 7, height = 3)
print(p)
dev.off()


# Volcano plot: highlight selected genes from GR signaling pathway in CD8+ Tem cells
###=================================Fig 4h======================================###

# Load differential expression results (CD8+ Tem cells, SR vs SS)
df <- read.csv("diff_CD8TEM.SRvs_SS_novel_type.csv", header = TRUE)

# Load selected gene list to be highlighted
# This list was manually curated as the intersection of DEGs (CD8+ Tem SR vs SS)
# and genes involved in the glucocorticoid receptor (GR) signaling pathway
gene_list <- read.csv("GR_pathway_diffgene.csv", header = TRUE)
gene <- as.character(gene_list$gene)  # Convert to character vector

# Use gene names as rownames in df for easier indexing
rownames(df) <- df$feature
gene_plot <- df[gene, ]  # Subset data for highlighted genes

# Load required libraries
library(ggplot2)
library(ggrepel)

# Customize thresholds
y_cutoff <- 75
max_overlaps <- 40

# Construct volcano plot
p <- ggplot(df, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_hline(aes(yintercept = 2), color = "#999999", linetype = "dashed", size = 1) +      # P-value cutoff line
  geom_vline(aes(xintercept = 0.1), color = "#999999", linetype = "dashed", size = 1) +    # Fold-change threshold
  geom_vline(aes(xintercept = -0.1), color = "#999999", linetype = "dashed", size = 1) +
  
  # Background genes (grey)
  geom_point(data = df[!(rownames(df) %in% rownames(gene_plot)), ],
             stroke = 0.5, size = 1, shape = 16, color = "grey", alpha = 0.4) +
  
  # Highlighted GR-pathway genes (green)
  geom_point(data = df[rownames(df) %in% rownames(gene_plot), ],
             stroke = 0.5, size = 2, shape = 16, color = "olivedrab3") +
  
  labs(x = "Log2 fold change",
       y = "-Log10(P-value)",
       title = "") +
  
  # Theme settings
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1, colour = "black"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12, color = "black"),
        plot.title = element_text(hjust = 0.5, size = 16),
        legend.position = "none") +
  
  # Add gene labels for GR-pathway genes
  geom_text_repel(data = gene_plot, aes(label = feature),
                  color = "black", size = 4, fontface = "italic",
                  arrow = arrow(ends = "first", length = unit(0.01, "npc")),
                  box.padding = 0.2,
                  point.padding = 0.3,
                  segment.color = 'black',
                  segment.size = 0.3,
                  force = 1, max.iter = 3e3,
                  max.overlaps = max_overlaps) +
  
  # Limit x and y axis ranges
  coord_cartesian(xlim = c(-0.5, 0.5), ylim = c(0, y_cutoff))

# Export the plot
pdf("velcano_SR_vs_SS_CD8TEM_GR_diffgene.pdf", width = 5, height = 4)
print(p)
dev.off()


##=============Heatmap of Z-Scored Expression for Selected DEGs in CD8⁺ Tem Cells ========================##
##=============Extended Data Fig. 5a–c, left panels ========================##

# Load required libraries
library(ComplexHeatmap)
library(circlize)

# Define selected genes of interest (differentially expressed genes)

# T cell activation-related genes (Extended Data Fig. 5a)
gene <- c("HMGA1", "HMGB1", "HMGB2", "CXCR6", "CCL3", "CCL3L1", "CCL23", "CD38", "CD27",
          "FCGR3A", "TNFRSF9", "COTL1", "TBX21", "EOMES", "GZMK", "GZMB", "GZMA")

# Cell cycle-related genes (Extended Data Fig. 5b)
# gene <- c("E2F4", "CCNB1", "CCNA2", "CENPE", "CCNB2", "CENPW", "MCMBP", "MCM5", "CDK9", "MKI67")

# Memory-like phenotype-related genes (Extended Data Fig. 5c)
# gene <- c("ITGB7", "ITGB1", "IL7R", "CISH", "HLA-DRA", "NCR3", "HLA-DRB1", "CX3CR1",
#           "LEF1", "CCR7", "S1PR1", "TCF7", "KLF2", "SELL")


# Step 1: Calculate average expression per group (log-normalized)
gene_cell_exp <- AverageExpression(
  CD8TEM,
  features = gene,
  group.by = "novel_type",
  slot = "data"  # Use log-normalized expression
)
gene_cell_exp <- as.data.frame(gene_cell_exp$RNA)

# Step 2: Perform Z-score normalization per gene (row-wise)
marker_exp <- t(scale(t(gene_cell_exp), center = TRUE, scale = TRUE))

# Step 3: Define heatmap color gradient
my_color_range <- colorRamp2(c(-2, 0, 2), c("#043161", "#F4F6F5", "#660020"))

# Step 4: Reorder columns to display 'SS' on the left and 'SR' on the right
marker_exp <- marker_exp[, c("non", "SS", "SR")]

# Step 5: Create ComplexHeatmap object
p <- Heatmap(
  marker_exp,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_column_names = TRUE,
  show_row_names = TRUE,
  column_title = NULL,
  heatmap_legend_param = list(title = " "),  # Legend without title
  col = my_color_range,
  border = "black",
  rect_gp = gpar(col = "black", lwd = 1),
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10)
)

# Step 6: Export heatmap to PDF
pdf(
  paste0(outdir, '/', prefix, '.CD8_TEM_heatmap.Diffgenes.T_cell_activation.pdf'),
  width = 3,
  height = 4
)
p
dev.off()


# Volcano plot: highlight selected genes from T cell activation, Cell_cycle, and Memory_like_phenotype
#in CD8+ Tem cells of SR vs. SS or SS vs. non
###====================Extended Data Fig. 5a–c, right panels=============================###

# Load differential expression results (CD8+ Tem cells, SR vs SS)
df <- read.csv("diff_CD8TEM.SRvs_SS_novel_type.csv", header = TRUE)
# Load differential expression results (CD8+ Tem cells, SS vs non)
#df <- read.csv("diff_CD8TEM.SSvs_non_novel_type.csv", header = TRUE)

# Load selected gene list to be highlighted
# These lists were manually curated as the intersection of DEGs (CD8+ Tem SR vs SS or SS vs non)
# and genes involved in the T cell activation, Cell_cycle, or Memory_like_phenotype
gene_list <- read.csv("DEG_list_T_cell_activation.csv", header = TRUE)
#gene_list <- read.csv("DEG_list_Cell_cycle.csv", header = TRUE)
#gene_list <- read.csv("DEG_list_Memory_like_phenotype.csv", header = TRUE)
gene <- as.character(gene_list$gene)  # Convert to character vector

# Use gene names as rownames in df for easier indexing
rownames(df) <- df$feature
gene_plot <- df[gene, ]  # Subset data for highlighted genes

# Load required libraries
library(ggplot2)
library(ggrepel)

# Adjust thresholds to optimize visual presentation
y_cutoff <- 75
max_overlaps <- 40

# Construct volcano plot
p <- ggplot(df, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_hline(aes(yintercept = 2), color = "#999999", linetype = "dashed", size = 1) +      # P-value cutoff line
  geom_vline(aes(xintercept = 0.1), color = "#999999", linetype = "dashed", size = 1) +    # Fold-change threshold
  geom_vline(aes(xintercept = -0.1), color = "#999999", linetype = "dashed", size = 1) +
  
  # Background genes (grey)
  geom_point(data = df[!(rownames(df) %in% rownames(gene_plot)), ],
             stroke = 0.5, size = 1, shape = 16, color = "grey", alpha = 0.4) +
  
  # Highlighted selected genes (green)
  geom_point(data = df[rownames(df) %in% rownames(gene_plot), ],
             stroke = 0.5, size = 2, shape = 16, color = "olivedrab3") +
  
  labs(x = "Log2 fold change",
       y = "-Log10(P-value)",
       title = "") +
  
  # Theme settings
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1, colour = "black"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12, color = "black"),
        plot.title = element_text(hjust = 0.5, size = 16),
        legend.position = "none") +
  
  # Add gene labels for selected genes
  geom_text_repel(data = gene_plot, aes(label = feature),
                  color = "black", size = 4, fontface = "italic",
                  arrow = arrow(ends = "first", length = unit(0.01, "npc")),
                  box.padding = 0.2,
                  point.padding = 0.3,
                  segment.color = 'black',
                  segment.size = 0.3,
                  force = 1, max.iter = 3e3,
                  max.overlaps = max_overlaps) +
  
  # Limit x and y axis ranges
  coord_cartesian(xlim = c(-0.5, 0.5), ylim = c(0, y_cutoff))

# Export the plot
pdf("velcano_SR_vs_SS_CD8TEM_diffgene.pdf", width = 5, height = 4)
print(p)
dev.off()


##============================Boxplot of GR pathway score==============================###
##============================Extended Data Fig. 6b right==============================###
SRTB<-TB[,TB$novel_type %in% c('SR')]
data<- FetchData(SRTB,vars = c("Identity","T_cell_activation","GR_PATHWAY"))

colors <- c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
            "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
            "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7")


comparisons <- list( c("CD8+TEM", "CD8+TN"),c("CD8+TEM", "CD8+TCM"),c("CD8+TEM", "CD8+TEFF"), 
                     c("CD8+TEM", "CD4+TN"), c("CD8+TEM", "CD4+TCM"), c("CD8+TEM", "Th17"),
                     c("CD8+TEM", "Th1"), c("CD8+TEM", "Treg"), c("CD8+TEM", "MAIT"),
                     c("CD8+TEM", "IFN_T"), c("CD8+TEM", "Vδ1_T"), c("CD8+TEM", "Vδ2_T"),
                     c("CD8+TEM", "NKT")
)


p <- ggplot(data, aes(x = Identity, y =GR_PATHWAY, fill = Identity, color = Identity)) +
  theme_bw() + RotatedAxis() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = 12, angle = 0, vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'none') +
  labs(x = NULL, y = NULL, title = "GR_PATHWAY") +
  geom_jitter_rast(col = "#00000033", pch = 19, cex = 2, position = position_jitter(0.2)) +
  geom_boxplot(position = position_dodge(0)) +
  scale_fill_manual(values = colors) +  # Specify fill colors
  geom_boxplot(position = position_dodge(0), color = 'black',
               outlier.colour = NA, outlier.fill = NA, outlier.shape = NA)+
  stat_compare_means(comparisons = comparisons, method = "wilcox.test", label = "p.signif") 

pdf(paste0(outdir, '/', prefix, 'GR_PATHWAY across cellsubset.pdf'), width = 10, height = 6)
print(p)
dev.off()

ggsave(
  filename = paste0(outdir, '/', prefix, '_GR_PATHWAY_across_cellsubset.png'),
  plot = p, width = 10, height = 6, units = "in", dpi = 300, bg = "white"
)

##============================Extended Data Fig. 6b left==============================###
##============================Feautre of GR pathway socre==============================###

p.legend<- FeaturePlot(SRTB,features = "GR_PATHWAY",
                       order = T,cols = viridis(256))
pdf(paste0(outdir,'/',prefix,'SR_GR_PATHWAY.featureplot.pdf'),width=6, height=5)
p.legend
dev.off()



# == Draw violin plot of NR3C1 expression across T cell subsets of SR== #
##============================Fig. 4e==============================###

# Load expression data for plotting
data <- FetchData(SRTcell, vars = c("Identity", "NR3C1"))

# Define custom colors for group fill
colors <- c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
            "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
            "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7")


# Set group order for proper plotting
data$Identity <- factor(data$Identity, levels = c("CD8+TN","CD8+TCM","CD8+TEM","CD8+TEFF", "MAIT", "CD4+TN","CD4+TCM",
                                                  "Th17","Th1","Treg", "IFN_T","Vδ1_T", "Vδ2_T", "NKT"))

# Create violin plot
p <- ggplot(data, aes(x = Identity, y = NR3C1, fill = Identity)) +
  geom_violin(trim = TRUE) +
  geom_jitter(position = position_jitter(0.2), size = 0.02, color = "black") +
  # Optionally add boxplot over violin
  # geom_boxplot(width = 0.2, position = position_dodge(0), color = 'black', outlier.colour = NA) +
  scale_fill_manual(values = colors) +
  labs(x = NULL, y = "Expression Level", title = "NR3C1") +
  theme(
    legend.position = "none",
    panel.spacing = unit(0, "lines"),
    axis.line = element_blank(),
    panel.background = element_rect(fill = NA, color = 'black'),
    strip.background = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.ticks.x = element_line(),
    axis.ticks.y = element_line(),
    axis.ticks.length = unit(0.2, 'cm'),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold")
  ) +
  # Add Wilcoxon test p-values for selected pairwise comparisons
  stat_compare_means(
    comparisons = list(c("CD8+TEM", "CD8+TN"), c("CD8+TEM", "CD8+TCM"), c("CD8+TEM", "CD8+TEFF"),
                       c("CD8+TEM", "CD4+TN"), c("CD8+TEM", "CD4+TCM"), c("CD8+TEM", "Th17"),
                       c("CD8+TEM", "Th1"), c("CD8+TEM", "Treg"), c("CD8+TEM", "MAIT"),
                       c("CD8+TEM", "IFN_T"), c("CD8+TEM", "Vδ1_T"), c("CD8+TEM", "Vδ2_T"),
                       c("CD8+TEM", "NKT")),
    method = "wilcox.test",
    label = "p.signif"
  )

# Display the plot
print(p)

# Save plot as PDF
ggsave(
  filename = "VlnPlot_Tcell_NR3C1.pdf",
  plot = p,
  device = "pdf",
  width = 10,
  height = 6
)



# == Violin plot of selected genes (e.g., CD28) across non, SS, and SR groups == #
# =========================== Extended Data Fig. 5e ============================ #

# Load expression data for plotting
data <- FetchData(CD8TEM, vars = c("novel_type", "CD28"))

# Define custom colors for group fill
colors <- c("#337AB7", "#F0AD4E", "#D9534F")  # non, SS, SR

# Set group order for proper plotting
data$novel_type <- factor(data$novel_type, levels = c("non", "SS", "SR"))

# Create violin plot
p <- ggplot(data, aes(x = novel_type, y = CD28, fill = novel_type)) +
  geom_violin(trim = TRUE) +
  geom_jitter(position = position_jitter(0.2), size = 0.02, color = "black") +
  # Optionally add boxplot over violin
  # geom_boxplot(width = 0.2, position = position_dodge(0), color = 'black', outlier.colour = NA) +
  scale_fill_manual(values = colors) +
  labs(x = NULL, y = "Expression Level", title = "CD28") +
  theme(
    legend.position = "none",
    panel.spacing = unit(0, "lines"),
    axis.line = element_blank(),
    panel.background = element_rect(fill = NA, color = 'black'),
    strip.background = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.ticks.x = element_line(),
    axis.ticks.y = element_line(),
    axis.ticks.length = unit(0.2, 'cm'),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold")
  ) +
  # Add Wilcoxon test p-values for selected pairwise comparisons
  stat_compare_means(
    comparisons = list(c("non", "SR"), c("SS", "SR")),
    method = "wilcox.test",
    label = "p.signif"
  )

# Display the plot
print(p)

# Save plot as PDF
ggsave(
  filename = "VlnPlot_CD8TEM_CD28.pdf",
  plot = p,
  device = "pdf",
  width = 2.5,
  height = 2.5
)


# ==================Extended Data Fig. 4a Expression analysis of RPS4Y1 across samples ================== #

library(Seurat)
library(ggplot2)

# Step 1: Get all sample IDs
sample_ids <- unique(Tcell$orig.ident)

# Step 2: Loop through each sample to compute expression ratio and save FeaturePlot
for (sample in sample_ids) {
  
  # Subset Seurat object for current sample
  sample_obj <- Tcell[, Tcell$orig.ident == sample]
  
  # Step 3: Extract raw expression values of RPS4Y1 (RNA assay, count slot)
  RPS4Y1_expr <- sample_obj[["RNA"]]@counts["RPS4Y1", ]
  
  # Step 4: Calculate percentage of cells expressing RPS4Y1
  expressing_cells <- sum(RPS4Y1_expr > 0)
  total_cells <- length(RPS4Y1_expr)
  expression_ratio <- expressing_cells / total_cells
  
  # Step 5: Output result
  cat("Sample:", sample, "- RPS4Y1 expression ratio:", round(expression_ratio * 100, 2), "%\n")
  
  # Step 6: Generate FeaturePlot for RPS4Y1 expression
  p <- FeaturePlot(sample_obj, features = "RPS4Y1", order = TRUE) +
    ggtitle(paste0("RPS4Y1 Expression in Sample ", sample))
  
  # Step 7: Save FeaturePlot to PDF
  ggsave(filename = paste0("FeaturePlot_RPS4Y1_", sample, ".pdf"),
         plot = p, width = 10, height = 10)
}



# =========================== Extended Data Fig. 4c & 4d: Code Availability ===========================
# Full pipeline including TCR-seq and scRNA-seq integration, preprocessing, and figure generation

# =========================== Load Required Libraries ===========================
library(Seurat)
library(scRepertoire)
library(dplyr)
library(ggplot2)
library(ggsci)
library(stringr)


# =========================== Load scRNA-seq Object ===========================
Tcell <- readRDS("Tcell_seurat_object_V1.rds")  # Replace with actual filename

# =========================== Load TCR Contig Files (Manually Read) ===========================
C1 <- read.csv('C1_filtered_contig_annotations.csv')
C2 <- read.csv('C2_filtered_contig_annotations.csv')
C6 <- read.csv('C6_filtered_contig_annotations.csv')
C7 <- read.csv('C7_filtered_contig_annotations.csv')
C9 <- read.csv('C9_filtered_contig_annotations.csv')
SS1 <- read.csv('SS1_filtered_contig_annotations.csv')
SS2 <- read.csv('SS2_filtered_contig_annotations.csv')
SS3 <- read.csv('SS3_filtered_contig_annotations.csv')
P1 <- read.csv('P1_filtered_contig_annotations.csv')
P2 <- read.csv('P2_filtered_contig_annotations.csv')
P3 <- read.csv('P3_filtered_contig_annotations.csv')
P4 <- read.csv('P4_filtered_contig_annotations.csv')
P5 <- read.csv('P5_filtered_contig_annotations.csv')
P7 <- read.csv('P7_filtered_contig_annotations.csv')
P8 <- read.csv('P8_filtered_contig_annotations.csv')
P1_remission <- read.csv('P1_remission_filtered_contig_annotations.csv')
P8_remission <- read.csv('P8_remission_filtered_contig_annotations.csv')

TCR_list <- list(C1,C2,C6,C7,C9,SS1,SS2,SS3,P1,P2,P3,P4,P5,P7,P8,P1_remission,P8_remission)

# =========================== Combine TCR Data Using scRepertoire ===========================
data_tcr <- combineTCR(
  TCR_list,
  ID = c("C1","C2","C6","C7","C9","SS1","SS2","SS3","P1",
         "P2","P3","P4","P5","P7","P8","P1_remission","P8_remission"),
  samples = c("non","non","non","non","non","SS","SS","SS","SR",
              "SR","SR","SR","SR","SR","SR","CR","CR"),
  removeNA = FALSE,
  removeMulti = FALSE,
  filterMulti = FALSE
)

# =========================== Standardize Barcode Format ===========================
# Strip prefixes, then reappend sample names to match Seurat cell barcodes
for (sample in names(data_tcr)) {
  data_tcr[[sample]][, 1] <- gsub(paste0("^", sample, "_"), "", data_tcr[[sample]][, 1])
}
samples <- c("C1","C2","C6","C7","C9","SS1","SS2","SS3","P1",
             "P2","P3","P4","P5","P7","P8","P1_remission","P8_remission")
for (i in 1:17) {
  data_tcr[[i]]$barcode <- paste0(data_tcr[[i]]$barcode, "_", samples[i])
}

# =========================== Integrate TCR and RNA Data ===========================
scTCR_RNA <- combineExpression(
  data_tcr,
  Tcell,
  cloneCall = "aa",
  cloneTypes = c(Single = 1, Small = 3, Medium = 10, Large = 30, Hyperexpanded = 100),
  proportion = FALSE
)

# =========================== Extended Data Fig. 4c ===========================
# CloneType distribution across T-cell subsets

scTCR_RNA$novel_type <- factor(scTCR_RNA$novel_type, levels = c("non", "SS", "SR", "latter"))
SRTCR_RNA <- scTCR_RNA[, scTCR_RNA$novel_type == "SR"]

levels_order <- c("CD8+TN", "CD8+TCM", "CD8+TEM", "CD8+TEFF", 
                  "CD4+TN", "CD4+TCM", "Th17", "Th1", "Treg", 
                  "MAIT", "IFN_T", "Vδ1_T", "Vδ2_T", "NKT")
SRTCR_RNA$ident <- factor(SRTCR_RNA$Identity, levels = levels_order)

p13 <- occupiedscRepertoire(SRTCR_RNA, x.axis = "ident", label = FALSE, proportion = TRUE) +
  scale_fill_npg() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

ggsave("ExtendedDataFig4c_cloneType_distribution.pdf", plot = p13, width = 7, height = 5)


# =========================== Extended Data Fig. 4d ===========================
# Clonal expansion UMAP overlay (clone size > 1)

scTCR_RNA$novel_type <- factor(scTCR_RNA$novel_type, levels = c("non", "SS", "SR", "latter"))

p11 <- clonalOverlay(
  scTCR_RNA,
  reduction = "umap",
  freq.cutpoint = 10,
  bins = 10,
  facet = "novel_type"
) +
  scale_color_manual(values = clustcol) +
  geom_density_2d(color = "black", size = 0.2)

pdf("ExtendedDataFig4d_clonalOverlay_plot.pdf", width = 7, height = 6)
print(p11)
dev.off()


# ===================== CellChat Setup and Object Generation ===================== #
# Description: This script subsets T cell Seurat object into "non", "SS", and "SR" groups,
# then creates CellChat objects with appropriate pre-processing steps for each.

# Load required libraries
library(CellChat)
library(Seurat)

# ===================== Step 1: Subset Seurat object by group ===================== #
# Extract "non", "SS", and "SR" samples from T cell Seurat object
non <- subset(Tcell, novel_type == 'non')
SS <- subset(Tcell, novel_type == 'SS')
SR <- subset(Tcell, novel_type == 'SR')

# ===================== Step 2: Define CellChat object creation function ===================== #
prepare_cellchat <- function(seurat_obj, assay_name = "RNA") {
  DefaultAssay(seurat_obj) <- assay_name
  expr_matrix <- GetAssayData(seurat_obj, layer = "data")
  meta_info <- seurat_obj@meta.data[, c("novel_type", "Identity")]
  colnames(meta_info) <- c("group", "labels")
  meta_info <- meta_info[colnames(expr_matrix), ]  # Ensure alignment
  chat <- createCellChat(object = expr_matrix, meta = meta_info, group.by = "labels")
  return(chat)
}

# ===================== Step 3: Create CellChat objects for each group ===================== #
non.cellchat <- prepare_cellchat(non)
SS.cellchat <- prepare_cellchat(SS)
SR.cellchat <- prepare_cellchat(SR)

# ===================== Step 4: Assign ligand-receptor database ===================== #
CellChatDB <- CellChatDB.human
non.cellchat@DB <- CellChatDB
SS.cellchat@DB <- CellChatDB
SR.cellchat@DB <- CellChatDB

# ===================== Step 5: Preprocessing for CellChat objects ===================== #
preprocess_cellchat <- function(chat_obj, workers = 32) {
  options(future.globals.maxSize = 100 * 1024^3)  # Allow up to 100 GB for future objects
  chat_obj <- subsetData(chat_obj)
  future::plan("multisession", workers = workers)
  chat_obj <- identifyOverExpressedGenes(chat_obj)
  chat_obj <- identifyOverExpressedInteractions(chat_obj)
  chat_obj <- computeCommunProb(chat_obj, type = "triMean")
  chat_obj <- filterCommunication(chat_obj, min.cells = 10)
  chat_obj <- computeCommunProbPathway(chat_obj)
  chat_obj <- aggregateNet(chat_obj)
  return(chat_obj)
}

# ===================== Step 6: Run preprocessing ===================== #
non.cellchat <- preprocess_cellchat(non.cellchat)
SS.cellchat <- preprocess_cellchat(SS.cellchat)
SR.cellchat <- preprocess_cellchat(SR.cellchat)

# ===================== Extended Data Fig. 4g: Interaction Heatmap ===================== #
# Plot interaction count heatmap for SR group
pdf("ExtendedDataFig4g_SR_InteractionCount_pheatmap1.pdf", width = 6, height = 5)
pheatmap::pheatmap(
  SR.cellchat@net$count,
  border_color = "black",
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  fontsize = 10,
  display_numbers = TRUE,
  number_color = "black",
  number_format = "%.0f"
)
dev.off()

# ===================== One-step CellChat function for other group analysis ===================== #
KS_cellchat <- function(input_obj,
                        assay = NULL,
                        group.by = NULL,
                        workers,
                        species = c('human', 'mouse'),
                        CellChatDB.use = NULL,
                        PPIuse = FALSE,
                        type = "triMean",
                        min.cells = 10) {
  
  cellchat.obj = createCellChat(input_obj, assay = assay, group.by = group.by)
  
  if (species == 'human') {
    CellChatDB <- CellChatDB.human
    ppi <- PPI.human
  } else {
    CellChatDB <- CellChatDB.mouse
    ppi <- PPI.mouse
  }
  
  if (is.null(CellChatDB.use)) {
    cellchat.obj@DB <- CellChatDB
  } else {
    CellChatDB <- subsetDB(CellChatDB, search = CellChatDB.use, key = "annotation")
    cellchat.obj@DB <- CellChatDB
  }
  
  cellchat.obj <- subsetData(cellchat.obj)
  future::plan("multisession", workers = workers)
  cellchat.obj <- identifyOverExpressedGenes(cellchat.obj)
  cellchat.obj <- identifyOverExpressedInteractions(cellchat.obj)
  
  if (!PPIuse) {
    cellchat.obj <- computeCommunProb(cellchat.obj, type = type)
  } else {
    cellchat.obj <- projectData(cellchat.obj, ppi)
    cellchat.obj <- computeCommunProb(cellchat.obj, raw.use = FALSE, type = type)
  }
  
  cellchat.obj <- filterCommunication(cellchat.obj, min.cells = min.cells)
  cellchat.obj <- computeCommunProbPathway(cellchat.obj)
  cellchat.obj <- aggregateNet(cellchat.obj)
  return(cellchat.obj)
}

# ===================== Re-run CellChat analysis on SS and SR with wrapper ===================== #
SS.cellchat <- KS_cellchat(SS, assay = 'RNA', group.by = "Identity", workers = 40, species = 'human')
SR.cellchat <- KS_cellchat(SR, assay = 'RNA', group.by = "Identity", workers = 40, species = 'human')

# ===================== Multi-group Comparison and Extended Data Fig. 4h ===================== #
# Merge SS and SR CellChat objects
object.list <- list(SS = SS.cellchat, SR = SR.cellchat)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

# Step 1: Identify DEGs between SR and SS groups
pos.dataset <- "SS"
features.name <- paste0(pos.dataset, ".merged")
cellchat <- identifyOverExpressedGenes(
  cellchat,
  group.dataset = "datasets",
  pos.dataset = pos.dataset,
  features.name = features.name,
  only.pos = FALSE,
  thresh.pc = 0.1,
  thresh.fc = 0.05,
  thresh.p = 0.05,
  group.DE.combined = TRUE
)

# Step 2: Map DEGs to ligand-receptor pairs
net <- netMappingDEG(cellchat, features.name = features.name, variable.all = TRUE)

# Step 3: Subset communication with upregulated signals in SR and downregulated in SS
net.up <- subsetCommunication(cellchat, net = net, datasets = "SR", ligand.logFC = 0.05)
net.down <- subsetCommunication(cellchat, net = net, datasets = "SS", ligand.logFC = -0.05)

# Step 4: Extract gene subsets (optional)
gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

# Step 5: Prepare LR pairs for plotting
pairLR.use.up <- net.up[, "interaction_name", drop = FALSE]
pairLR.use.down <- net.down[, "interaction_name", drop = FALSE]

# Step 6: Generate bubble plots targeting CD8+TEM
gg1 <- netVisual_bubble(
  cellchat,
  pairLR.use = pairLR.use.up,
  targets.use = "CD8+TEM",
  comparison = c(1, 2),
  angle.x = 90,
  remove.isolate = TRUE,
  title.name = paste0("Up-regulated signaling in ", names(object.list)[2])
)

gg2 <- netVisual_bubble(
  cellchat,
  pairLR.use = pairLR.use.down,
  targets.use = "CD8+TEM",
  comparison = c(1, 2),
  angle.x = 90,
  remove.isolate = TRUE,
  title.name = paste0("Down-regulated signaling in ", names(object.list)[2])
)

# Step 7: Combine and export plot
pdf("ExtendedDataFig4h_LR_Diff_BubblePlot_CD8TEM.pdf", width = 10, height = 5)
(gg1 | gg2) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")
dev.off()
