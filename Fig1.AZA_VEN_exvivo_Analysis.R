##########################################################################
### Process RNA-seq Data and Perform Drug Profiling Analysis
##########################################################################

# Normalize raw counts to TPM (Transcripts Per Million)
count <- read.csv("../02.data/File_8_RNA_seq_Raw_Reads_163S_4Healthy.csv", header = TRUE)
colnames(count)[1] <- "id2"
gene_length <- read.csv("../02.data/gencode.v22_geneLength.csv", header = TRUE)
count <- merge(gene_length, count, by = "id2")
row.names(count) <- count$id2
count <- count[, -c(1:5)]

row.names(gene_length) <- gene_length$id2
gene_length <- gene_length[row.names(count), ]
Length <- as.vector(gene_length$length)

# TPM normalization
TPM <- t(t(count / t(Length)) * 1e6 / colSums(count / t(Length)))
tpm <- cbind(gene_length, TPM)
tpm <- tpm[, -c(1:3)]
#write.csv(tpm, "../02.data/File_8_RNA_seq_TPM_163S_4Healthy.csv")

# Load drug profiling data
library(readxl)
drug_response <- read.csv("../02.data/File_3_Drug_response_sDSS_164S_17Healthy.csv", header = TRUE, check.names = FALSE)

# Identify overlapped samples between RNA-seq and Drug profiling datasets
samples_overlap <- intersect(colnames(tpm), colnames(drug_response))
tpm <- tpm[, c("Symbol", "Type", samples_overlap)]
drug_response <- drug_response[, c("Drug_name", samples_overlap)]

##########################################################################
### Unsupervised Clustering of Drug Profiling Data (Venetoclax and Azacitidine)
##########################################################################

library(pheatmap)
drug_data <- read.csv("../02.data/File_3_Drug_response_sDSS_164S_17Healthy.csv", header = TRUE)
drug_data <- subset(drug_data, Drug_name %in% c("Azacitidine", "Venetoclax"))

sample_ids <- read.csv("../02.data/Functional_Precision_Medicine_RNAseq_samples.csv", header = TRUE)
overlapping_ids <- data.frame(intersect(colnames(drug_data), sample_ids$sample))
colnames(overlapping_ids) <- "sample"

drug_data <- drug_data[, c("Drug_name", as.character(overlapping_ids$sample))]
row.names(drug_data) <- drug_data$Drug_name
drug_data <- drug_data[, -1]
drug_data <- drug_data[, colSums(is.na(drug_data)) == 0]

# Generate heatmap
core_matrix <- as.matrix(drug_data)
bk <- unique(c(seq(-2, 2, length = 100)))
heatmap_output <- pheatmap(core_matrix,
                           cutree_cols = 2,
                           color = colorRampPalette(c("#4292C6", "#F0F0F0", "#A50F15"))(100),
                           scale = "row", border_color = "black", angle_col = 90,
                           show_rownames = TRUE, show_colnames = FALSE, na.rm = TRUE,
                           clustering_method = "ward.D2",
                           cluster_cols = TRUE, cluster_rows = TRUE)

# Sample clustering and ordering
cluster_order <- cutree(heatmap_output$tree_col, k = 2)
ordered_samples <- t(core_matrix[, heatmap_output$tree_col$order])
write.csv(ordered_samples, "../03.out/figure1/01.FPMTB_exVivo_Resistant_Sensitive_HeatmapOrder.csv")

# Assign cluster group (Resistant or Sensitive)
cluster_assignments <- data.frame(cutree(heatmap_output$tree_col, k = 2))
colnames(cluster_assignments) <- "Cluster"
data <- cbind(t(core_matrix), cluster_assignments)
data$group <- ifelse(data$Cluster == 1, "Resistant", "Sensitive")
table(data$group)
write.csv(data, "../03.out/figure1/02.FPMTB_exVivo_Resistant_Sensitive_Raw.csv")

# Heatmap with annotation
cluster_assignments$Group <- ifelse(cluster_assignments$Cluster == 1, "Resistant", "Sensitive")
annotation_col <- data.frame(Group = factor(cluster_assignments$Group))
rownames(annotation_col) <- colnames(core_matrix)

ann_colors <- list(Group = c(Sensitive = "#7c9d97", Resistant = "#e9b383"))

pdf("../03.out/figure1/01.FPMTB_exVivo_drug_heatmap.pdf", width = 10, height = 6)
pheatmap(core_matrix,
         cutree_cols = 3,
         color = colorRampPalette(c("#4292C6", "#F0F0F0", "#A50F15"))(100),
         scale = "row", border_color = "black", angle_col = 90,
         show_rownames = TRUE, show_colnames = FALSE, na.rm = TRUE,
         clustering_method = "ward.D",
         annotation_colors = ann_colors,
         annotation_col = annotation_col,
         cluster_cols = TRUE, cluster_rows = TRUE)
dev.off()

# Clustering samples into three groups
cluster_assignments <- data.frame(cutree(heatmap_output$tree_col, k = 3))
colnames(cluster_assignments) <- "Cluster"
data <- cbind(t(core_matrix), cluster_assignments)
data$group <- ifelse(data$Cluster == 1, "Resistant", "Sensitive")
table(data$Cluster)
write.csv(data, "../03.out/figure1/03.FPMTB_exVivo_Resistant_Sensitive_Raw.csv")

##########################################################################
### Boxplot for Drug Sensitivity Between Resistant and Sensitive Groups
##########################################################################

library(ggplot2)
library(ggpubr)
library(ggrastr)
library(tidyverse)

# Resistant vs. Sensitive groups
drug_sensitivity <- read.csv("../03.out/figure1/01.FPMTB_exVivo_Resistant_Sensitive_Raw.csv", header = TRUE, row.names = 1)
drug_sensitivity <- drug_sensitivity[, c("Azacitidine", "Venetoclax", "group")]
drug_sensitivity <- gather(drug_sensitivity, -group, key = "drug", value = "sDSS")

p <- ggplot(drug_sensitivity, aes(x = group, y = sDSS, color = group)) + 
  geom_boxplot_jitter() + geom_quasirandom(width = 0.15, alpha = 0.5) +
  facet_wrap(~drug) + 
  theme_pubr() + 
  xlab("") + ylab("Drug sensitivity (sDSS)") +
  theme(strip.text.x = element_text(size = 13)) +
  stat_compare_means(comparisons = list(c('Resistant', 'Sensitive')), method = 'wilcox.test', label = "p.signif") + 
  scale_color_manual(values = c("#e9b383", "#7c9d97")) +  
  theme(legend.position = 'None')

p
ggsave(p, filename = "../03.out/figure1/02.FPMTB_Resistant_Sensitive_DrugsDSS_group.pdf", width = 6, height = 4)

# Clusters (1, 2, 3)
drug_sensitivity_cluster <- read.csv("../03.out/figure1/01.FPMTB_exVivo_Resistant_Sensitive_Raw.csv", header = TRUE, row.names = 1)
drug_sensitivity_cluster <- drug_sensitivity_cluster[, c("Azacitidine", "Venetoclax", "Cluster")]
drug_sensitivity_cluster <- gather(drug_sensitivity_cluster, -Cluster, key = "drug", value = "sDSS")

drug_sensitivity_cluster$Cluster <- as.factor(drug_sensitivity_cluster$Cluster)

p1 <- ggplot(drug_sensitivity_cluster, aes(x = Cluster, y = sDSS, color = Cluster)) + 
  geom_boxplot_jitter() + geom_quasirandom(width = 0.15, alpha = 0.5) +
  facet_wrap(~drug) + 
  theme_pubr() + 
  xlab("") + ylab("Drug sensitivity (sDSS)") +
  theme(strip.text.x = element_text(size = 13)) +
  stat_compare_means(comparisons = list(c('1', '2'), c('1', '3'), c('2', '3')), method = 'wilcox.test', label = "p.signif") + 
  scale_color_manual(values = c("#e9b383", "#7c9d97", "#9cb0c3")) +  
  theme(legend.position = 'None')

p1
ggsave(p1, filename = "../03.out/figure1/02.FPMTB_Resistant_Sensitive_DrugsDSS_Cluster.pdf", width = 6, height = 4)

##########################################################################
### Differential Expression Analysis and Volcano Plot
##########################################################################

library(DESeq2)
source("RunDESeq2.R")
sample_order <- read.csv("../03.out/figure1/01.FPMTB_exVivo_Resistant_Sensitive_Raw.csv", header = TRUE, row.names = 1)
sample_order <- sample_order[order(sample_order$group), ]
table(sample_order$group)

tpm_data <- read.csv("../02.data/File_8_RNA_seq_TPM_163S_4Healthy.csv", header = TRUE)
raw_counts <- read.csv("../02.data/File_8_RNA_seq_Raw_Reads_163S_4Healthy.csv", header = TRUE)
merged_data <- merge(tpm_data[, 1:3], raw_counts, by = "X")
merged_data <- subset(merged_data, Type == "protein_coding")
filtered_tpm <- tpm_data[rowMeans(tpm_data[, -c(1:3)]) > 0.1, ]
row.names(merged_data) <- merged_data$X
merged_data <- merged_data[, -c(1:3)]

# Intersect samples
overlap_samples <- intersect(row.names(sample_order), colnames(merged_data))
sample_order <- sample_order[overlap_samples, ]
filtered_counts <- merged_data[, as.character(row.names(sample_order))]
filtered_tpm <- filtered_tpm[, c("ID", "Symbol", "Type", as.character(row.names(sample_order)))]

RunDESeq2(count_mat = filtered_counts, n.cont = 43, n.treat = 55,
          prefix = "../03.out/figure1/FPMTB_DESeq2_out_coding_filteredLow", sort.p = FALSE,
          merge.normalized = TRUE, normalized_mat = filtered_tpm)

# Volcano plot for Differential Expression Analysis
source("plotVolcanoV3.R")
de_results <- read.table("../03.out/figure1/FPMTB_DESeq2_out_coding_filteredLow_with_normalized_mat.txt", header = TRUE)
de_genes <- subset(de_results, abs(log2FoldChange) > 0.5 & padj < 0.05)

de_genes_up <- subset(de_results, log2FoldChange > 0.5 & padj < 0.05)
de_genes_down <- subset(de_results, log2FoldChange < -0.5 & padj < 0.05)

de_results$log2FoldChange[de_results$log2FoldChange < -4] <- -4
de_results$padj[de_results$padj < 1e-10] <- 1e-10

p1 <- plotVolcano(mat = de_results, gene.col = "Symbol", x.col = "log2FoldChange", y.col = "padj", 
                  labx = "log2FoldChange", laby = "-Log10FDR",
                  x_cut1 = 0.5, x_cut2 = 1, y_cut1 = 0.05, y_cut2 = 0.01, x.lim = 4, y.lim = 10, 
                  label = FALSE, title = "AZA+VEN (Ex vivo)", selected_genes = "NA") + theme_classic()
p1
ggsave(p1, filename = "../03.out/figure1/03.FPMTB_volcano_plot.pdf", width = 4.5, height = 4)


##########################################################################
### Generate Heatmap for RNA-seq Data
##########################################################################

# Load required package
library(pheatmap)

# Load data
heatmap_data <- read.csv("figure1/FPMTB_heatmap.csv", header = TRUE, row.names = 1)

# Prepare the matrix for heatmap generation
core_matrix <- as.matrix(heatmap_data)
breaks_seq <- unique(c(seq(-3, 3, length = 100)))

# Define sample group annotations
sample_groups <- data.frame(
  Group = factor(c(rep("Resistant", 43), rep("Sensitive", 55)))
)
rownames(sample_groups) <- colnames(core_matrix)

# Define colors for annotation groups
annotation_colors <- list(
  Group = c(Sensitive = "#7c9d97", Resistant = "#e9b383")
)

# Create and save the heatmap
pdf("figure1/Out_FPMTB_heatmap.pdf", width = 5.5, height = 6)
pheatmap(core_matrix,
         breaks = breaks_seq,
         cutree_cols = 2,
         cutree_rows = 2,
         color = colorRampPalette(c("#4292C6", "#F0F0F0", "#A50F15"))(100),
         scale = "row", 
         border_color = "black", 
         angle_col = 90,
         show_rownames = TRUE, 
         show_colnames = FALSE, 
         na.rm = TRUE,
         clustering_method = "ward.D",
         annotation_colors = annotation_colors,
         annotation_col = sample_groups,
         cluster_cols = FALSE, 
         cluster_rows = TRUE)
dev.off()


##########################################################################
### Gene Set Enrichment Analysis (GSEA)
##########################################################################

library(clusterProfiler)
library(enrichplot)
library(GseaVis)

DEgene_data <- read.csv("figure1/FPMTB_DESeq2_out.csv", header = TRUE)
DE_list <- dplyr::distinct(DEgene_data[, c(9, 3)], ID, .keep_all = TRUE) %>%
  na.omit() %>%
  mutate(rank = rank(Corr, ties.method = "random"))

geneList <- sort(DE_list$Corr, decreasing = TRUE)
names(geneList) <- DE_list$ID

gene_sets <- read.gmt("LSC_MONO_genesets.gmt")
set.seed(123456)
gsea_result <- GSEA(geneList, TERM2GENE = gene_sets, pvalueCutoff = 1, maxGSSize = 800, 
                    pAdjustMethod = "BH", seed = TRUE)

pdf("figure1/Out_GSEA_LSC_UP.pdf", width = 4.5, height = 3.5)
gseaplot(gsea_result, by = "runningScore", title = gsea_result$Description[1], geneSetID = 1)
dev.off()

pdf("figure1/Out_GSEA.LSC-core17.pdf", width = 4.5, height = 3.5)
gseaplot(gsea.enrich, by = "runningScore", 
         title = gsea.enrich$Description[4], geneSetID = 4)
dev.off()

pdf("figure1/Out_GSEA.Eppert et al. LSC Up.pdf", width = 4.5, height = 3.5)
gseaplot(gsea.enrich, by = "runningScore", 
         title = gsea.enrich$Description[6], geneSetID = 6)
dev.off()

pdf("figure1/Out_GSEA.MONOCYTE.pdf", width = 4.5, height = 3.5)
gseaplot(gsea.enrich, by = "runningScore", 
         title = gsea.enrich$Description[5], geneSetID = 5)
dev.off()

########################################
### CIBERSORTx Results Visualization
########################################

# Load required libraries
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(tidyverse)
library(ggrastr)
library(ggbeeswarm)

# Load CIBERSORTx results and filter necessary columns
cibersortx_results <- read.csv("figure1/CIBERSORTx_Results_FPMTB.csv", header = TRUE, row.names = 1, check.names = FALSE)
cibersortx_results <- cibersortx_results[, c("LSPC-Primed", "Mono-like", "group")] 
# Available columns for future use: "LSPC-Quiescent", "LSPC-Cycle", "GMP-like", "ProMono-like", "cDC-like"

# Transform data from wide to long format
cibersortx_long <- cibersortx_results %>% 
  gather(-group, key = "cell_type", value = "Score")

# Generate boxplot for cellular abundance comparison (Resistant vs. Sensitive)
p <- ggplot(cibersortx_long, aes(x = group, y = Score, color = group)) + 
  geom_boxplot_jitter() + 
  geom_quasirandom(width = 0.15, alpha = 0.5) +
  facet_wrap(~cell_type, ncol = 1) + 
  theme_pubr() + 
  xlab("") + 
  ylab("Cellular abundance (CIBERSORTx)") +
  theme(strip.text.x = element_text(size = 13)) +
  stat_compare_means(comparisons = list(c('Resistant', 'Sensitive')), method = 'wilcox.test', label = "p.signif") + 
  scale_color_manual(values = c("#e9b383", "#7c9d97")) +  
  theme(legend.position = 'None')

# Display plot
p

# Save plot to file
ggsave(p, filename = "figure1/Out_FPMTB_Resistant_Sensitive_cibersortx.pdf", width = 2.5, height = 4.5)


#####################################
## GO Enrichment Dot Plot
#####################################

# Load GO enrichment data
go_enrichment_data <- read.csv("figure1/Metascape_DEGs_plot.csv", header = TRUE)

# Reorder factor levels based on GO description
go_enrichment_data$Description <- factor(go_enrichment_data$Description, levels = rev(go_enrichment_data$Description))

# Generate dot plot for GO enrichment results
p <- ggplot(go_enrichment_data, aes(x = NA, y = Description, color = Z.score, size = negLogP)) + 
  geom_point(stroke = 1) + 
  labs(x = "", y = "") +
  theme_classic() +
  scale_color_gradient(low = "#DADAEB", high = "#7582c1") +
  scale_size("negLogP", limits = c(5, 18), range = c(2, 6)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45))

# Display plot
p

# Save plot to file
ggsave(p, filename = "figure1/Out_DEGs_enrichment_dotplot.pdf", width = 4, height = 5)