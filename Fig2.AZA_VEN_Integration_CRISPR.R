####################################
### Table 1 for RJAML-cohort 1
####################################

library(readxl)
library(tableone)

# Load data and filter for relevant responses
data <- read_xlsx("01.110_AZA_VEN_RNAseq_ClinicalMut.xlsx", sheet = 1)
# data$Gender <- ifelse(data$Gender == 1, "Male", "Female")
# data$Secondary_AML <- ifelse(data$Secondary_AML == 1, "Yes", "No")
# data$Complex_Karyotype <- ifelse(data$Complex_Karyotype == 1, "Yes", "No")
table(data$Best_response)

# Retain mutations mutated in >= 3 patients
variables <- c("Age", "BM_blast", "Best_response", "Gender", "Secondary_AML", "Complex_Karyotype", "ELN")
grouping <- "Response"

# Define categorical variables
factorVars <- c("Best_response", "Gender", "Secondary_AML", "Complex_Karyotype", "ELN")

# Create Table 1 object
table1 <- CreateTableOne(vars = variables, strata = grouping, data = data, factorVars = factorVars, addOverall = TRUE)
table1
write.csv(print(table1), "figure2/rjaml-cohort1_baseline.csv")

# library(table1)
# table1(~ Age + Gender + BM_blast + Secondary_AML + Cpmplex_Karyotype + ELN + `FLT3-ITD` + NPM1|Response, data = data)
# # Print the summary of Table 1 with non-normal option for continuous variables
# print(table1, showAllLevels = TRUE, nonnormal = c("Age", "BM_blast"))



#####################################################
### Kaplan-Meier (KM) Plots for RJAML Cohort 1
#####################################################

library(ggplot2)
library(survminer)
library(survival)

# Load data
data <- read_xlsx("01.RF8_SurvivalAnalysis.xlsx")
group_counts <- table(data$Response)
legend_labs <- c(paste0("Non-CR/CRi (N = ", group_counts["0"], ")"),
                 paste0("CR/CRi (N = ", group_counts["1"], ")"))

# Event-Free Survival (EFS) - Calculate C-index and HR
fit_efs <- survfit(Surv(EFS_Time, EFS_Status) ~ Response, data = data)
cox_model_efs <- coxph(Surv(EFS_Time, EFS_Status) ~ Response, data = data)
c_index_efs <- summary(cox_model_efs)$concordance[1]
hr_efs <- exp(coef(cox_model_efs))[1]
hr_ci_efs <- exp(confint(cox_model_efs))[1, ]
hr_label_efs <- paste0("HR: ", round(hr_efs, 2), " (95% CI: ", round(hr_ci_efs[1], 2), " - ", round(hr_ci_efs[2], 2), ")")
c_index_label_efs <- paste0("C-index: ", round(c_index_efs, 3))

# Generate KM plot for EFS
p1 <- ggsurvplot(fit_efs, data = data, pval = TRUE, conf.int = FALSE,
                 risk.table = TRUE, font.legend = 13,
                 legend.title = "Response",
                 legend = c(0.75, 0.85),
                 xlab = 'Survival time (years)', 
                 ylab = 'Survival probability (EFS)',
                 legend.labs = legend_labs) + 
  labs(caption = paste0(c_index_label_efs, "\n", hr_label_efs))

# Save EFS KM plot
pdf("figure2/01.RF8_EFS.pdf", width = 5, height = 5.5, onefile = FALSE)
p1
dev.off()

# Overall Survival (OS) - Calculate C-index and HR
fit_os <- survfit(Surv(OS_Time, OS_Status) ~ Response, data = data)
cox_model_os <- coxph(Surv(OS_Time, OS_Status) ~ Response, data = data)
c_index_os <- summary(cox_model_os)$concordance[1]
hr_os <- exp(coef(cox_model_os))[1]
hr_ci_os <- exp(confint(cox_model_os))[1, ]
hr_label_os <- paste0("HR: ", round(hr_os, 2), " (95% CI: ", round(hr_ci_os[1], 2), " - ", round(hr_ci_os[2], 2), ")")
c_index_label_os <- paste0("C-index: ", round(c_index_os, 3))

# Generate KM plot for OS
p2 <- ggsurvplot(fit_os, data = data, pval = TRUE, conf.int = FALSE,
                 risk.table = TRUE, font.legend = 13,
                 legend.title = "Response",
                 legend = c(0.75, 0.85),
                 xlab = 'Survival time (years)', 
                 ylab = 'Survival probability (OS)',
                 legend.labs = legend_labs) + 
  labs(caption = paste0(c_index_label_os, "\n", hr_label_os))

# Save OS KM plot
pdf("figure2/01.RF8_OS.pdf", width = 5, height = 5.5, onefile = FALSE)
p2
dev.off()

##########################################################################
### Differential Expression Analysis for RJAML AZA+VEN Cohort (In Vivo)
##########################################################################

library(dplyr)
library(DESeq2)
source("RunDESeq2.R")

# Load TPM and raw count data
tpm <- read.table("../02.data/rjaml_TPM_matrix.txt", header = TRUE, check.names = FALSE)
mat <- read.table("../02.data/rjaml_count_matrix.txt", header = TRUE, row.names = 1, check.names = FALSE)
sample_info <- read.csv("../02.data/rjaml_rnaseq_group.csv", header = TRUE, check.names = FALSE)

# Filter data for relevant samples
tpm <- tpm[, c("ID", "Symbol", "Type", as.character(sample_info$SampleID))]
mat <- mat[, as.character(sample_info$SampleID)]
tpm <- tpm[which(rowMeans(tpm[, -c(1:3)]) > 0.5), ]

# Run DESeq2 analysis
RunDESeq2(count_mat = mat, n.cont = 33, n.treat = 77,
          prefix = "../03.out/figure2/rjaml_DESeq2_out", sort.p = FALSE,
          merge.normalized = TRUE,
          normalized_mat = tpm)

# Volcano plot for Differential Expression Analysis (In vivo)
source("plotVolcanoV3.R")
de_results <- read.table("../03.out/figure2/rjaml_DESeq2_out_with_normalized_mat.txt", header = TRUE)
de_results_filtered <- de_results[, 1:9]
degs_in_vivo <- subset(de_results_filtered, abs(log2FoldChange) > 0.5 & padj < 0.05)

# Cap extreme values
de_results_filtered$log2FoldChange[de_results_filtered$log2FoldChange < -4] <- -4
de_results_filtered$log2FoldChange[de_results_filtered$log2FoldChange > 4] <- 4
de_results_filtered$padj[de_results_filtered$padj < 1e-8] <- 1e-8

# Generate volcano plot
volcano_plot <- plotVolcano(mat = de_results_filtered, gene.col = "Symbol", x.col = "log2FoldChange", y.col = "padj", 
                            labx = "log2FoldChange", laby = "-Log10FDR",
                            x_cut1 = 0.5, x_cut2 = 1, y_cut1 = 0.05, y_cut2 = 0.01, x.lim = 4, y.lim = 8, 
                            label = FALSE, title = "AZA+VEN (In vivo)", selected_genes = "NA") + theme_classic()

# Save volcano plot
ggsave(volcano_plot, filename = "../03.out/figure2/04.RJAML_valcanoplot.pdf", width = 4.5, height = 4)


####################################
### Venn Diagram for DEGs
####################################

library(ggvenn)
library(patchwork)

# Function to create Venn diagram and perform hypergeometric test
create_venn_plot <- function(gene_list1, gene_list2, labels, colors, N) {
  # Remove NAs and duplicates
  gene_list1 <- na.omit(unique(gene_list1))
  gene_list2 <- na.omit(unique(gene_list2))
  
  # Check for non-empty lists
  if (length(gene_list1) == 0 || length(gene_list2) == 0) {
    warning("One or both gene lists are empty. Cannot create Venn diagram.")
    return(NULL)
  }
  
  # Create Venn plot
  venn_plot <- ggvenn(
    list(gene_list1, gene_list2),
    fill_color = colors,
    stroke_size = 0.5,
    set_name_size = 4,
    auto_scale = FALSE
  )
  
  # Calculate overlap and perform hypergeometric test
  overlap <- length(intersect(gene_list1, gene_list2))
  p_value <- phyper(overlap - 1, length(gene_list2), N - length(gene_list2), length(gene_list1), lower.tail = FALSE)
  
  # Annotate plot with p-value
  venn_plot <- venn_plot + annotate("text", x = 0, y = 1, label = paste0("P-value = ", signif(p_value, 3)))
  
  return(venn_plot)
}

# Define total number of genes in the background
N <- 20000

# Generate Venn plots for DEGs, up-regulated, and down-regulated genes
p0 <- create_venn_plot(degs1$Symbol, degs_in_vivo$Symbol, c("AZA+VEN (Ex vivo)", "AZA+VEN (In vivo)"), c("#868686FF", "#CD534CFF"), N)
p1 <- create_venn_plot(degs_up$Symbol, degs_up2$Symbol, c("AZA+VEN (Ex vivo) UP", "AZA+VEN (In vivo) UP"), c("#0073C2FF", "#EFC000FF"), N)
p2 <- create_venn_plot(degs_down$Symbol, degs_down2$Symbol, c("AZA+VEN (Ex vivo) DN", "AZA+VEN (In vivo) DN"), c("#868686FF", "#CD534CFF"), N)

# Combine plots
combined_plots <- p0 | p1 | p2
combined_plots

# Save combined Venn plots
ggsave(combined_plots, filename = "../03.out/figure2/05.vennPlots.pdf", width = 9, height = 3)


#################################################
### Integration of DEGs (Ex vivo and In vivo)
#################################################

library(ggplot2)
library(ggrepel)

# Identify overlapping genes between ex vivo and in vivo DEGs
overlap_genes <- intersect(degs1$Symbol, degs2$Symbol)

# Extract relevant columns for ex vivo and in vivo DEGs
degs_exvivo <- subset(m2, Symbol %in% overlap_genes)[, c("ID", "Symbol", "log2FoldChange", "padj")]
degs_invivo <- subset(m4, Symbol %in% overlap_genes)[, c("log2FoldChange", "padj")]

# Combine the two DEG datasets
degs <- cbind(degs_exvivo, degs_invivo)
colnames(degs) <- c("ID", "Symbol", "ExVivo_log2FC", "ExVivo_padj", "InVivo_log2FC", "InVivo_padj")

# Calculate -log10(padj) for ex vivo and in vivo
degs$ExVivo_neg_log10_padj <- -log10(degs$ExVivo_padj)
degs$InVivo_neg_log10_padj <- -log10(degs$InVivo_padj)

# Define quadrants based on log2 fold changes
degs$quadrant <- with(degs, ifelse(ExVivo_log2FC > 0.5 & InVivo_log2FC > 0.5, "x > 0.5 & y > 0.5",
                                   ifelse(ExVivo_log2FC < -0.5 & InVivo_log2FC < -0.5, "x < -0.5 & y < -0.5",
                                          ifelse(ExVivo_log2FC > 0.5 & InVivo_log2FC < -0.5, "x > 0.5 & y < -0.5", 
                                                 "x < -0.5 & y > 0.5"))))

# Label selected genes
genes_to_label <- c("BCL2L1", "PINK1", "CEBPE", "CSMD3", "HPSE2", "CD177", "ROBO2", "TEAD2", "CLEC2L", "TPSAB1", 
                    "USP44", "SLC9A3", "CEBPA", "MANEAL")
degs$label <- ifelse(degs$Symbol %in% genes_to_label, degs$Symbol, NA)

# Create scatter plot for DEGs
p <- ggplot(degs, aes(x = ExVivo_log2FC, y = InVivo_log2FC)) +
  geom_point(aes(size = ExVivo_neg_log10_padj + InVivo_neg_log10_padj, color = quadrant), alpha = 0.6) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "black") +
  geom_hline(yintercept = c(-0.5, 0.5), linetype = "dashed", color = "black") +
  scale_color_manual(values = c("x > 0.5 & y > 0.5" = "#006CAF", 
                                "x < -0.5 & y < -0.5" = "#D078A7",
                                "x > 0.5 & y < -0.5" = "#A09D9B", 
                                "x < -0.5 & y > 0.5" = "#00976A")) +
  geom_label_repel(aes(label = label), size = 4, box.padding = 0.35, point.padding = 0.5, 
                   segment.color = 'grey50', arrow = arrow(length = unit(0.02, "npc")), 
                   segment.size = 0.5, label.size = 0.5, color = "black", fill = "white") +
  labs(x = "Log2 FoldChange \n(Ex vivo, CR/CRi vs. Non-CR/CRi)", 
       y = "Log2 FoldChange \n(In vivo, CR/CRi vs. Non-CR/CRi)", 
       size = "-log10(padj)", 
       color = "Quadrant") +
  theme_bw()

# Save the plot
ggsave(p, filename = "../03.out/figure2/06.inteDotplot.pdf", width = 5.8, height = 4)


#################################################
### GSEA Analysis for CRISPRi Drug Sensitivity
#################################################

library(tidyverse)
library(dplyr)
library(GseaVis)
library(ggpubr)
library(patchwork)
library(clusterProfiler)
library(enrichplot)

# Load and prepare the CRISPRi dataset
DE_list <- read.csv("../02.data/CRISPR_Ven.csv", header = TRUE)
DE_list <- DE_list[, c(1, 3)]
colnames(DE_list) <- c("ID", "Corr")
DE_list <- dplyr::distinct(DE_list, ID, .keep_all = TRUE)
DE_list <- na.omit(DE_list)

# Rank genes based on correlation values
DE_list <- DE_list %>% mutate(rank = rank(Corr, ties.method = "random"))
geneList <- DE_list$Corr
names(geneList) <- DE_list$ID
geneList <- sort(geneList, decreasing = TRUE)

# Create gene sets for GSEA
geneSet <- data.frame(
  term = c(rep("AZA_Ven_Sensitive_DN", length(intersect(degs_down$Symbol, degs_down2$Symbol))),
           rep("AZA_Ven_Sensitive_UP", length(intersect(degs_up$Symbol, degs_up2$Symbol)))),
  gene = c(intersect(degs_down$Symbol, degs_down2$Symbol),
           intersect(degs_up$Symbol, degs_up2$Symbol))
)

# Check the gene sets created
table(geneSet$term)

# Perform GSEA
set.seed(123456)
gsea.enrich <- GSEA(geneList, TERM2GENE = geneSet, pvalueCutoff = 1, nPermSimple = 1000,
                    maxGSSize = 1100, pAdjustMethod = "BH", seed = TRUE, eps = 0)

# Save GSEA results
pdf("../03.out/figure2/07.CRISPRi_gsea_AZA.pdf", width = 6, height = 6.5)
gseaplot(gsea.enrich, pvalue_table = TRUE, geneSetID = c("AZA_Ven_Sensitive_DN"))
dev.off()


##########################################################################
### Differential Expression Analysis for shPINK1
##########################################################################

# library(GenomicFeatures)
# library(stringr)
# dir <- "D:\\01.Projects\\01.project\\35_AZA_VEN\\00.RNAseq_Model\\figure2"
# setwd(dir)
# files <- list.files(dir, pattern = ".count", recursive=TRUE)
# count.files <- list()
# 
# for(i in 1:length(files)){
#   count.files[[i]] <- read.table(files[[i]], header=FALSE, sep="\t", row.names=1)
# }
# count.matrix <- do.call(cbind, count.files)
# colnames <- gsub(".count", "", files)
# colnames(count.matrix) <- colnames
# count.matrix <- count.matrix[1:60483, ]
# ## Save the count table
# write.table(count.matrix, "shPINK1_count.matrix.txt", sep="\t", row.names=T, col.names=NA, quote=FALSE)
# 
# gtffile <- "D:\\01.Projects\\01.project\\gencode.v22.annotation.gtf.gz"
# txdb <- makeTxDbFromGFF(gtffile, format="gtf")
# ebg <- exonsBy(txdb, by="gene")
# ebgList <- sum(width(reduce(ebg)))
# genes <- intersect(rownames(count.matrix), names(ebgList))
# Length <- as.vector(ebgList[genes])
# Length <- as.vector(Length)
# ## Normalize
# TPM <- t(t(count.matrix / t(Length)) * 1e6 / colSums(count.matrix / t(Length)))
# TPM <- data.frame(TPM)
# TPM$ID <- row.names(TPM)
# ## gene anno
# anno <- read.table("geneAnnotation.txt", header = F)
# colnames(anno) <- c("ID", "Symbol", "Type")
# tpm <- merge(anno, TPM, by = "ID")
# write.table(tpm, "shPINK1_ID_matched_TPM_matrix.txt")

####
#### PCA analysis
####
library(ggplot2)
library(dplyr)
library(ggfortify)
library(plotly)
library(scatterplot3d)

tpm <- read.table("figure2\\shPINK1_ID_matched_TPM_matrix.txt", header = TRUE, check.names = FALSE)
## MOLM13 cell lines
dat <- tpm[, 4:9]
row.names(dat) <- tpm$ID
dat <- dat[which(rowMeans(dat[, -c(1:3)]) > 0), ]
scaled_matrix <- scale(t(dat))

pca_result <- prcomp(scaled_matrix, scale. = TRUE)
# Create a grouping vector
groups <- c(rep("Scramble", 2), rep("shPINK1-1", 2), rep("shPINK1-2", 2))

pca_2d_df <- as.data.frame(pca_result$x)
pca_2d_df$Group <- groups

# Define colors for groups
group_colors <- c("red", "green", "blue")
names(group_colors) <- unique(groups)

# Add variance explained for axis labels
explained_variance <- summary(pca_result)$importance[2, ]
x_label <- paste0("PC1 (", round(explained_variance[1] * 100, 2), "%)")
y_label <- paste0("PC2 (", round(explained_variance[2] * 100, 2), "%)")

pdf("figure2/shPINK1_M13_2D_PCA.pdf", width = 5.3, height = 4)

ggplot(pca_2d_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3) +
  scale_color_manual(values = group_colors) +
  labs(
    x = paste0("PC1 (", round(explained_variance[1] * 100, 2), "%)"),
    y = paste0("PC2 (", round(explained_variance[2] * 100, 2), "%)"),
    title = "MOLM13 cells"
  ) +
  theme_classic() +
  theme(legend.title = element_blank()) +
  guides(color = guide_legend(title = "Groups"))
dev.off()


## OA3 cell lines
dat <- tpm[, 10:15]
row.names(dat) <- tpm$ID
dat <- dat[which(rowMeans(dat[, -c(1:3)]) > 0), ]
scaled_matrix <- scale(t(dat))

pca_result <- prcomp(scaled_matrix, scale. = TRUE)
# Create a grouping vector
groups <- c(rep("Scramble", 2), rep("shPINK1-1", 2), rep("shPINK1-2", 2))

pca_2d_df <- as.data.frame(pca_result$x)
pca_2d_df$Group <- groups

# Define colors for groups
group_colors <- c("red", "green", "blue")
names(group_colors) <- unique(groups)

# Add variance explained for axis labels
explained_variance <- summary(pca_result)$importance[2, ]
x_label <- paste0("PC1 (", round(explained_variance[1] * 100, 2), "%)")
y_label <- paste0("PC2 (", round(explained_variance[2] * 100, 2), "%)")

pdf("figure2/shPINK1_OA3_2D_PCA.pdf", width = 5.3, height = 4)
ggplot(pca_2d_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3) +
  scale_color_manual(values = group_colors) +
  labs(
    x = paste0("PC1 (", round(explained_variance[1] * 100, 2), "%)"),
    y = paste0("PC2 (", round(explained_variance[2] * 100, 2), "%)"),
    title = "OA3 cells"
  ) +
  theme_classic() +
  theme(legend.title = element_blank()) +
  guides(color = guide_legend(title = "Groups"))
dev.off()


library(dplyr)
library(DESeq2)

source("D:\\03.R_functions\\JPComUse\\R\\RunDESeq2.R")
# Load TPM and raw count data
tpm <- read.table("figure2\\shPINK1_ID_matched_TPM_matrix.txt", header = TRUE, check.names = FALSE)
mat <- read.table("figure2\\shPINK1_count.matrix.txt", header = TRUE, row.names = 1, check.names = FALSE)
tpm <- subset(tpm, Type == "protein_coding")
## MOLM13 sh1 vs. control
samples <- c("M0_rep1", "M0_rep2", "M1_rep1", "M1_rep2")
# Filter low-expressed genes
tpm1 <- tpm[which(rowMeans(tpm[, -c(1:3)]) > 0), c("ID", "Symbol", "Type", samples)]
mat1 <- mat[tpm1$ID, samples]
# Run DESeq2 analysis
RunDESeq2(count_mat = mat1, n.cont = 2, n.treat = 2,
          prefix = "figure2\\M13_sh1_control_out", sort.p = FALSE,
          merge.normalized = TRUE,
          normalized_mat = tpm1)

## MOLM13 sh2 vs. control
samples <- c("M0_rep1", "M0_rep2", "M2_rep1", "M2_rep2")
# Filter low-expressed genes
tpm1 <- tpm[which(rowMeans(tpm[, -c(1:3)]) > 0), c("ID", "Symbol", "Type", samples)]
mat1 <- mat[tpm1$ID, samples]
# Run DESeq2 analysis
RunDESeq2(count_mat = mat1, n.cont = 2, n.treat = 2,
          prefix = "figure2\\M13_sh2_control_out", sort.p = FALSE,
          merge.normalized = TRUE,
          normalized_mat = tpm1)

## MOLM13 sh1 sh2 vs. control
samples <- c("M0_rep1", "M0_rep2", "M1_rep1", "M1_rep2", "M2_rep1", "M2_rep2")
# Filter low-expressed genes
tpm1 <- tpm[which(rowMeans(tpm[, -c(1:3)]) > 0), c("ID", "Symbol", "Type", samples)]
mat1 <- mat[tpm1$ID, samples]
# Run DESeq2 analysis
RunDESeq2(count_mat = mat1, n.cont = 2, n.treat = 4,
          prefix = "figure2\\M13_sh12_control_out", sort.p = FALSE,
          merge.normalized = TRUE,
          normalized_mat = tpm1)


## OA3 sh1 vs. control
samples <- c("O0_rep1", "O0_rep2", "O1_rep1", "O1_rep2")
# Filter low-expressed genes
tpm1 <- tpm[which(rowMeans(tpm[, -c(1:3)]) > 0), c("ID", "Symbol", "Type", samples)]
mat1 <- mat[tpm1$ID, samples]
# Run DESeq2 analysis
RunDESeq2(count_mat = mat1, n.cont = 2, n.treat = 2,
          prefix = "figure2\\OA3_sh1_control_out", sort.p = FALSE,
          merge.normalized = TRUE,
          normalized_mat = tpm1)

## OA3 sh4 vs. control
samples <- c("O0_rep1", "O0_rep2", "O4_rep1", "O4_rep2")
# Filter low-expressed genes
tpm1 <- tpm[which(rowMeans(tpm[, -c(1:3)]) > 0), c("ID", "Symbol", "Type", samples)]
mat1 <- mat[tpm1$ID, samples]
# Run DESeq2 analysis
RunDESeq2(count_mat = mat1, n.cont = 2, n.treat = 2,
          prefix = "figure2\\OA3_sh4_control_out", sort.p = FALSE,
          merge.normalized = TRUE,
          normalized_mat = tpm1)

## OA3 sh1 sh4 vs. control
samples <- c("O0_rep1", "O0_rep2", "O1_rep1", "O1_rep2", "O4_rep1", "O4_rep2")
# Filter low-expressed genes
tpm1 <- tpm[which(rowMeans(tpm[, -c(1:3)]) > 0), c("ID", "Symbol", "Type", samples)]
mat1 <- mat[tpm1$ID, samples]
# Run DESeq2 analysis
RunDESeq2(count_mat = mat1, n.cont = 2, n.treat = 4,
          prefix = "figure2\\OA3_sh14_control_out", sort.p = FALSE,
          merge.normalized = TRUE,
          normalized_mat = tpm1)

####################################
######### Valcano plots
####################################
library(patchwork)
source("D:\\03.R_functions\\JPComUse\\R\\plotVolcano.R")
resis_sen <- read.csv("figure2\\selected_genes_resistant_sensitive.csv", header = T)  
# gsym  Pathway
# Gene1 KEGG_RIBOSOME
# Gene2 KEGG_DNA_REPLICATION
custom_colors <- brewer.pal(4, "BrBG")  # PiYG palette

## MOLM13
## shRNA1 vs Scramble (control)
deseq_output1 <- read.table("figure2//M13_sh1_control_out_with_normalized_mat.txt", header = T)
selected_data1 <- deseq_output1[, 1:9]
# Cap log2 fold changes at 5 to avoid extreme values affecting visualization
selected_data1$log2FoldChange[selected_data1$log2FoldChange > 5] <- 5
# Set a minimum threshold for adjusted p-values to avoid extremely small values
selected_data1$padj[selected_data1$padj < 1e-60] <- 1e-60
degs <- subset(selected_data1, padj < 0.05)
overlap_genes <- intersect(resis_sen$gsym, degs$Symbol)
selected_genes <- subset(resis_sen, gsym %in% overlap_genes)

# Generate the volcano plot
volcano_plot1 <- plotVolcano(
  data = selected_data1, gene.col = "Symbol",
  logFC.col = "log2FoldChange", pval.col = "padj",
  selected_genes = selected_genes,
  logFCcut = 0.5, logFCcut2 = 1, logFCcut3 = 2,
  pvalCut = 0.05, pvalCut2 = 0.0001, pvalCut3 = 0.00001,
  threshold_colors = custom_colors,
  plot_mode = "advanced", x_range = c(-5, 5), y_range = c(0, 60),
  title = "MOLM13: shPINK1-1 vs. Scramble"
)

## shRNA2 vs Scramble (control)
deseq_output2 <- read.table("figure2//M13_sh2_control_out_with_normalized_mat.txt", header = T)
selected_data2 <- deseq_output2[, 1:9]
selected_data2$log2FoldChange[selected_data2$log2FoldChange > 5] <- 5
selected_data2$log2FoldChange[selected_data2$log2FoldChange < -5] <- -5
selected_data2$padj[selected_data2$padj < 1e-60] <- 1e-60
degs <- subset(selected_data2, padj < 0.05)
overlap_genes <- intersect(resis_sen$gsym, degs$Symbol)
selected_genes <- subset(resis_sen, gsym %in% overlap_genes)

volcano_plot2 <- plotVolcano(
  data = selected_data2, gene.col = "Symbol",
  logFC.col = "log2FoldChange", pval.col = "padj",
  selected_genes = selected_genes,
  logFCcut = 0.5, logFCcut2 = 1, logFCcut3 = 2,
  pvalCut = 0.05, pvalCut2 = 0.0001, pvalCut3 = 0.00001,
  threshold_colors = custom_colors,
  plot_mode = "advanced", x_range = c(-5, 5), y_range = c(0, 60),
  title = "MOLM13: shPINK1-2 vs. Scramble"
)

## shRNA1-2 vs Scramble (control)
deseq_output3 <- read.table("figure2//M13_sh12_control_out_with_normalized_mat.txt", header = T)
selected_data3 <- deseq_output3[, 1:9]
selected_data3$log2FoldChange[selected_data3$log2FoldChange > 5] <- 5
selected_data3$log2FoldChange[selected_data3$log2FoldChange < -5] <- -5
selected_data3$padj[selected_data3$padj < 1e-50] <- 1e-50
degs <- subset(selected_data3, padj < 0.05)
overlap_genes <- intersect(resis_sen$gsym, degs$Symbol)
selected_genes <- subset(resis_sen, gsym %in% overlap_genes)

volcano_plot3 <- plotVolcano(
  data = selected_data3, gene.col = "Symbol",
  logFC.col = "log2FoldChange", pval.col = "padj",
  selected_genes = selected_genes,
  logFCcut = 0.5, logFCcut2 = 1, logFCcut3 = 2,
  pvalCut = 0.05, pvalCut2 = 0.0001, pvalCut3 = 0.00001,
  threshold_colors = custom_colors,
  plot_mode = "advanced", x_range = c(-5, 5), y_range = c(0, 50),
  title = "MOLM13: shPINK1-1-2 vs. Scramble"
)

plots <- volcano_plot1 | volcano_plot2 | volcano_plot3
ggsave("figure2\\MOLM13_shPINK1_volcano_plot.pdf", plot = plots, width = 17, height = 5)


## OA3
## shRNA1 vs Scramble (control)
# gsym  Pathway
# Gene1 KEGG_RIBOSOME
# Gene2 KEGG_DNA_REPLICATION
custom_colors <- brewer.pal(4, "PiYG")  # PiYG BrBG palette

deseq_output1 <- read.table("figure2//OA3_sh1_control_out_with_normalized_mat.txt", header = T)
selected_data1 <- deseq_output1[, 1:9]
# Cap log2 fold changes at 5 to avoid extreme values affecting visualization
selected_data1$log2FoldChange[selected_data1$log2FoldChange > 5] <- 5
# Set a minimum threshold for adjusted p-values to avoid extremely small values
selected_data1$padj[selected_data1$padj < 1e-30] <- 1e-30
degs <- subset(selected_data1, padj < 0.05)
overlap_genes <- intersect(resis_sen$gsym, degs$Symbol)
selected_genes <- subset(resis_sen, gsym %in% overlap_genes)

# Generate the volcano plot
volcano_plot1 <- plotVolcano(
  data = selected_data1, gene.col = "Symbol",
  logFC.col = "log2FoldChange", pval.col = "padj",
  selected_genes = selected_genes,
  logFCcut = 0.5, logFCcut2 = 1, logFCcut3 = 2,
  pvalCut = 0.05, pvalCut2 = 0.0001, pvalCut3 = 0.00001,
  threshold_colors = custom_colors,
  plot_mode = "advanced", x_range = c(-5, 5), y_range = c(0, 30),
  title = "OA3: shPINK1-1 vs. Scramble"
)

## shRNA2 vs Scramble (control)
deseq_output2 <- read.table("figure2//OA3_sh4_control_out_with_normalized_mat.txt", header = T)
selected_data2 <- deseq_output2[, 1:9]
selected_data2$log2FoldChange[selected_data2$log2FoldChange > 5] <- 5
selected_data2$log2FoldChange[selected_data2$log2FoldChange < -5] <- -5
selected_data2$padj[selected_data2$padj < 1e-60] <- 1e-60
degs <- subset(selected_data2, padj < 0.05)
overlap_genes <- intersect(resis_sen$gsym, degs$Symbol)
selected_genes <- subset(resis_sen, gsym %in% overlap_genes)

volcano_plot2 <- plotVolcano(
  data = selected_data2, gene.col = "Symbol",
  logFC.col = "log2FoldChange", pval.col = "padj",
  selected_genes = selected_genes,
  logFCcut = 0.5, logFCcut2 = 1, logFCcut3 = 2,
  pvalCut = 0.05, pvalCut2 = 0.0001, pvalCut3 = 0.00001,
  threshold_colors = custom_colors,
  plot_mode = "advanced", x_range = c(-5, 5), y_range = c(0, 60),
  title = "MOLM13: shPINK1-2 vs. Scramble"
)

## shRNA1-2 vs Scramble (control)
deseq_output3 <- read.table("figure2//M13_sh12_control_out_with_normalized_mat.txt", header = T)
selected_data3 <- deseq_output3[, 1:9]
selected_data3$log2FoldChange[selected_data3$log2FoldChange > 5] <- 5
selected_data3$log2FoldChange[selected_data3$log2FoldChange < -5] <- -5
selected_data3$padj[selected_data3$padj < 1e-40] <- 1e-40
degs <- subset(selected_data3, padj < 0.05)
overlap_genes <- intersect(resis_sen$gsym, degs$Symbol)
selected_genes <- subset(resis_sen, gsym %in% overlap_genes)

volcano_plot3 <- plotVolcano(
  data = selected_data3, gene.col = "Symbol",
  logFC.col = "log2FoldChange", pval.col = "padj",
  selected_genes = selected_genes,
  logFCcut = 0.5, logFCcut2 = 1, logFCcut3 = 2,
  pvalCut = 0.05, pvalCut2 = 0.0001, pvalCut3 = 0.00001,
  threshold_colors = custom_colors,
  plot_mode = "advanced", x_range = c(-5, 5), y_range = c(0, 40),
  title = "MOLM13: shPINK1-1-2 vs. Scramble"
)

plots <- volcano_plot1 | volcano_plot2 | volcano_plot3
ggsave("figure2\\OA3_shPINK1_volcano_plot.pdf", plot = plots, width = 17, height = 5)

####################################
######### GSEA analysis
####################################

library(ggplot2)
library(clusterProfiler)
library(enrichplot)
library(tidyverse)
library(dplyr)
library(GseaVis)

DEgene_mat <- read.table("figure2//M13_sh12_control_out_with_normalized_mat.txt", header = T)
DE_list <- DEgene_mat[,c(9,3)]
colnames(DE_list) <- c("ID","Corr")
DE_list <- dplyr::distinct(DE_list, ID, .keep_all=TRUE)
DE_list <- na.omit(DE_list)

DE_list <- DE_list %>%
  mutate(rank = rank(Corr,  ties.method = "random")) 

geneList <- DE_list$Corr
names(geneList)=DE_list$ID
geneList <- sort(geneList,decreasing = T)
# geneSet <- read.gmt("..//GOBP_OXIDATIVE_PHOSPHORYLATION.v2024.1.Hs.gmt") #   c2.cp.kegg.v2023.1.Hs.symbols.gmt
geneSet <- read.gmt("figure2//AZA_VEN_responsive_genes.gmt")
set.seed(123456)

gsea.enrich <- GSEA(geneList,TERM2GENE = geneSet, pvalueCutoff = 1, nPermSimple = 1000, 
                    maxGSSize = 800, pAdjustMethod = "BH", seed = TRUE)
gsea.out <- gsea.enrich@result


###### 
pdf("figure2//GSEA_OA3_sh12.pdf", 6, 3.3)
gseaNb(object = gsea.enrich,
       geneSetID = c("Resistant_genes","Sensitive_genes"), 
       addPval = T,
       subPlot = 2,
       curveCol= c('#7582c1', '#f99f1c'), 
       htCol= c( "#7582c1", "#f99f1c"), 
       rankCol= c( "#7582c1", "white", "#f99f1c"))
dev.off()


#####################################
##   Dot plot for GO enrichment    ##
#####################################

library(ggplot2)
gsea.out <- read.csv("figure2\\gsea_dotplots.csv", header = T)
gsea.out$Ajusted.Pvalue <- ifelse(gsea.out$Ajusted.Pvalue == "<0.001", 0.001, gsea.out$Ajusted.Pvalue)
gsea.out$Ajusted.Pvalue <- as.numeric(gsea.out$Ajusted.Pvalue)
gsea.out$negLog10P.adj <- -log10(gsea.out$Ajusted.Pvalue)
gsea.out$p_label <- ifelse(gsea.out$Ajusted.Pvalue < 0.05, "*", NA)
gsea.out$p_color <- ifelse(gsea.out$Ajusted.Pvalue < 0.05 & gsea.out$NES > 0, "#40004B", 
                      ifelse(gsea.out$P.value < 0.05 & gsea.out$NES < 0, "#1B7837", NA))
# gsea.out$geneset <- factor(gsea.out$geneset, levels = c("Sensitive", "Resistant"))
gsea.out$sh <- factor(gsea.out$sh, levels = rev(c("shPINK1-1", "shPINK1-2", "shPINK1-1-2")))

p <- ggplot(gsea.out, aes(x=geneset, y=sh, fill=NES)) +
  theme_classic() +
  facet_grid(.~Cell_line) +
  geom_point(shape = 21, aes(size = 0.5), color="white") +
  geom_text(aes(label=p_label, color = p_color), size=5, vjust=0.76) +
  scale_fill_gradient2(low = "#007575", high = "#FFAA00", mid="#F7F7F7", midpoint=0) +
  scale_color_manual(values=c("#D3D3D3", "#000000")) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    strip.text.y = element_text(angle = 0),
    strip.text.x = element_text(angle = 0),
    strip.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA)
  ) +
  guides(size="none", color="none") +
  labs(x="", y="",fill="NES")
p

ggsave(p, filename = "figure2\\GSEA_dotplot.pdf", width = 3.5, height = 3.2)





















