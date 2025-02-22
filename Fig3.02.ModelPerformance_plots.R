####################################
### Model Performance Plot
####################################

# Set working directory
setwd("D:\\01.Projects\\01.project\\35_AZA_VEN\\00.RNAseq_Model")

# Load necessary libraries
library(ggplot2)
library(tidyverse)
library(RColorBrewer)

# Load the results
data <- read.csv("figure3/ML_Results_GeneExp.csv", header = TRUE, row.names = 1, check.names = FALSE)

# Melt the data for ggplot
melted_data <- data %>%
  pivot_longer(cols = c(classif.auc, classif.acc), 
               names_to = "metric", 
               values_to = "value") %>%
  mutate(learner_id = toupper(gsub("classif.", "", learner_id)),
         metric = toupper(gsub("classif.", "", metric)))

# Filter out non-finite values
melted_data <- melted_data %>% filter(is.finite(value))

# Get the median AUC for sorting
median_auc <- melted_data %>%
  filter(metric == "AUC") %>%
  group_by(learner_id) %>%
  summarise(median_value = mean(value, na.rm = TRUE)) %>%
  arrange(desc(median_value)) %>%
  pull(learner_id)

# Update factor levels
melted_data$learner_id <- factor(melted_data$learner_id, levels = median_auc)
melted_data$metric <- factor(melted_data$metric, levels = c("AUC", "ACC"))

# Define color palette
palette_colors <- colorRampPalette(rev(brewer.pal(6, "Blues")))(length(unique(melted_data$learner_id)))

# Create boxplot
p <- ggplot(melted_data, aes(x = learner_id, y = value, fill = learner_id)) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "text", aes(label = round(after_stat(y), 3)), 
               vjust = -0.5, color = "black", size = 3.5) + 
  facet_wrap(.~metric, ncol = 1) +
  scale_fill_manual(values = palette_colors) +
  labs(title = "", x = "", y = "Model Performance") +
  theme_classic() +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, hjust = 1))

# Save plot
ggsave(p, filename = "figure3/01.Model_performances.pdf", width = 6, height = 5)

####################################
### Feature Importance Plot
####################################

# Load feature importance data
dat <- read.csv("figure3/RF_FeatureImportance.csv", header = TRUE)

# Create feature importance bar plot
p1 <- ggplot(dat, aes(x = reorder(Feature, Importance), y = Importance, fill = Importance)) +
  geom_bar(stat = "identity", width = 0.7) +
  coord_flip() +
  labs(title = "", x = "", y = "Feature Importance") +
  scale_fill_gradient(low = "#C6DBEF", high = "#2171B5") +
  theme_bw()

# Save plot
ggsave(p1, filename = "figure3/02.Features_Importance.pdf", width = 3.6, height = 3.4)

####################################
## Heatmap for Confusion Matrices
####################################

# Load test set data
data <- read.csv("figure3/RF8_RJAML_cohort1_test.csv", header = TRUE, row.names = 1)
data$Predicted_Response <- ifelse(data$RF8.prob.CR > 0.5, "Responder", "Non-responder")

# Create confusion matrix for test set
confusion_matrix <- table(data$Response, data$Predicted_Response)

# Convert to data frame
confusion_matrix_df <- as.data.frame(confusion_matrix)
colnames(confusion_matrix_df) <- c("True_Response", "Predicted_Response", "Count")

# Calculate percentages
total_counts <- rowSums(confusion_matrix)
confusion_matrix_df$Percent <- with(confusion_matrix_df, Count / total_counts[True_Response])

# Plot confusion matrix heatmap for the test set
p2 <- ggplot(confusion_matrix_df, aes(x = Predicted_Response, y = True_Response, fill = Percent)) +
  geom_tile(color = "white") +
  geom_text(aes(label = paste(Count, "\n(", sprintf("%.2f", Percent), ")", sep = "")), size = 4) +
  scale_fill_gradient(low = "#F7FBFF", high = "#2171B5", limits = c(0, 1)) +
  labs(title = "Test Set", x = "Predicted Response", y = "True Response") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Save plot
ggsave(p2, filename = "figure3/03.ConfusionMat_Test.pdf", width = 4.7, height = 3.8)

# Load entire set data
data <- read.csv("figure3/RF8_RJAML_cohort1_entire.csv", header = TRUE, row.names = 1)
data$Predicted_Response <- ifelse(data$RF8.prob.CR > 0.5, "Responder", "Non-responder")

# Create confusion matrix for entire set
confusion_matrix <- table(data$Response, data$Predicted_Response)

# Convert to data frame
confusion_matrix_df <- as.data.frame(confusion_matrix)
colnames(confusion_matrix_df) <- c("True_Response", "Predicted_Response", "Count")

# Calculate percentages
total_counts <- rowSums(confusion_matrix)
confusion_matrix_df$Percent <- with(confusion_matrix_df, Count / total_counts[True_Response])

# Plot confusion matrix heatmap for the entire set
p3 <- ggplot(confusion_matrix_df, aes(x = Predicted_Response, y = True_Response, fill = Percent)) +
  geom_tile(color = "white") +
  geom_text(aes(label = paste(Count, "\n(", sprintf("%.2f", Percent), ")", sep = "")), size = 4) +
  scale_fill_gradient(low = "#F7FBFF", high = "#2171B5", limits = c(0, 1)) +
  labs(title = "Entire Set", x = "Predicted Response", y = "True Response") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Save plot
ggsave(p3, filename = "figure3/04.ConfusionMat_Entire.pdf", width = 4.7, height = 3.8)

####################################
## Density Plot for Response Probability
####################################

# Load entire set data again
data <- read.csv("figure3/RF8_RJAML_cohort1_entire.csv", header = TRUE, row.names = 1)

# Create density plot for response probability
p <- ggdensity(data, x = "RF8.prob.CR", add = "mean", rug = TRUE, color = "Response",
               xlab = "Response Probability (RF8)", ylab = "Density",
               fill = "Response", palette = "npg")

# Save plot
ggsave(p, filename = "figure3/RF8_RJAML_cohort1_entire_density.pdf", width = 4.5, height = 4.3)


######################################################
### Violin and Boxplots for Response Probability (RF8)
######################################################

# Load necessary libraries
library(ggplot2)
library(reshape2)
library(ggpubr)
library(patchwork)

# Define a function to generate plots for each dataset
generate_plot <- function(file_path, dataset_label) {
  data <- read.csv(file_path, header = TRUE, row.names = 1)
  data$dataset <- dataset_label
  
  p <- ggplot(data, aes(x = Response, y = RF8.prob.CR, fill = Response)) +
    geom_violin(trim = FALSE) +
    geom_boxplot(width = 0.1, position = position_dodge(0.9), outlier.shape = NA) +
    scale_fill_manual(values = c("#93C89A", "#FFCC98")) +
    facet_wrap(.~ dataset) +
    labs(title = "",
         x = "True Response",
         y = "Response probability (RF8)") +
    theme_classic() +
    theme(legend.position = "none") +
    stat_compare_means(aes(label = ..p.format..), comparisons = list(c("CR", "NonCR")))
  
  return(p)
}

# Generate plots for each dataset
p4 <- generate_plot("figure3/RF8_RJAML_cohort1_test.csv", "Test set")
p5 <- generate_plot("figure3/RF8_RJAML_cohort1_entire.csv", "Entire set")
p6 <- generate_plot("figure3/RF8_FPMTB.csv", "FPMTB set")
p7 <- generate_plot("figure3/RF8_BeatAML.csv", "BeatAML set")
p8 <- generate_plot("figure3/RF8_Cohort2.csv", "Cohort 2")

# Combine all plots into a single layout using patchwork
plots <- p4 | p5 | p6 | p7 | p8

# Save the combined plot as a PDF file
ggsave(plots, filename = "figure3/RF8_violinplot.pdf", width = 10, height = 3)


####################################
## Barplot for Response
####################################
# Load necessary libraries
library(ggstatsplot)
library(ggplot2)
library(patchwork)

# Test set bar plot
data_test <- read.csv("figure3/RF8_RJAML_cohort1_test.csv", header = TRUE, row.names = 1)
data_test$Predicted_Response <- ifelse(data_test$RF8.prob.CR > 0.5, "Responder", "Non-responder")
p_test <- ggbarstats(data = data_test, x = Response, y = Predicted_Response) + 
  theme_classic() +
  ggtitle("Test Set")

# Entire set bar plot
data_entire <- read.csv("figure3/RF8_RJAML_cohort1_entire.csv", header = TRUE, row.names = 1)
data_entire$Predicted_Response <- ifelse(data_entire$RF8.prob.CR > 0.5, "Responder", "Non-responder")
p_entire <- ggbarstats(data = data_entire, x = Response, y = Predicted_Response) + 
  theme_classic() +
  ggtitle("Entire Set")

# Combine plots
combined_barplot <- p_test | p_entire
ggsave(combined_barplot, filename = "figure3/RF8_barplot_Test_Entire.pdf", width = 7, height = 4.5)


###################################################
## Integration of Ex Vivo Drug-Response in BeatAML
###################################################

library(reshape2)
library(ggridges)
library(ggplot2)
library(tidyr)
library(tidyverse)
library(ggpubr)
library(patchwork)

# Read and filter drug AUC data for Azacitidine and Venetoclax
drug_auc <- read.table("../02.data/figure3/01_BeatAML_Inhibitor_AUC_values_v4_dbgap.txt", header = TRUE, sep = "\t", check.names = FALSE)
drug_select <- subset(drug_auc, inhibitor %in% c("Azacytidine", "Venetoclax"))
drug_select <- drug_select[, c("Sample", "inhibitor", "auc")]

# Reshape data for heatmap
dat <- dcast(drug_select, Sample ~ inhibitor)
colnames(dat) <- c("Sample", "Azacitidine", "Venetoclax")
id <- read.csv("../02.data/figure3/02_RNAseq_Sample_ids.csv", header = TRUE, check.names = FALSE)
dat <- merge(dat, id, by = "Sample")
row.names(dat) <- dat$Sample

# Filter samples based on specific criteria
dat <- subset(dat, diseaseStageAtSpecimenCollection %in% c("Initial Diagnosis", "Relapse"))
dat <- subset(dat, rnaSeq == "y")
dat$NA_number <- rowSums(is.na(dat[, 2:3]))
dat <- subset(dat, NA_number == 0)
dat <- dat[, c("Azacitidine", "Venetoclax")]
head(dat)

####################################
## Heatmap for Drug Resistance
####################################
library(pheatmap)
core_mat <- dat
core_mat1 <- log2(t(as.matrix(core_mat)) + 1)

# Customize heatmap breaks and color palette
bk <- unique(c(seq(-2, 2, length = 100)))
out <- pheatmap(core_mat1,
                breaks = bk, 
                cutree_cols = 2,
                color = colorRampPalette(c("#8C510A", "#F5F5F5", "#01665E"))(100), 
                scale = "row", border_color = "black", angle_col = 90,
                show_rownames = TRUE, show_colnames = FALSE,
                clustering_method = "ward.D",
                cluster_cols = TRUE, cluster_rows = TRUE)

# Process clustered sample data
clust_sample <- colnames(core_mat1[, out$tree_col[["order"]]])
clust_sample <- data.frame(clust_sample)
sample_cluster <- data.frame(cutree(out$tree_col, k = 3))
colnames(sample_cluster) <- "Cluster"

# Assign labels for resistant and sensitive groups
data <- data.frame(core_mat)
data <- cbind(data, sample_cluster)
data$group <- ifelse(data$Cluster == 1, "Resistant", "Sensitive")
write.csv(data, "../03.out/figure3/BeatAML_exVivo_Resistant_Sensitive_Raw.csv")

# Annotate heatmap
sample_cluster$Label <- ifelse(sample_cluster$Cluster == 1, "Resistant", "Sensitive")
annotation_col <- data.frame(Label = factor(sample_cluster$Label))
rownames(annotation_col) <- colnames(core_mat1)

ann_colors <- list(
  Label = c(Sensitive = "#7c9d97", Resistant = "#e9b383")
)

# Save heatmap to
pdf("../03.out/figure3/01.BeatAML_exVivo_drug_heatmap.pdf", width = 10, height = 2)
pheatmap(core_mat1,
         breaks = bk, 
         cutree_cols = 3,
         color = colorRampPalette(c("#8C510A", "#F5F5F5", "#01665E"))(100), 
         scale = "row", border_color = "black", angle_col = 90,
         show_rownames = TRUE, show_colnames = FALSE,
         clustering_method = "ward.D",
         annotation_colors = ann_colors,
         annotation_col = annotation_col,
         cluster_cols = TRUE, cluster_rows = TRUE)
dev.off()


###############################################################
## Boxplot for Drug Sensitivity (Azacitidine & Venetoclax)
###############################################################

library(ggplot2)
library(ggpubr)
library(tidyr)

# Read drug sensitivity data
mat <- read.csv("../03.out/figure3/BeatAML_exVivo_Resistant_Sensitive_Raw.csv", header = TRUE, row.names = 1)
mat <- mat[, c("group", "Azacitidine", "Venetoclax")]
mat <- gather(mat, key = "drug", value = "AUC", -group)

# Boxplot for Azacitidine and Venetoclax
p_box <- ggplot(mat, aes(x = group, y = AUC, color = group)) +
  geom_boxplot_jitter() + 
  geom_quasirandom(width = 0.15, alpha = 0.75) +
  facet_wrap(.~drug, ncol = 2) + 
  theme_pubr() + 
  xlab("") + ylab("Drug Resistance (AUC)") +
  theme(strip.text.x = element_text(size = 13)) +
  stat_compare_means(label = "p.signif", comparisons = list(c('Resistant', 'Sensitive')), method = 'wilcox.test') + 
  scale_color_manual(values = c("#e9b383", "#7c9d97")) + 
  theme(legend.position = "none")

# Save boxplot to file
ggsave(p_box, filename = "../03.out/figure3/02.BeatAML_Drug_Sensitivity_Boxplot.pdf", width = 8, height = 4)

###############################################################
## Boxplot for Drug Sensitivity Across Clusters
###############################################################
# Prepare data for cluster-based analysis
mat_cluster <- read.csv("../03.out/figure3/BeatAML_exVivo_Resistant_Sensitive_Raw.csv", header = TRUE, row.names = 1)
mat_cluster <- gather(mat_cluster[, c("Azacitidine", "Venetoclax", "Cluster")], key = "drug", value = "sDSS", -Cluster)
mat_cluster$Cluster <- as.factor(mat_cluster$Cluster)

# Boxplot for cluster-based drug sensitivity
p_cluster <- ggplot(mat_cluster, aes(x = Cluster, y = sDSS, color = Cluster)) + 
  geom_boxplot_jitter() + 
  geom_quasirandom(width = 0.15, alpha = 0.5) +
  facet_wrap(.~drug) + 
  theme_pubr() + 
  xlab("") + ylab("Drug Sensitivity (AUC)") +
  theme(strip.text.x = element_text(size = 13)) +
  stat_compare_means(comparisons = list(c('1', '2'), c('1', '3'), c('2', '3')), 
                     method = 'wilcox.test', label = "p.signif") + 
  scale_color_manual(values = c("#e9b383", "#7c9d97", "#9cb0c3")) + 
  theme(legend.position = "none")

# Save cluster-based boxplot to file
ggsave(p_cluster, filename = "../03.out/figure3/03.BeatAML_Cluster_Drug_Sensitivity_Boxplot.pdf", width = 8, height = 4)
