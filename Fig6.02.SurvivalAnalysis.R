# Set working directory
setwd("D:/01.Projects/01.project/35_AZA_VEN/00.RNAseq_Model")

####################################
# Monotonic Response Prediction Plot
####################################

# Load necessary libraries
library(ggplot2)
library(scales)

# Load data for plotting
df <- read.csv("figure6/source_data_forplot.csv", header = TRUE)

# Create the plot
p1 <- ggplot(df, aes(x = RF8_score, y = Prob_mean * 100)) + 
  geom_line(color = "red") + 
  geom_ribbon(aes(ymin = Prob_lower * 100, ymax = Prob_upper * 100), fill = "red", alpha = 0.2) +  
  labs(x = "RF8 Score", y = "Response Probability (%)") + 
  theme_classic() + 
  scale_y_continuous(breaks = seq(0, 100, by = 10)) +  
  scale_x_continuous(breaks = seq(0, 1, by = 0.1)) + 
  theme(
    panel.grid.major = element_line(size = 0.25, linetype = "dashed", colour = "grey"),
    panel.grid.minor = element_line(size = 0.25, linetype = "dashed", colour = "grey"),
    panel.grid.major.x = element_blank()  
  )

# Save the plot 
ggsave(p1, filename = "figure6/01.RF8_MonotonicPlot.pdf", width = 5, height = 4)

####################################
# Monotonic Survival Prediction
####################################

library(survival)
library(survminer)
library(readxl)

# Load survival analysis data
data <- read_xlsx("figure6/01.RF8_SurvivalAnalysis.xlsx")
data$RF8.group <- factor(data$RF8.group, levels = c("High", "Intermediate", "Low"))

# Count samples in each group
n_high <- sum(data$RF8.group == "High")
n_intermediate <- sum(data$RF8.group == "Intermediate")
n_low <- sum(data$RF8.group == "Low")

# Fit the survival model (Event-Free Survival)
fit_efs <- survfit(Surv(EFS_Time, EFS_Status) ~ RF8.group, data = data)

# Perform Cox proportional hazards model for pairwise comparisons
cox_model_efs <- coxph(Surv(EFS_Time, EFS_Status) ~ RF8.group, data = data)
summary_cox_efs <- summary(cox_model_efs)
hr_values_efs <- exp(summary_cox_efs$coefficients[, "coef"])  # Hazard ratios
p_values_efs <- summary_cox_efs$coefficients[, "Pr(>|z|)"]    # P-values

# Extract HR and p-values for group comparisons
hr_high_intermediate_efs <- hr_values_efs["RF8.groupIntermediate"]
pval_high_intermediate_efs <- p_values_efs["RF8.groupIntermediate"]

hr_high_low_efs <- hr_values_efs["RF8.groupLow"]
pval_high_low_efs <- p_values_efs["RF8.groupLow"]

# Create labels for hazard ratios and p-values
label_high_intermediate_efs <- paste0("HR (Intermediate vs High): ", round(hr_high_intermediate_efs, 2), 
                                      ", P = ", format.pval(pval_high_intermediate_efs, digits = 3))
label_high_low_efs <- paste0("HR (Low vs High): ", round(hr_high_low_efs, 2), 
                             ", P = ", format.pval(pval_high_low_efs, digits = 3))

# Create legend labels with sample counts
legend_labels_efs <- c(
  paste0("High (n = ", n_high, ")"),
  paste0("Intermediate (n = ", n_intermediate, ")"),
  paste0("Low (n = ", n_low, ")")
)

# Plot Kaplan-Meier curves for Event-Free Survival
p2 <- ggsurvplot(fit_efs, data = data, pval = TRUE, conf.int = FALSE,
                 palette = c("#7c9d97", "#9cb0c3", "#e9b383"),
                 risk.table = FALSE, font.legend = 13,
                 legend.title = "RF8",
                 legend = c(0.75, 0.85),
                 xlab = 'Survival Time (years)', 
                 ylab = 'Survival Probability (EFS)',
                 legend.labs = legend_labels_efs)

# Annotate the plot with HR and p-values
p2$plot <- p2$plot +
  annotate("text", x = max(data$OS_Time) * 0.8, y = 0.4, label = label_high_intermediate_efs, size = 4, hjust = 0) +
  annotate("text", x = max(data$OS_Time) * 0.8, y = 0.3, label = label_high_low_efs, size = 4, hjust = 0)

# Save the EFS plot as a PDF
pdf("figure6/02.RF8_EFS.pdf", width = 5, height = 4.6, onefile = FALSE)
p2
dev.off()

####################################
# Overall Survival Prediction (OS)
####################################

# Fit the survival model (Overall Survival)
fit_os <- survfit(Surv(OS_Time, OS_Status) ~ RF8.group, data = data)

# Perform Cox proportional hazards model for pairwise comparisons
cox_model_os <- coxph(Surv(OS_Time, OS_Status) ~ RF8.group, data = data)
summary_cox_os <- summary(cox_model_os)
hr_values_os <- exp(summary_cox_os$coefficients[, "coef"])
p_values_os <- summary_cox_os$coefficients[, "Pr(>|z|)"]

# Extract HR and p-values for group comparisons
hr_high_intermediate_os <- hr_values_os["RF8.groupIntermediate"]
pval_high_intermediate_os <- p_values_os["RF8.groupIntermediate"]

hr_high_low_os <- hr_values_os["RF8.groupLow"]
pval_high_low_os <- p_values_os["RF8.groupLow"]

# Create labels for hazard ratios and p-values
label_high_intermediate_os <- paste0("HR (Intermediate vs High): ", round(hr_high_intermediate_os, 2), 
                                     ", P = ", format.pval(pval_high_intermediate_os, digits = 3))
label_high_low_os <- paste0("HR (Low vs High): ", round(hr_high_low_os, 2), 
                            ", P = ", format.pval(pval_high_low_os, digits = 3))

# Plot Kaplan-Meier curves for Overall Survival
p3 <- ggsurvplot(fit_os, data = data, pval = TRUE, conf.int = FALSE,
                 palette = c("#7c9d97", "#9cb0c3", "#e9b383"),
                 risk.table = FALSE, font.legend = 13,
                 legend.title = "RF8",
                 legend = c(0.75, 0.85),
                 xlab = 'Survival Time (years)', 
                 ylab = 'Survival Probability (OS)',
                 legend.labs = legend_labels_efs)

# Annotate the plot with HR and p-values
p3$plot <- p3$plot +
  annotate("text", x = max(data$OS_Time) * 0.8, y = 0.4, label = label_high_intermediate_os, size = 4, hjust = 0) +
  annotate("text", x = max(data$OS_Time) * 0.8, y = 0.3, label = label_high_low_os, size = 4, hjust = 0)

# Save the OS plot as a PDF
pdf("figure6/02.RF8_OS.pdf", width = 5, height = 4.6, onefile = FALSE)
p3
dev.off()

####################################
# Brier Score Analysis
####################################

# Load necessary libraries for Brier score analysis
library(pec)
library(survival)
library(readxl)

# Load data
data <- read_xlsx("figure6/01.RF8_SurvivalAnalysis.xlsx")

# Define models for Event-Free Survival (EFS)
models_efs <- list(
  "RF8" = coxph(Surv(EFS_Time, EFS_Status) ~ RF8.prob.CR, data = data, x = TRUE, y = TRUE),
  "ELN2022" = coxph(Surv(EFS_Time, EFS_Status) ~ ELN2022, data = data, x = TRUE, y = TRUE),
  "gene4_signature" = coxph(Surv(EFS_Time, EFS_Status) ~ gene4_signature, data = data, x = TRUE, y = TRUE),
  "molecular_predictor" = coxph(Surv(EFS_Time, EFS_Status) ~ molecular_predictor, data = data, x = TRUE, y = TRUE),
  "mayo_predictor" = coxph(Surv(EFS_Time, EFS_Status) ~ mayo_predictor, data = data, x = TRUE, y = TRUE)
)

# Calculate Brier scores for EFS
brier_efs <- pec(models_efs, data = data, formula = Surv(EFS_Time, EFS_Status) ~ 1)
print(brier_efs)

# Plot Brier scores for EFS
pdf("figure6/01.BrierScores_EFS.pdf", width = 4, height = 4.3)
plot(
  brier_efs,
  xlim = c(0, 3),
  col = c("black", "#E74342", "#3469AF", "#38AC6B", "#9A78B4", "#73C8EF"),
  lwd = 2.5,
  xlab = "Survival Time (years)",
  ylab = "Prediction Error (Brier Score)",
  main = "Brier Score for Event-Free Survival (EFS)",
  legend = TRUE
)
grid(col = "gray", lty = "dotted", lwd = 0.8)
dev.off()

# Define models for Overall Survival (OS)
models_os <- list(
  "RF8" = coxph(Surv(OS_Time, OS_Status) ~ RF8.prob.CR, data = data, x = TRUE, y = TRUE),
  "ELN2022" = coxph(Surv(OS_Time, OS_Status) ~ ELN2022, data = data, x = TRUE, y = TRUE),
  "gene4_signature" = coxph(Surv(OS_Time, OS_Status) ~ gene4_signature, data = data, x = TRUE, y = TRUE),
  "molecular_predictor" = coxph(Surv(OS_Time, OS_Status) ~ molecular_predictor, data = data, x = TRUE, y = TRUE),
  "mayo_predictor" = coxph(Surv(OS_Time, OS_Status) ~ mayo_predictor, data = data, x = TRUE, y = TRUE)
)

# Calculate Brier scores for OS
brier_os <- pec(models_os, data = data, formula = Surv(OS_Time, OS_Status) ~ 1)
print(brier_os)

# Plot Brier scores for OS
pdf("figure6/01.BrierScores_OS.pdf", width = 4, height = 4.3)
plot(
  brier_os,
  xlim = c(0, 3),
  col = c("black", "#E74342", "#3469AF", "#38AC6B", "#9A78B4", "#73C8EF"),
  lwd = 2.5,
  xlab = "Survival Time (years)",
  ylab = "Prediction Error (Brier Score)",
  main = "Brier Score for Overall Survival (OS)",
  legend = TRUE
)
grid(col = "gray", lty = "dotted", lwd = 0.8)
dev.off()


####################################
### C-index Comparison
####################################

# Load required libraries
library(survcomp)
library(survival)
library(readxl)
library(ggplot2)

# Load data
data <- read_xlsx("figure6/01.RF8_SurvivalAnalysis.xlsx")
features <- c("RF8.prob.CR", "ELN2022", "gene4_signature", "molecular_predictor", "mayo_predictor")

# Initialize an empty data frame to store C-index results
c_index_df <- data.frame(
  Feature = character(),
  Outcome = character(),
  C_index = numeric(),
  Lower_CI = numeric(),
  Upper_CI = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each feature to calculate C-index and 95% CI for both OS and EFS
for (feature in features) {
  
  # Calculate C-index for OS
  cox_model_os <- coxph(Surv(OS_Time, OS_Status) ~ data[[feature]], data = data)
  predictions_os <- predict(cox_model_os)
  c_index_result_os <- concordance.index(predictions_os, surv.time = data$OS_Time, 
                                         surv.event = data$OS_Status)
  
  # Store OS results in the dataframe
  c_index_df <- rbind(c_index_df, data.frame(
    Feature = feature,
    Outcome = "OS",
    C_index = round(c_index_result_os$c.index, 3),
    Lower_CI = round(c_index_result_os$lower, 3),
    Upper_CI = round(c_index_result_os$upper, 3)
  ))
  
  # Calculate C-index for EFS
  cox_model_efs <- coxph(Surv(EFS_Time, EFS_Status) ~ data[[feature]], data = data)
  predictions_efs <- predict(cox_model_efs)
  c_index_result_efs <- concordance.index(predictions_efs, surv.time = data$EFS_Time, 
                                          surv.event = data$EFS_Status)
  
  # Store EFS results in the dataframe
  c_index_df <- rbind(c_index_df, data.frame(
    Feature = feature,
    Outcome = "EFS",
    C_index = round(c_index_result_efs$c.index, 3),
    Lower_CI = round(c_index_result_efs$lower, 3),
    Upper_CI = round(c_index_result_efs$upper, 3)
  ))
}

# Display the C-index dataframe
c_index_df

# Set factor levels for Feature to maintain consistent ordering in the plot
c_index_df$Feature <- factor(c_index_df$Feature, levels = features)

# Plot C-index results
p4 <- ggplot(c_index_df, aes(x = Feature, y = C_index, fill = Feature)) +
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.2, size = 0.8) +
  geom_bar(stat = "identity", width = 1) +
  facet_wrap(~ Outcome) +  
  labs(y = "C-index", x = "", title = "C-index with 95% CI for OS and EFS") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Set1") 
p4

# Save the C-index plot
ggsave(p4, filename = "figure6/03.cindex.pdf", width = 7.4, height = 5)

####################################
### KM-Plots for Cohort 2
####################################

library(ggplot2)
library(survminer)
library(survival)

# Load survival data for Cohort 2
data <- read.csv("figure6/RF8_Cohort2_survival.csv", header = TRUE)

# Create predicted groups based on RF8 probability
data$Predicted_group <- ifelse(data$RF8.prob.CR > 0.9, "0.9-1",
                               ifelse(data$RF8.prob.CR > 0.6, "0.6-0.9",
                                      ifelse(data$RF8.prob.CR > 0.3, "0.3-0.6", "0.1-0.3")))
group_counts <- table(data$Predicted_group)

# Create labels for the legend with group sizes
legend_labs <- c(paste0("0.1-0.3 (N = ", group_counts["0.1-0.3"], ")"),
                 paste0("0.3-0.6 (N = ", group_counts["0.3-0.6"], ")"),
                 paste0("0.6-0.9 (N = ", group_counts["0.6-0.9"], ")"),
                 paste0("0.9-1 (N = ", group_counts["0.9-1"], ")"))

# EFS - Calculate C-index and HR
fit_efs <- survfit(Surv(EFS_Time, EFS_Status) ~ Predicted_group, data = data)
cox_model_efs <- coxph(Surv(EFS_Time, EFS_Status) ~ Predicted_group, data = data) 
c_index_efs <- summary(cox_model_efs)$concordance[1]
hr_efs <- exp(coef(cox_model_efs))[1]
hr_ci_efs <- exp(confint(cox_model_efs))[1,]
hr_label_efs <- paste0("HR: ", round(hr_efs, 2), " (95% CI: ", round(hr_ci_efs[1], 2), " - ", round(hr_ci_efs[2], 2), ")")
c_index_label_efs <- paste0("C-index: ", round(c_index_efs, 3))

# EFS Plot
p5 <- ggsurvplot(fit_efs, data = data, pval = TRUE, conf.int = FALSE,
                 risk.table = FALSE, font.legend = 13,
                 legend.title = "RF8 score",
                 legend = c(0.75, 0.35),
                 xlab = 'Survival time (days)', 
                 ylab = 'Survival probability (EFS)',
                 legend.labs = legend_labs) +
  labs(caption = paste0(c_index_label_efs, "\n", hr_label_efs))
p5

# Save EFS plot
pdf("figure6/02.cohort2.RF8_EFS.pdf", width = 5, height = 4.6)
p5
dev.off()

# OS - Calculate C-index and HR
fit_os <- survfit(Surv(OS_Time, OS_Status) ~ Predicted_group, data = data)
cox_model_os <- coxph(Surv(OS_Time, OS_Status) ~ Predicted_group, data = data) 
c_index_os <- summary(cox_model_os)$concordance[1]
hr_os <- exp(coef(cox_model_os))[1]
hr_ci_os <- exp(confint(cox_model_os))[1,]
hr_label_os <- paste0("HR: ", round(hr_os, 2), " (95% CI: ", round(hr_ci_os[1], 2), " - ", round(hr_ci_os[2], 2), ")")
c_index_label_os <- paste0("C-index: ", round(c_index_os, 3))

# OS Plot
p6 <- ggsurvplot(fit_os, data = data, pval = TRUE, conf.int = FALSE,
                 risk.table = FALSE, font.legend = 13,
                 legend.title = "RF8 score",
                 legend = c(0.75, 0.35),
                 xlab = 'Survival time (days)', 
                 ylab = 'Survival probability (OS)',
                 legend.labs = legend_labs) +
  labs(caption = paste0(c_index_label_os, "\n", hr_label_os))
p6

# Save OS plot
pdf("figure6/02.cohort2.RF8_OS.pdf", width = 5, height = 4.6)
p6
dev.off()

####################################
### Barplot for CR/CRi
####################################

# Load necessary libraries for barplot
library(readxl)
library(patchwork)

# Load response data for Cohort 2
data <- read.csv("figure6/RF8_Cohort2_survival.csv", header = TRUE)

# Create predicted groups based on RF8 probability
data$Predicted_group <- ifelse(data$RF8.prob.CR > 0.9, "0.9-1",
                               ifelse(data$RF8.prob.CR > 0.6, "0.6-0.9",
                                      ifelse(data$RF8.prob.CR > 0.3, "0.3-0.6", "0.1-0.3")))

# Summarize CR/CRi rates by predicted group
response_summary <- aggregate(Response ~ Predicted_group, data, function(x) {
  mean(x == "CR") * 100
})
colnames(response_summary)[2] <- "CR_Percentage"

# Create the barplot
p7 <- ggplot(response_summary, aes(x = Predicted_group, y = CR_Percentage, fill = Predicted_group)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = paste0(round(CR_Percentage, 1), "%")), vjust = -0.5) +
  labs(x = "Groups stratified by RF8 scores", y = "CR/CRi (%)") +
  ylim(0, 100) +
  theme_classic() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 0, hjust = 1)) 

# Display the barplot
p7

# Save the barplot
ggsave(p7, filename = "figure6/04.barplot_CR_rate.pdf", width = 5, height = 4)


