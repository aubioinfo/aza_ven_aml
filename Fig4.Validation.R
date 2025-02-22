####################################
### Table 1 for RJAML-cohort 2
####################################

library(readxl)
library(tableone)

data <- read_xlsx("02.32_AZA_VEN_RNAseq_Clinical.xlsx", sheet = 1)
table(data$Best_response)

# Retain mutations mutated in >= 3 patients
variables <- c("Age", "BM_blast", "Best_response", "Gender", "Secondary_AML", "Complex_Karyotype", "ELN")
grouping <- "Response"

# Define categorical variables
factorVars <- c("Best_response", "Gender", "Secondary_AML", "Complex_Karyotype", "ELN")

# Create Table 1 object
table1 <- CreateTableOne(vars = variables, strata = grouping, data = data, factorVars = factorVars, addOverall = TRUE)
table1
write.csv(print(table1), "figure4/rjaml-cohort2_baseline.csv")

################################################
## Predict on External Cohort 2
################################################

# Load the trained model
rf_tuned <- readRDS("figure3/AZA_VEN_RandomForest8.rds")

# Load the new data for Cohort 2
newDat <- read_xlsx("../rjaml_test_all.xlsx")
newDat <- newDat[, -1]  # Remove the first column (Patient ID)

# Select features and log-transform
data_selected <- newDat[, c("Response", selected_features)]
data_selected[, -1] <- log2(data_selected[, -1] + 1)

# Predict class labels and probabilities for Cohort 2
new_predictions <- predict(rf_tuned, newdata = data_selected)
new_probabilities <- predict(rf_tuned, newdata = data_selected, type = "prob")

# Add predictions to the dataset
newDat$Predicted_Response <- new_predictions
newDat$RF8.prob.CR <- new_probabilities[, "CR"]
write.csv(newDat, "figure3/RF8_Cohort2.csv")

# Compute and plot ROC for Cohort 2
roc_obj <- roc(newDat$Response, newDat$RF8.prob.CR)
roc_obj$auc

pdf("figure3/RF8_Cohort2_ROC.pdf", width = 4.2, height = 4)
plot(smooth(roc_obj), legacy.axes = TRUE, lwd = 2,
     main = "Cohort 2 (N = 32)",
     xlab = "False Positive Rate (1 - Specificity)", 
     ylab = "True Positive Rate (Sensitivity)", 
     cex.main = 1.5, cex.lab = 1.2, cex.axis = 1.1)
abline(h = seq(0, 1, by = 0.1), v = seq(0, 1, by = 0.1), col = "gray80", lty = 2)
legend("bottomright", legend = c(paste("AUC =", round(roc_obj$auc, 2))),
       col = "#1f77b4", lwd = 2, bty = "n", cex = 1.2)
dev.off()

# Compute and plot Precision-Recall (PR) curve for Cohort 2
newDat$Response <- factor(newDat$Response, levels = c("NonCR", "CR"))
response_binary <- as.numeric(newDat$Response) - 1
pr_curve <- pr.curve(scores.class0 = newDat$RF8.prob.CR, weights.class0 = response_binary, curve = TRUE)

pdf("figure3/RF8_Cohort2_AUPRC.pdf", width = 4.2, height = 4.5)
plot(pr_curve, main = "Cohort 2 (N = 32)", color = "#1f77b4", auc.main = FALSE)
abline(h = seq(0, 1, by = 0.1), v = seq(0, 1, by = 0.1), col = "gray80", lty = 2)
legend("bottomright", legend = c(paste("AUPRC =", round(pr_curve$auc.integral, 2))),
       col = "#1f77b4", lwd = 2, bty = "n", cex = 1.2)
dev.off()

################################################
## Predict on External FPMTB Set (Ex Vivo)
################################################

# Load the new data for the FPMTB set
newDat <- read_xlsx("Validation_Exvivo_FPMTB_84Genes_43R_55S.xlsx")
newDat <- newDat[, -1]  # Remove the first column (Patient ID)

# Select features and log-transform
data_selected <- newDat[, c("Response", selected_features)]
data_selected[, -1] <- log2(data_selected[, -1] + 1)

# Predict class labels and probabilities for the FPMTB set
new_predictions <- predict(rf_tuned, newdata = data_selected)
new_probabilities <- predict(rf_tuned, newdata = data_selected, type = "prob")

# Add predictions to the dataset
newDat$Predicted_Response <- new_predictions
newDat$RF8.prob.CR <- new_probabilities[, "CR"]
write.csv(newDat, "figure3/RF8_FPMTB.csv")

# Compute and plot ROC for the FPMTB set
roc_obj <- roc(newDat$Response, newDat$RF8.prob.CR)
roc_obj$auc

pdf("figure3/RF8_FPMTB_ROC.pdf", width = 4.2, height = 4)
plot(smooth(roc_obj), legacy.axes = TRUE, lwd = 2,
     main = "FPMTB set (N = 98)",
     xlab = "False Positive Rate (1 - Specificity)", 
     ylab = "True Positive Rate (Sensitivity)", 
     cex.main = 1.5, cex.lab = 1.2, cex.axis = 1.1)
abline(h = seq(0, 1, by = 0.1), v = seq(0, 1, by = 0.1), col = "gray80", lty = 2)
legend("bottomright", legend = c(paste("AUC =", round(roc_obj$auc, 2))),
       col = "#1f77b4", lwd = 2, bty = "n", cex = 1.2)
dev.off()

# Compute and plot Precision-Recall (PR) curve for the FPMTB set
newDat$Response <- factor(newDat$Response, levels = c("NonCR", "CR"))
response_binary <- as.numeric(newDat$Response) - 1
pr_curve <- pr.curve(scores.class0 = newDat$RF8.prob.CR, weights.class0 = response_binary, curve = TRUE)

pdf("figure3/RF8_FPMTB_AUPRC.pdf", width = 4.2, height = 4.5)
plot(pr_curve, main = "FPMTB set (N = 98)", color = "#1f77b4", auc.main = FALSE)
abline(h = seq(0, 1, by = 0.1), v = seq(0, 1, by = 0.1), col = "gray80", lty = 2)
legend("bottomright", legend = c(paste("AUPRC =", round(pr_curve$auc.integral, 2))),
       col = "#1f77b4", lwd = 2, bty = "n", cex = 1.2)
dev.off()

################################################
## Predict on External BeatAML Set (Ex Vivo)
################################################

# Load the new data for the BeatAML set
newDat <- read_xlsx("Validation_Exvivo_BeatAML_84Genes_170R_88S.xlsx")
newDat <- newDat[, -1]  # Remove the first column (Patient ID)

# Select features and log-transform
data_selected <- newDat[, c("Response", selected_features)]
data_selected[, -1] <- log2(data_selected[, -1] + 1)

# Predict class labels and probabilities for the BeatAML set
new_predictions <- predict(rf_tuned, newdata = data_selected)
new_probabilities <- predict(rf_tuned, newdata = data_selected, type = "prob")

# Add predictions to the dataset
newDat$Predicted_Response <- new_predictions
newDat$RF8.prob.CR <- new_probabilities[, "CR"]
write.csv(newDat, "figure3/RF8_BeatAML.csv")

# Compute and plot ROC for the BeatAML set
roc_obj <- roc(newDat$Response, newDat$RF8.prob.CR)
roc_obj$auc

pdf("figure3/RF8_BeatAML_ROC.pdf", width = 4.2, height = 4)
plot(smooth(roc_obj), legacy.axes = TRUE, lwd = 2,
     main = "BeatAML set (N = 258)",
     xlab = "False Positive Rate (1 - Specificity)", 
     ylab = "True Positive Rate (Sensitivity)", 
     cex.main = 1.5, cex.lab = 1.2, cex.axis = 1.1)
abline(h = seq(0, 1, by = 0.1), v = seq(0, 1, by = 0.1), col = "gray80", lty = 2)
legend("bottomright", legend = c(paste("AUC =", round(roc_obj$auc, 2))),
       col = "#1f77b4", lwd = 2, bty = "n", cex = 1.2)
dev.off()

# Compute and plot Precision-Recall (PR) curve for the BeatAML set
newDat$Response <- factor(newDat$Response, levels = c("NonCR", "CR"))
response_binary <- as.numeric(newDat$Response) - 1
pr_curve <- pr.curve(scores.class0 = newDat$RF8.prob.CR, weights.class0 = response_binary, curve = TRUE)

pdf("figure3/RF8_BeatAML_AUPRC.pdf", width = 4.2, height = 4.5)
plot(pr_curve, main = "BeatAML set (N = 258)", color = "#1f77b4", auc.main = FALSE)
abline(h = seq(0, 1, by = 0.1), v = seq(0, 1, by = 0.1), col = "gray80", lty = 2)
legend("bottomright", legend = c(paste("AUPRC =", round(pr_curve$auc.integral, 2))),
       col = "#1f77b4", lwd = 2, bty = "n", cex = 1.2)
dev.off()
