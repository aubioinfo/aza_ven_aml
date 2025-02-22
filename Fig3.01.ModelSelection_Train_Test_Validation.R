################################################
## Benchmark ML Methods for Prediction
################################################

# Set working directory
setwd("D:\\01.Projects\\01.project\\35_AZA_VEN\\00.RNAseq_Model")

library(mlr3)
library(mlr3learners)
library(mlr3verse)
library(mlr3tuning)
library(mlr3measures)
library(data.table)
library(mlr3benchmark)
library(Boruta)
library(caret)
library(randomForest)
library(pROC)
library(PRROC)
library(readxl)

# Set seed for reproducibility
set.seed(123456)

################################################
## Load Dataset and Prepare Classification Task
################################################
# Read gene expression data
data <- readxl::read_excel("TrainingTest_RJAML_Cohort1_84_genes.xlsx")
data <- data[, -1]  # Remove the first column if unnecessary (PatientID column)
data$Response <- as.factor(data$Response)  # Convert 'Response' to factor for classification

# Create classification task for mlr3
task <- TaskClassif$new(id = "AML", backend = data, target = "Response")

################################################
## Define Learners (Models) and Resampling Strategy
################################################
# Define a list of learners (classification models)
learners <- list(
  lrn("classif.rpart", predict_type = "prob"),       # Decision Tree
  lrn("classif.ranger", predict_type = "prob"),      # Random Forest
  lrn("classif.xgboost", predict_type = "prob"),     # XGBoost
  lrn("classif.log_reg", predict_type = "prob"),     # Logistic Regression
  lrn("classif.svm", predict_type = "prob"),         # Support Vector Machine
  lrn("classif.nnet", predict_type = "prob"),        # Neural Network
  lrn("classif.kknn", predict_type = "prob"),        # K-Nearest Neighbors
  lrn("classif.naive_bayes", predict_type = "prob"), # Naive Bayes
  lrn("classif.glmnet", predict_type = "prob"),      # ElasticNet
  lrn("classif.cv_glmnet", predict_type = "prob")    # Cross-validated ElasticNet
)

# Define a 5-fold cross-validation resampling strategy repeated 100 times
resampling_repeated <- rsmp("repeated_cv", folds = 5, repeats = 100)

################################################
## Benchmarking Setup and Execution
################################################
# Create benchmarking design grid
design <- benchmark_grid(
  tasks = task,
  learners = learners,
  resamplings = resampling_repeated
)

# Execute benchmark
bmr <- benchmark(design)

# Define performance measures
measures <- list(msr("classif.auc"), msr("classif.acc"))

# Aggregate performance results
mean_score <- bmr$aggregate(measures)

# Extract benchmark scores for each iteration
score <- bmr$score(measures)
score <- score[, c("learner_id", "iteration", "classif.auc", "classif.acc")]

# Calculate mean and median AUC
mean_auc <- mean(score$classif.auc, na.rm = TRUE)
median_auc <- median(score$classif.auc, na.rm = TRUE)

# Save benchmark results to a CSV file
write.csv(score, "figure3/ML_Results_GeneExp.csv")

################################################
## Boruta Feature Selection
################################################

# Reload dataset and prepare for feature selection
data <- read_excel("TrainingTest_RJAML_Cohort1_84_genes.xlsx")
data <- data[, -1]  # Remove the first column if unnecessary
data$Response <- as.factor(data$Response)  # Convert 'Response' to factor

# Perform Boruta feature selection
set.seed(123)
boruta_result <- Boruta(Response ~ ., data = data, pValue = 0.01, maxRuns = 1000, doTrace = 2)
boruta_selected <- TentativeRoughFix(boruta_result)

# Extract selected features (excluding tentative features)
selected_features <- getSelectedAttributes(boruta_selected, withTentative = FALSE)
selected_features

#######################################################
## Train a Random Forest Model with Selected Features
#######################################################

# Reload dataset and prepare with selected features
data <- read_excel("TrainingTest_RJAML_Cohort1_84_genes.xlsx")
data <- data[, -1]  # Remove PatientID column
data_selected <- data[, c("Response", selected_features)]  # Select important features

# Log-transform feature values (excluding the Response column)
data_selected[, -1] <- log2(data_selected[, -1] + 1)
data_selected$Response <- as.factor(data_selected$Response)  # Ensure 'Response' is a factor

# Split the dataset: 70% for training, 30% for testing
set.seed(12345)
trainIndex <- createDataPartition(data_selected$Response, p = 0.7, list = FALSE)
trainData <- data_selected[trainIndex, ]
testData <- data_selected[-trainIndex, ]

# Define training control with 5-fold cross-validation repeated 50 times
control <- trainControl(method = "repeatedcv", number = 5, repeats = 50,
                        classProbs = TRUE, summaryFunction = twoClassSummary)

# Define hyperparameter tuning grid
tuneGrid <- expand.grid(.mtry = c(1, 3, 5))

# Train a Random Forest model with hyperparameter tuning
rf_tuned <- train(
  Response ~ ., 
  data = trainData, 
  method = "rf", 
  trControl = control, 
  tuneGrid = tuneGrid, 
  ntree = 50,  # Number of trees
  metric = "ROC"  # Use ROC as the evaluation metric
)

# Save the trained Random Forest model to an RDS file
saveRDS(rf_tuned, file = "figure3/AZA_VEN_RandomForest8.rds")

# Extract variable importance from the model
importance <- varImp(rf_tuned, scale = FALSE)

# Convert importance to a data frame for easier visualization
importance_df <- data.frame(
  Feature = rownames(importance$importance),
  Importance = importance$importance$Overall
)

# Sort features by importance
importance_df <- importance_df[order(importance_df$Importance, decreasing = TRUE), ]

# Save feature importance to a CSV file
write.csv(importance_df, "figure3/RF_FeatureImportance.csv")

#########################################################
## Predict Class Labels on the Test Set and Plot ROC
#########################################################

library(pROC)
library(PRROC)

# Predict class labels and probabilities on the test set
test_predictions <- predict(rf_tuned, newdata = testData)
test_probabilities <- predict(rf_tuned, newdata = testData, type = "prob")

# Generate confusion matrix
confusionMatrix(test_predictions, testData$Response)

# Add predicted probabilities for CR
testData$RF8.prob.CR <- test_probabilities[, "CR"]
write.csv(testData, "figure3/RF8_RJAML_cohort1_test.csv")

# Plot ROC curve for the test set
roc_obj <- roc(testData$Response, testData$RF8.prob.CR)

# Save ROC plot to PDF
pdf("figure3/RF8_Test_ROC.pdf", width = 4.2, height = 4)
plot(smooth(roc_obj), legacy.axes = TRUE, lwd = 2,                          
     main = "Test set (N = 32)",
     xlab = "False Positive Rate (1 - Specificity)", 
     ylab = "True Positive Rate (Sensitivity)", 
     cex.main = 1.5, cex.lab = 1.2, cex.axis = 1.1)
abline(h = seq(0, 1, by = 0.1), v = seq(0, 1, by = 0.1), col = "gray80", lty = 2)
legend("bottomright", legend = c(paste("AUC =", round(roc_obj$auc, 2))),
       col = "#1f77b4", lwd = 2, bty = "n", cex = 1.2)
dev.off()

# Compute Precision-Recall curve and AUPRC for the test set
testData$Response <- factor(testData$Response, levels = c("NonCR", "CR"))
response_binary <- as.numeric(testData$Response) - 1
pr_curve <- pr.curve(scores.class0 = testData$RF8.prob.CR, weights.class0 = response_binary, curve = TRUE)

# Save AUPRC plot to PDF
pdf("figure3/RF8_Test_AUPRC.pdf", width = 4.2, height = 4.5)
plot(pr_curve, main = "Test set (N = 32)", color = "#1f77b4", auc.main = FALSE)
abline(h = seq(0, 1, by = 0.1), v = seq(0, 1, by = 0.1), col = "gray80", lty = 2)
legend("bottomright", legend = c(paste("AUPRC =", round(pr_curve$auc.integral, 2))),
       col = "#1f77b4", lwd = 2, bty = "n", cex = 1.2)
dev.off()

#########################################################
## Predict Class Labels on the Training Set and Plot ROC
#########################################################

# Predict class labels and probabilities on the training set
train_predictions <- predict(rf_tuned, newdata = trainData)
train_probabilities <- predict(rf_tuned, newdata = trainData, type = "prob")

# Generate confusion matrix
confusionMatrix(train_predictions, trainData$Response)

# Add predicted probabilities for CR
trainData$RF8.prob.CR <- train_probabilities[, "CR"]

# Plot ROC curve for the training set
roc_obj_train <- roc(trainData$Response, trainData$RF8.prob.CR)
roc_obj_train$auc

# Save ROC plot to PDF
pdf("figure3/RF8_Training_ROC.pdf", width = 4.2, height = 4)
plot(roc_obj_train, legacy.axes = TRUE, lwd = 2,                          
     main = "Training set (N = 78)",
     xlab = "False Positive Rate (1 - Specificity)", 
     ylab = "True Positive Rate (Sensitivity)", 
     cex.main = 1.5, cex.lab = 1.2, cex.axis = 1.1)
abline(h = seq(0, 1, by = 0.1), v = seq(0, 1, by = 0.1), col = "gray80", lty = 2)
legend("bottomright", legend = c(paste("AUC =", round(roc_obj_train$auc, 2))),
       col = "#1f77b4", lwd = 2, bty = "n", cex = 1.2)
dev.off()

#########################################################
## Predict Class Labels on the Entire Dataset and Plot ROC
#########################################################

# Predict class labels and probabilities on the entire set
entire_predictions <- predict(rf_tuned, newdata = data_selected)
entire_probabilities <- predict(rf_tuned, newdata = data_selected, type = "prob")

# Generate confusion matrix
confusionMatrix(entire_predictions, data_selected$Response)

# Add predicted probabilities and class labels
data_selected$RF8.prob.CR <- entire_probabilities[, "CR"]
data_selected$Predicted_Response <- entire_predictions
write.csv(data_selected, "figure3/RF8_RJAML_cohort1_entire.csv")

# Plot ROC curve for the entire dataset
roc_obj_entire <- roc(data_selected$Response, data_selected$RF8.prob.CR)
roc_obj_entire$auc

# Save ROC plot to PDF
pdf("figure3/RF8_Entire_ROC.pdf", width = 4.2, height = 4)
plot(roc_obj_entire, legacy.axes = TRUE, lwd = 2,                          
     main = "Entire set (N = 110)",
     xlab = "False Positive Rate (1 - Specificity)", 
     ylab = "True Positive Rate (Sensitivity)", 
     cex.main = 1.5, cex.lab = 1.2, cex.axis = 1.1)
plot(roc_obj_entire, print.thres = "best", col = "#1f77b4", print.thres.cex = 1.2, add = TRUE)
abline(h = seq(0, 1, by = 0.1), v = seq(0, 1, by = 0.1), col = "gray80", lty = 2)
legend("bottomright", legend = c(paste("AUC =", round(roc_obj_entire$auc, 2))),
       col = "#1f77b4", lwd = 2, bty = "n", cex = 1.2)
dev.off()

# Compute Precision-Recall curve and AUPRC for the entire dataset
data_selected$Response <- factor(data_selected$Response, levels = c("NonCR", "CR"))
response_binary_entire <- as.numeric(data_selected$Response) - 1
pr_curve_entire <- pr.curve(scores.class0 = data_selected$RF8.prob.CR, weights.class0 = response_binary_entire, curve = TRUE)

# Save AUPRC plot to PDF
pdf("figure3/RF8_Entire_AUPRC.pdf", width = 4.2, height = 4.5)
plot(pr_curve_entire, main = "Entire set (N = 110)", color = "#1f77b4", auc.main = FALSE)
abline(h = seq(0, 1, by = 0.1), v = seq(0, 1, by = 0.1), col = "gray80", lty = 2)
legend("bottomright", legend = c(paste("AUPRC =", round(pr_curve_entire$auc.integral, 2))),
       col = "#1f77b4", lwd = 2, bty = "n", cex = 1.2)
dev.off()

#########################################################
## AUC and AUPRC Heatmap for the Entire and Test Sets
#########################################################
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(gridExtra)

# Prepare data for plotting: AUC and AUPRC for both Entire and Test datasets
auc_entire_data <- auc_auprc_entire[, c("Feature", "AUC")]
auprc_entire_data <- auc_auprc_entire[, c("Feature", "AUPRC")]

auc_test_data <- auc_auprc_test[, c("Feature", "AUC")]
auprc_test_data <- auc_auprc_test[, c("Feature", "AUPRC")]

# Melt the data for ggplot
auc_entire_melt <- melt(auc_entire_data, id.vars = "Feature")
auprc_entire_melt <- melt(auprc_entire_data, id.vars = "Feature")

auc_test_melt <- melt(auc_test_data, id.vars = "Feature")
auprc_test_melt <- melt(auprc_test_data, id.vars = "Feature")

# Function to generate heatmaps
plot_heatmap <- function(data, title) {
  ggplot(data, aes(x = variable, y = Feature, fill = value)) +
    geom_tile(color = "white") +
    geom_text(aes(label = round(value, 2)), size = 4) +
    scale_fill_gradientn(colors = rev(brewer.pal(9, "BrBG")), limits = c(0.2, 1)) +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    ggtitle(title)
}

# Generate heatmaps for the Entire and Test datasets
auc_entire_plot <- plot_heatmap(auc_entire_melt, "AUC (Entire)")
auprc_entire_plot <- plot_heatmap(auprc_entire_melt, "AUPRC (Entire)")

auc_test_plot <- plot_heatmap(auc_test_melt, "AUC (Test)")
auprc_test_plot <- plot_heatmap(auprc_test_melt, "AUPRC (Test)")

# Combine the four heatmaps into a 2x2 grid
pdf("figure3/auc_auprc_plots.pdf", width = 6, height = 6)
grid.arrange(auc_entire_plot, auprc_entire_plot, auc_test_plot, auprc_test_plot, ncol = 2)
dev.off()


