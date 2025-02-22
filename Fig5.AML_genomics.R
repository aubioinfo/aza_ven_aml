library(tidyverse)
library(ggpubr)
library(ggsci)
library(ggridges)
library(cowplot)
library(ggrastr)
library(readxl)
library(ggrastr)

# Single mutation
hcat <- read_xlsx("figure5//rjaml_mutations_rf8.xlsx")
mut_variables <- c("RF8.prob.CR", "TP53",	"ASXL1",	"BCOR",	"BCORL1",
                   "CBL",	"CEBPA",	"ETV6",	"DDX41",	"DNMT3A",	"EZH2",	
                   "FLT3-ITD",	"FLT3-TKD",	"IDH1",	"IDH2",	"KRAS",	"NF1","NPM1",
                   "NRAS",	"PHF6",	"PTPN11",	"RUNX1",	"SF3B1",	"SMC1A",	"SRSF2",
                   "STAG2",	"TET2",	"U2AF1",	"WT1",	"CBFB::MYH11",	"RUNX1::RUNX1T1"
)
hcat <- hcat[, mut_variables]

# Transform data for plotting
alteration <- hcat %>%
  gather(-RF8.prob.CR, key='Alteration', value='Status') %>%
  filter(Status > 0)

p1 <- alteration %>%
  ggplot(aes(y=fct_reorder(Alteration, RF8.prob.CR), x=RF8.prob.CR, fill=stat(x))) +
  geom_density_ridges_gradient(
    jittered_points = TRUE, scale=1.25,
    position=position_points_jitter(width=0.05, height=0),
    point_shape='|', point_size=3, point_alpha=0.3) +
  scale_fill_gradient2(low="grey", high="#9970AB", midpoint=0.5, name="") +
  theme_pubr() +
  ylab('') +
  xlab('RF8 score') +
  theme(legend.position='right',
        axis.title.x=element_text(size=15))
ggsave(p1, filename = "02.mutations_rf8.pdf", width = 6.5, height = 11)

######
## Obtain mutation pairs in AML patients
library(dplyr)
library(tidyr)
library(readr)

combos <- read_xlsx("figure5//mutations_pairs_rf8.xlsx") %>% 
  column_to_rownames(var = "PatientID")

# Select columns that contain 'mut' and replace values
df <- combos %>%
  select(contains("mut")) %>%
  mutate(across(everything(), ~ replace_na(as.character(.), "NA"))) %>%
  # mutate(across(everything(), ~ recode(., "Positive" = "1", "Negative" = "0", "NA" = "0"))) %>%
  mutate(across(everything(), as.integer))

# Reorder by frequency and filter rows and columns
df <- df %>%
  select(order(colSums(df), decreasing = TRUE)) %>%
  filter(rowSums(.) >= 2) %>%
  select(where(~ sum(.) > 1))

# Generate combinations of columns
combinations <- combn(colnames(df), 2)

# Compute product of column pairs and bind them together
df_comb <- purrr::map_dfc(seq_len(ncol(combinations)), function(i) {
  combo <- combinations[, i]
  col1 <- df[[combo[1]]]
  col2 <- df[[combo[2]]]
  setNames(as.data.frame(col1 * col2), paste(combo, collapse = " + "))
})

row.names(df_comb) <- row.names(df)

# Filter columns by frequency and reset index
combo2 <- df_comb %>%
  select(where(~ sum(.) > 1)) %>%
  rownames_to_column(var = "PatientID")
# write_csv(combo2, "figure5/Two_mutations_combination_riskScore_final.csv",
#           col_names = TRUE)

combo2 <- combo2 %>% gather(-PatientID, key='Mutation', value='Present') %>%
  filter(Present > 0) %>%
  mutate(dnmt3a_dup = if_else((Mutation %>% str_detect("DNMT3A \\(")) & (Mutation %>% str_detect("DNMT3A_mut")), 'TRUE', 'FALSE'),
         flt3_dup = if_else((Mutation %>% str_detect("FLT3 \\(")) & (Mutation %>% str_detect("FLT3_mut")), 'TRUE', 'FALSE')) %>%
  filter((dnmt3a_dup == 'FALSE') & (flt3_dup == 'FALSE')) %>%
  select(-dnmt3a_dup, -flt3_dup) %>%
  mutate(Mutation = Mutation %>% str_replace_all('_mut', ''))


a <- combo2 %>% group_by(Mutation) %>%
  summarise(Total = sum(Present)) %>%
  ungroup() %>% arrange(-Total) #%>% head(50)
a

write_csv(a, "figure5/Two_mutations_combination_Freq_riskScore_final.csv")


## Plot mutation pairs
dat <- read.csv("figure5/Two_mutations_combination_riskScore_final.csv",
                header = T, check.names = FALSE)

# Initialize empty data frame for statistics
stat <- data.frame()

# Calculate statistics for each mutation
for (i in 3:ncol(dat)) {
  tmp <- dat[, c(2, i)]
  colnames(tmp)[2] <- "Mutation"
  
  tmp_0 <- filter(tmp, Mutation == 0)
  tmp_1 <- filter(tmp, Mutation == 1)
  
  testout <- wilcox.test(tmp_0$RiskScore, tmp_1$RiskScore, alternative = "two.sided")
  
  tmp_data <- data.frame(
    Negative_median = median(tmp_0$RiskScore, na.rm = TRUE),
    Positive_median = median(tmp_1$RiskScore, na.rm = TRUE),
    P.value = testout$p.value,
    Mutation = colnames(dat)[i]
  )
  
  stat <- rbind(stat, tmp_data)
}

# Load frequency data and merge with statistics
freq <- read_csv("figure5/Two_mutations_combination_Freq_riskScore_final.csv")
stat <- merge(stat, freq, by = "Mutation")

# write_csv(stat, "figure5/Two_mutations_Mut_vs.WT_stat.csv")

# Plot two mutations
plot <- read_csv("figure5/Two_mutations_Mut_vs.WT_stat.csv") %>%
  filter(Total >= 5)

mutants <- dat %>%
  select(RiskScore, all_of(plot$Mutation))

alteration <- mutants %>%
  gather(key = 'Alteration', value = 'Status', -RiskScore) %>%
  filter(Status > 0)

p9 <- alteration %>%
  ggplot(aes(y = fct_reorder(Alteration, RiskScore), x = RiskScore, fill = stat(x))) +
  geom_density_ridges_gradient(
    jittered_points = TRUE, scale = 1.25,
    position = position_points_jitter(width = 0.05, height = 0),
    point_shape = '|', point_size = 3, point_alpha = 0.3
  ) +
  scale_fill_gradient2(low = "#01665E", high = "#8C510A", midpoint = 0.5, name = "") +
  theme_pubr() +
  ylab('') +
  xlab('Risk score') +
  theme(
    legend.position = 'right',
    axis.title.x = element_text(size = 15)
  )

# Print plot
print(p9)
ggsave(p9, filename = "figure5/04.two_Mutations_RF8.pdf", width = 7, height = 8)

##############################
# Barplot for real response
##############################

### CR NCR VS. Mut WT
library(ggstatsplot)
library(patchwork)
dat <- read.csv("figure5/multi-cohorts_4gene_data.csv", header = TRUE)
p1 <- ggbarstats(dat, y = TP53, x = Response_group, palette = "Paired", label = "counts")
p2 <- ggbarstats(dat, y = DNMT3A, x = Response_group, palette = "Paired", label = "counts")
p3 <- ggbarstats(dat, y = IDH2, x = Response_group, palette = "Paired", label = "counts")
p4 <- ggbarstats(dat, y = CEBPA, x = Response_group, palette = "Paired", label = "counts")

plots <- p1|p2|p3|p4
ggsave(plots, filename = "figure5/05.barplot_CR_NCR.pdf", width = 14, height = 4.5)


####################################
### Table one for multi-cohorts
####################################
library(tableone)
data <- read.csv("figure5/multi-cohorts_4gene_Clinicaldata.csv",
                 header = T, check.names = FALSE)
myVars <- c("Age", "Gender", "Secondary_AML", "Complex_Karyotype", 
            "FLT3-ITD", "NPM1", "KRAS", "NRAS", "CEBPA", "DNMT3A",
            "IDH2", "TP53")
catVars <- c("Gender", "Secondary_AML", "Complex_Karyotype",
             "FLT3-ITD", "NPM1", "KRAS", "NRAS", "CEBPA", "DNMT3A",
             "IDH2", "TP53")

tab <- CreateTableOne(vars = myVars, strata = c("Response_group", "Cohort"),
                      addOverall = TRUE, data = data, factorVars = catVars)
tab2 <- print(tab, showAllLevels = TRUE)
write.csv(tab2, file = "figure5/multi-cohorts_4gene_Clinicaldata_out.csv")

####################################
### Complexheatmap for mutation
####################################

# Load required libraries
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(readxl)
data <- read_xlsx("rjaml_mutations_rf8.xlsx")
# Define clinical and mutation variables
clinic_variables <- c("PatientID", "Group", "Age", "BM_blast", "Gender", 
                      "Secondary_AML", "Complex_Karyotype", "ELN")
mut_variables <- c("PatientID", "Group", "TP53",	"ASXL1",	"BCOR",	"BCORL1",
                   "CBL",	"CEBPA",	"ETV6",	"DDX41",	"DNMT3A",	"EZH2",	
                   "FLT3-ITD",	"FLT3-TKD",	"IDH1",	"IDH2",	"KRAS",	"NF1","NPM1",
                   "NRAS",	"PHF6",	"PTPN11",	"RUNX1",	"SF3B1",	"SMC1A",	"SRSF2",
                   "STAG2",	"TET2",	"U2AF1",	"WT1",	"CBFB::MYH11",	"RUNX1::RUNX1T1"
)

# Extract clinical and mutation data
clinical_data <- data[, clinic_variables]
mutation_data <- data[, mut_variables]

# Set PatientID as row names and remove from both mutation and clinical data
rownames(mutation_data) <- mutation_data$PatientID
mutation_data <- mutation_data[, -1]  # Remove PatientID column from mutation data

rownames(clinical_data) <- clinical_data$PatientID
clinical_data <- clinical_data[, -1]  # Remove PatientID column from clinical data

# Extract mutation matrix and transpose it
mutation_matrix <- mutation_data[, -1]  # Remove Response column
mutation_binary <- apply(mutation_matrix, 2, function(x) ifelse(x == 1, "Mutation", ""))
mutation_binary_transposed <- t(mutation_binary)

# Define custom colors for clinical variables
col_age <- colorRamp2(c(40, 80), c("#F7F7F7", "#E69F00"))
col_blast <- colorRamp2(c(20, 100), c("#F7F7F7", "#56B4E9"))

clinical_annotation <- HeatmapAnnotation(
  Group = clinical_data$Group,
  Age = clinical_data$Age,
  BM_blast = clinical_data$BM_blast,
  Gender = clinical_data$Gender,
  Secondary_AML = clinical_data$Secondary_AML, 
  Complex_Karyotype = clinical_data$Complex_Karyotype, 
  ELN = clinical_data$ELN, 
  show_annotation_name = TRUE,
  col = list(
    Age = col_age,
    BM_blast = col_blast,
    Gender = c("0" = "#F0E442", "1" = "#0072B2"),
    Secondary_AML = c("0" = "#56B4E9", "1" = "#E69F00"),
    Complex_Karyotype = c("0" = "#999999", "1" = "#D55E00"),
    ELN = c("Favorable" = "#7c9d97", "Intermediate" = "#9cb0c3", "Adverse" = "#e9b383"),
    Group = c("High" = "#7c9d97", "Low" = "#e9b383")
  )
)

# Define alteration function
alter_fun <- list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "mm"), h - unit(1, "mm"), gp = gpar(fill = '#ebebe1ff', col = NA))
  }, 
  Mutation = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "mm"), h - unit(1, "mm"), gp = gpar(fill = "#377EB8", col = NA))  # Mutation = Dark Blue
  }
)

# Color function for mutation heatmap
col <- c("Mutation" = "#377EB8")

# Generate oncoPrint heatmap
pdf("03.oncoplot_map.pdf", width = 16, height = 9.5)
oncoPrint(
  mutation_binary_transposed,
  alter_fun = alter_fun,
  col = col,
  column_split = clinical_data$Group,  # Split columns by Response (CR/CRi vs Non-CR/CRi)
  bottom_annotation = clinical_annotation,  # Add clinical annotations to the heatmap
  show_column_names = TRUE,
  show_row_names = TRUE,
  heatmap_legend_param = list(title = "Alterations")
)
dev.off()


#####################################
############## fisher test ##############
#####################################
# Load necessary libraries
library(tidyverse)
library(ggplot2)
library(ggsci)
library(scales)
library(ggbreak)

# Read input data
fisher.dat <- read_xlsx("rjaml_mutations_rf8.xlsx")
mut_variables <- c("RF8.prob.CR", "TP53",	"ASXL1",	"BCOR",	"BCORL1",
                   "CBL",	"CEBPA",	"ETV6",	"DDX41",	"DNMT3A",	"EZH2",	
                   "FLT3-ITD",	"FLT3-TKD",	"IDH1",	"IDH2",	"KRAS",	"NF1","NPM1",
                   "NRAS",	"PHF6",	"PTPN11",	"RUNX1",	"SF3B1",	"SMC1A",	"SRSF2",
                   "STAG2",	"TET2",	"U2AF1",	"WT1",	"CBFB::MYH11",	"RUNX1::RUNX1T1"
)
fisher.dat <- fisher.dat[, mut_variables]
# Set options
options(scipen = 100)

# Initialize matrix to store Fisher test results
stat1 <- matrix(NA, nrow = ncol(fisher.dat) - 1, ncol = 4)

# Perform Fisher test for each column (gene)
for (i in 2:ncol(fisher.dat)) {
  m1 <- fisher.dat[, c(1, i)]
  m1 <- m1[complete.cases(m1), ]
  m1$RF8.prob.CR <- ifelse(m1$RF8.prob.CR >= median(m1$RF8.prob.CR), "High", "Low")
  
  # Create contingency table
  m2 <- table(m1)
  m3 <- matrix(c(m2[1], m2[2], m2[3], m2[4]), ncol = 2, nrow = 2)
  
  # Perform Fisher test
  s1 <- fisher.test(m3)
  
  # Store results
  stat1[i - 1, ] <- c(round(s1$p.value, 5), as.numeric(round(s1$estimate, 3)), round(s1$conf.int[1], 3), round(s1$conf.int[2], 3))
}

# Convert results to data frame and add gene names
stat1 <- data.frame(stat1)
stat1$Gene <- colnames(fisher.dat)[2:ncol(fisher.dat)]
colnames(stat1) <- c("P.value", "Odds_Ratio", "Lower95%", "Upper95%", "Gene")

# Reorder columns and sort by p-value
stat <- stat1[, c("Gene", "P.value", "Odds_Ratio", "Lower95%", "Upper95%")]
stat <- stat[order(stat$P.value), ]

# Adjust p-values for false discovery rate
# stat$fdr <- p.adjust(stat$P.value, method = "fdr", n = nrow(stat))
# Write results to CSV
write.csv(stat, "Mutation_fisher_out.csv")

###################################
######  dotplot Fisher test  ######
###################################

# Read results data
tab <- read.csv("Mutation_fisher_out.csv", header = TRUE, row.names = 1, check.names = FALSE)

# Subset data for significant results
tab <- subset(tab, P.value < 0.05)
tab <- tab[order(tab$Odds_Ratio, decreasing = TRUE), ]
# Order genes by odds ratio
tab$Gene <- factor(tab$Gene, levels = rev(tab$Gene))

# Create dot plot
forest <- ggplot(data = tab, aes(x = Odds_Ratio, y = Gene)) +
  geom_point(aes(color = -log2(P.value), size = 0.75), shape = 19) +
  geom_errorbarh(aes(xmin = `Lower95%`, xmax = `Upper95%`), height = 0, color = "black") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
  theme_classic() +
  labs(x = "Odds Ratio (OR)", y = "") +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 1, color = "black"),
    panel.grid.major = element_line(size = 0.1, linetype = "dashed", color = "grey"),
    panel.grid.minor = element_line(size = 0.1, linetype = "dashed", color = "grey")
  )

# Customize plot color and scale
p2 <- forest + scale_color_gradient2(low = "white", mid = "white", high = "#0072B5FF") + scale_x_log10()

ggsave(p2, filename = "Mutations_Fisher.pdf", width = 4.5, height = 3.5)

## Barplot

library(ggstatsplot)
library(patchwork)
data <- read_xlsx("rjaml_mutations_rf8.xlsx")

p1 <- ggbarstats(data = data, x = TP53, y = Group)
p2 <- ggbarstats(data = data, x = DNMT3A, y = Group)
p3 <- ggbarstats(data = data, x = IDH2, y = Group)
p4 <- ggbarstats(data = data, x = CEBPA, y = Group)
plots <- p1|p2|p3|p4
ggsave(plots, filename = "mutation_freq_bar.pdf", width = 10, height = 4)




##  Multivariable LogisticReg
library(forestmodel)
library(ggplot2)
library(ggsci)
library(scales)
library(ggbreak)

#---------------------
# Load the Data
#---------------------
dat <- read_xlsx("rjaml_mutations_rf8.xlsx")
dat$Response <- ifelse(dat$Response == 'NonCR', 1, 0)
dat$ELN <- ifelse(dat$ELN == "Adverse", 1, 0)
# Define forest plot panels
forest.panels <- list(
  list(width = 0.01),
  list(width = 0.1, display = ~variable, fontface = "bold", heading = "Variable"),
  list(width = 0.1, display = ~level),
  list(width = 0.05, display = ~n, hjust = 1, heading = "N"),
  list(width = 0.03, item = "vline", hjust = 0.5),
  list(width = 0.55, item = "forest", hjust = 0.5, heading = "Odds ratio", linetype = "dashed", line_x = 0),
  list(width = 0.03, item = "vline", hjust = 0.5),
  list(width = 0.12, display = ~ ifelse(reference, "Reference", sprintf("%0.3f (%0.3f, %0.3f)", trans(estimate), trans(conf.low), trans(conf.high))), display_na = NA),
  list(width = 0.03, display = ~ ifelse(reference, "", format.pval(p.value, digits = 1, eps = 0.001)), display_na = NA, hjust = 1, heading = "p"),
  list(width = 0.03)
)

# Define function to create forest plots
create_forest_plot <- function(data) {
  formula <- as.formula(paste0('Response ~ ', paste(c('RF8.prob.CR', 'Age', 'Gender', 'BM_blast',	'Secondary_AML',	'Complex_Karyotype', 'ELN', 'TP53', 'IDH2', 'DNMT3A', 'CEBPA'), collapse = ' + ')))
  model <- glm(formula, data = data, family = binomial)
  print(forest_model(model, panels = forest.panels, factor_separate_line = FALSE))
}

pdf("forestplot.multivariable.pdf", width = 8, height = 5)
create_forest_plot(dat)
dev.off()

