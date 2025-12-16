########
#
# Analysis of Bulk RNA-Seq Data
#
# Simon Grund Sorensen, Jakob Skou Pedersen, Søren Besenbacher, Farhad Zamani
# Aarhus University
#
########

################################################################################
#### PART 1: Loading Data and Initial Visualizations ####
################################################################################

################################################################################
#### SECTION 1: Loading Packages ####
################################################################################

# Load required packages for data manipulation and modeling:
library(tidyverse)
library(tidymodels)

################################################################################
#### SECTION 2: Loading and Preparing the Data ####
################################################################################

# Load the TCGA formatted expression data:
d <- readRDS("Data/TCGA_Formatted_data.rds")

# We select two cancer types to compare (Breast Cancer vs. Colorectal Adenocarcinoma):
d <- filter(d, Cancertype %in% c("BRCA", "COAD"))

# Clean the Sample ID (remove extra characters):
d <- mutate(d, SampleID = str_sub(Sample, end = -14))

# Rearrange columns: put SampleID first, remove old Sample and Cancertype columns:
d <- dplyr::select(d, SampleID, Cancertype, everything()) %>%
  dplyr::select(-c(Sample, Cancertype))

# Load the metadata file:
meta <- read_delim("Data/TCGA_Metadata.tsv")

# Clean column names and select relevant metadata variables:
meta <- dplyr::rename(meta, SampleID = `Sample ID`) %>%
  dplyr::select(
    SampleID,
    Study  = `TCGA PanCanAtlas Cancer Type Acronym`,
    Age    = `Diagnosis Age`,
    Gender = Sex,
    Type   = `Sample Type`
  ) %>%
  filter(SampleID %in% d$SampleID)

# How many samples do we have metadata for?
table(d$SampleID %in% meta$SampleID)
# Note: 43 samples are missing metadata – that is OK for now.

################################################################################
#### SECTION 3: Exploring the RNA-Seq Count Data ####
################################################################################

# Check out the RNA read counts:
glimpse(d)

# EXERCISE A: Data Structure
# What is in the rows? What is in the columns?
# Write your answer here:
#
#

################################################################################
#### SECTION 4: Summarizing Gene Expression ####
################################################################################

# Reformat to long format and count the number of samples with expression > 0
# for each gene:
Summary_d <- d %>%
  pivot_longer(cols = -SampleID) %>%
  group_by(name) %>%
  mutate(n_above_0 = sum(value > 0)) %>%
  distinct(name, n_above_0) %>%
  dplyr::rename(gene_id = name) %>%
  ungroup()

# EXERCISE B: Understanding the Code
# Comment on each line of the code block above. What does each step do?
# Discuss with a partner or write your explanation here:
#
#

# Visualize the distribution of expressed genes:
ggplot(Summary_d, aes(x = n_above_0)) +
  geom_histogram() +
  theme_minimal(base_size = 30) +
  geom_rug(sides = "b")

# EXERCISE C: Zero-Expression Genes
# How many genes have 0 expression across all samples?
# Write your answer here:
#
#

################################################################################
#### SECTION 5: Filtering Low-Expression Genes ####
################################################################################

# We keep only genes that are expressed in at least 25% of samples:
genes_to_include <- filter(Summary_d, n_above_0 > 0.25 * ncol(d))$gene_id

# Subset the data to include only these genes:
d <- select(d, any_of(c("SampleID", "Cancertype")), any_of(genes_to_include))

# EXERCISE D: Remaining Genes
# How many genes remain after filtering?
# Write your answer here:
#
#

################################################################################
#### SECTION 6: Normalization and Log-Transformation ####
################################################################################

# We normalize the counts using log-transformed Counts Per Million (logCPM).
# For each gene, we calculate the sum of reads across all samples in millions.

d <- d %>%
  pivot_longer(cols = -SampleID) %>%
  group_by(name) %>%
  mutate(
    sample_summed_counts_in_million = sum(value) / 1e6
  ) %>%
  ungroup() %>%
  mutate(
    CPM    = (value + 1) / sample_summed_counts_in_million,
    logCPM = log(CPM)
  )

# Remove intermediate column to clean up:
d$sample_summed_counts_in_million <- NULL

# Check the distribution of log-transformed counts:
ggplot(d, aes(x = logCPM)) +
  geom_histogram() +
  theme_minimal(base_size = 30)

################################################################################
#### SECTION 7: Adding Metadata and Visualization ####
################################################################################

# Add metadata annotation to the expression data:
d <- left_join(meta, d)

# Make a boxplot of logCPM for the two Study groups:
ggplot(d, aes(x = Study, y = logCPM)) +
  geom_boxplot() +
  theme_minimal(base_size = 30) +
  labs(title = "Gene Expression by Cancer Type",
       x = "Cancer Type",
       y = "logCPM")   



# EXERCISE E: Interpreting the Boxplot
# What does this graph tell us about the two cancer types?
# Write your answer here:
#
#

################################################################################
#### SECTION 8: Comparing Two Samples ####
################################################################################

# Let's plot the expression of two individual samples against each other:
s1 <- filter(d, SampleID == "TCGA-AA-3527-01")
s2 <- filter(d, SampleID == "TCGA-D8-A1JK-01")

ggplot(data.frame(x = s1$logCPM, y = s2$logCPM), aes(x = x, y = y)) +
  geom_point() +
  labs(x = "TCGA-AA-3527-01", y = "TCGA-D8-A1JK-01") +
  theme_minimal(base_size = 30)

# EXERCISE F: Sample Comparison
# What does this graph tell us about the similarity of gene expression
# between these two samples?
# Write your answer here:
#
#

# EXERCISE G: Improved Visualization
# Create a nicer version of this plot using ggplot().
# Make sure the axis labels are the sample names and consider adding a
# regression line with geom_smooth(method = "lm").
#
# Your code here:
#
#

################################################################################
#### SECTION 9: Exporting the Formatted Data ####
################################################################################

# Save the processed data for use in Part 2:
write.table(d, "Data/TCGA_RNA_seq_filtered_and_formatted.tsv",
            sep = "\t", col.names = TRUE, row.names = FALSE)


################################################################################
#### PART 2: Differential Expression Analysis ####
################################################################################

################################################################################
#### SECTION 10: Reload and Prepare Data ####
################################################################################

# Load the formatted data (if not continuing directly from Part 1):
d <- read_delim("Data/TCGA_RNA_seq_filtered_and_formatted.tsv")

# Encode Study and Gender as factors:
d$Study  <- as.factor(d$Study)
d$Gender <- as.factor(d$Gender)

################################################################################
#### SECTION 11: Testing a Single Gene ####
################################################################################

# We use linear models to find differentially expressed genes between the
# two cancer types.

# First, let's look at a single gene: BRCA1
tmp <- filter(d, name == "BRCA1")

# Visualize the expression difference:
ggplot(tmp, aes(x = Study, y = logCPM)) +
  geom_boxplot() +
  theme_minimal(base_size = 30)

# EXERCISE H: Visual Comparison
# Do you think the two groups show a significant difference?
# Write your answer here:
#
#

# Perform a t-test:
t_fit <- t.test(
  filter(tmp, Study == "BRCA")$logCPM,
  filter(tmp, Study != "BRCA")$logCPM
)
t_fit$p.value

# CHECKPOINT: We could also use a Wilcoxon test or a simple t-test.
# Using a linear model framework allows us to easily add covariates later.

################################################################################
#### SECTION 12: Linear Regression with Tidymodels ####
################################################################################

# Set up a linear regression workflow and fit it:
lm_form_fit <- linear_reg() %>%
  set_engine("lm") %>%
  fit(logCPM ~ Study, data = tmp)

fit1 <- tidy(lm_form_fit)
fit1

# EXERCISE I: Comparing Results
# Do you get the same result from the linear model as from the t-test?
# Write your answer here:
#
#

################################################################################
#### SECTION 13: Adding Covariates ####
################################################################################

# Let's add Age and Gender to the model formula:
lm_form_fit <- linear_reg() %>%
  set_engine("lm") %>%
  fit(logCPM ~ Age + Gender + Study, data = tmp)

fit2 <- tidy(lm_form_fit)
fit2

# EXERCISE J: Interpreting the Model
# How would you interpret fit2? Why has the p-value of Study changed?
# Write your answer here:
#
#

################################################################################
#### SECTION 14: Testing All Genes ####
################################################################################

# Now let's test all genes for differential expression.
# For each gene, we save the estimate, std.error, and p-value.

genes <- unique(d$name)
head(genes)

# Predefine model:
lm_form_fit <- linear_reg() %>%
  set_engine("lm")

# Loop over all genes and store the coefficients:
for (g in genes) {
  tmp <- filter(d, name == g)
  tmp_lm <- lm_form_fit %>%
    fit(logCPM ~ Study, data = tmp)

  # Extract coefficients:
  tmp_summary <- summary(tmp_lm$fit)
  coef      <- tmp_summary$coefficients[2]
  std_error <- tmp_summary$coefficients[4]
  p_value   <- tmp_summary$coefficients[8]

  # Store in data frame:
  tmp_out <- data.frame(
    Gene        = g,
    Coefficient = coef,
    Std_error   = std_error,
    p_value     = p_value
  )

  if (g == genes[1]) {
    out <- tmp_out
  } else {
    out <- bind_rows(out, tmp_out)
  }
}

fit1_all_genes <- out
fit1_all_genes <- arrange(fit1_all_genes, p_value)
head(fit1_all_genes)

# EXERCISE K: Most Significant Gene
# Which gene is the most significant when comparing the two cancer types?
# (HINT: Use slice_min() with p_value as the ordering variable)
# Write your answer here:
#
#

# EXERCISE L: Direction of Change
# Is this gene up- or down-regulated in the disease group?
# Write your answer here:
#
#

# EXERCISE M: Gene Annotation
# What is the common name of this gene? (Hint: look it up online)
# Write your answer here:
#
#

################################################################################
#### SECTION 15: Multiple Testing Correction ####
################################################################################

# EXERCISE N: Multiple Testing Problem
# What is the problem with performing so many statistical tests?
# Write your answer here:
#
#

# Account for multiple testing by calculating the false-discovery rate (FDR):
fit1_all_genes$q_value <- p.adjust(fit1_all_genes$p_value, method = "fdr")

# EXERCISE O: Top Hit Evaluation
# Do you think the top hit is interesting? Why or why not?
# Write your answer here:
#
#

################################################################################
#### SECTION 16: Adding Covariates to All Genes ####
################################################################################

# Now we repeat the loop, but also include Age and Gender in the model.
# Note: We adjust the coefficient extraction indices accordingly.

for (g in genes) {
  tmp <- filter(d, name == g)
  tmp_lm <- lm_form_fit %>%
    fit(logCPM ~ Age + Gender + Study, data = tmp)

  # Extract coefficients (Study is now the 4th term):
  tmp_summary <- summary(tmp_lm$fit)
  coef      <- tmp_summary$coefficients[4]
  std_error <- tmp_summary$coefficients[8]
  p_value   <- tmp_summary$coefficients[16]

  # Store in data frame:
  tmp_out <- data.frame(
    Gene        = g,
    Coefficient = coef,
    Std_error   = std_error,
    p_value     = p_value
  )

  if (g == genes[1]) {
    out <- tmp_out
  } else {
    out <- bind_rows(out, tmp_out)
  }
}

fit2_all_genes <- out
head(fit2_all_genes)

# EXERCISE P: Comparing Results
# Is the same gene still the most significant after accounting for Age and Gender?
# Write your answer here:
#
#

################################################################################
#### KEY TAKEAWAYS ####
################################################################################

# 1. Normalization (e.g., logCPM) is essential for comparing gene expression
#    across samples with different sequencing depths.
#
# 2. Linear models are flexible: you can add covariates to control for
#    confounders like age and gender.
#
# 3. Multiple testing correction (e.g., FDR) is crucial when testing thousands
#    of genes to avoid false positives.
#
# 4. Exploratory visualization helps you understand patterns before modeling.
