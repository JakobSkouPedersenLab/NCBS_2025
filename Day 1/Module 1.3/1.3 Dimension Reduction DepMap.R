#####
#
# Dimension Reduction in R - PCA and t-SNE
#
# Simon Grund Sorensen, Jakob Skou Pedersen, Søren Besenbacher, Farhad Zamani
# Aarhus University
# 
#####

################################################################################
#### PART 1: Principal Component Analysis (PCA) ####
################################################################################

################################################################################
#### SECTION 1: Introduction to Dimension Reduction ####
################################################################################

# Load required packages
library(tidyverse)  # For data manipulation and visualization
library(broom)      # For tidying model outputs
library(Rtsne)      # For t-SNE dimension reduction

# WHAT IS DIMENSION REDUCTION?
# When we have many variables (dimensions), it's hard to visualize patterns.
# Dimension reduction techniques compress many variables into 2-3 dimensions
# while preserving important patterns. This helps us:
# - Visualize high-dimensional data
# - Remove noise and redundancy
# - Identify clusters or groups
# - Prepare data for machine learning

# Two main techniques:
# 1. PCA (Principal Component Analysis) - linear method
# 2. t-SNE (t-Distributed Stochastic Neighbor Embedding) - non-linear method

################################################################################
#### SECTION 2: Loading and Preparing Data ####
################################################################################

# Load the CHD dataset
depmap_data <- read.csv("Data/depmap_export_2025-12-01 13_24_31.776975_subsetted.csv")
glimpse(depmap_data)

# Select column lineage_1 and all expression columns
depmap_data <- depmap_data %>% select(lineage_1, starts_with("Expression."))

# remove columns that are all NA
depmap_data <- depmap_data %>% select(where(~ !all(is.na(.))))

# Remove rows with missing values (required for PCA)
depmap_data <- depmap_data %>% filter(complete.cases(.))

# PCA works best with numeric variables. Let's select only numeric health metrics:
chd_numeric <- depmap_data %>% 
  select( -lineage_1)

# Verify we have only numeric data:
glimpse(chd_numeric)

################################################################################
#### SECTION 3: Principal Component Analysis (PCA) ####
################################################################################

# PCA finds new axes (principal components) that capture maximum variation.
# scale = TRUE standardizes variables so they're on the same scale.

pca_result <- prcomp(chd_numeric, scale = TRUE)

# Let's see what PCA found:
summary(pca_result)

# The augment() function adds PCA coordinates back to the original data
depmap_with_pca <- augment(pca_result, depmap_data)
glimpse(depmap_with_pca)

# Note: New columns .fittedPC1, .fittedPC2, etc. contain the PC coordinates

################################################################################
#### SECTION 4: Visualizing PCA Results ####
################################################################################

# EXERCISE A: Create a PCA plot showing PC1 vs PC2, colored by heart disease status

# First, let's look at what columns we have:
# names(depmap_with_pca)

# Make the plot (you fill in the blanks):
ggplot(depmap_with_pca, aes(x = .fittedPC1, y = .fittedPC2, color = lineage_1)) +
  geom_point() +
  labs(title = "PCA of DepMap Expression Data",
       x = "First Principal Component",
       y = "Second Principal Component")

# CHECKPOINT: Do you see any separation between different cancer lineages?

################################################################################
#### SECTION 5: Understanding Variance Explained ####
################################################################################

# How much information does each PC capture?
# We can calculate the proportion of variance explained by each component.

# Create a dataframe with variance information:
variance_explained <- data.frame(
  PC = seq_len(ncol(chd_numeric)),  # PC numbers: 1, 2, 3, 4, 5
  variance = pca_result$sdev^2      # Variance of each PC
)

# Calculate proportions:
variance_explained <- variance_explained %>%
  mutate(proportion_variance = variance / sum(variance))

# Look at the results:
variance_explained

# EXERCISE B: Read the code above and explain what each line does:
# Line 1 (data.frame): ?
# Line 2 (PC =): ?
# Line 3 (variance =): ?
# Line 4 (mutate): ?

# EXERCISE C: Create a bar plot showing variance explained by each PC
ggplot(variance_explained, aes(x = PC, y = proportion_variance)) +
  geom_bar(stat = "identity") +
  labs(title = "Variance Explained by Each Principal Component",
       y = "Proportion of Variance Explained")

# QUESTION: How many PCs do we need to capture most of the variation?

################################################################################
#### SECTION 6: Another Look at the PCA Plot ####
################################################################################

# EXERCISE D: Create another PCA plot (PC1 vs PC2) colored by lineage_1
# This is similar to Exercise A - good practice!
ggplot(depmap_with_pca, aes(x = ?, y = ?, color = ?)) +
  geom_point(size = 2, alpha = 0.6) +
  labs(title = "PCA: Looking for cancer lineage")

# EXERCISE E: Interpret the plot
# Does PCA separate patients by the lineage of their cancer? Does it work for all?
# Write your observations as a comment below:
# 


################################################################################
#### SECTION 7: Exporting PCA Results ####
################################################################################

# EXERCISE F: Save the data with PCA coordinates
# Replace the ? marks with: "\t" (tab separator) and TRUE
write.table(depmap_with_pca, "Data/depmap_w_pca_SAVED.tsv", 
            sep = ?, col.names = ?, row.names = FALSE)


################################################################################
#### PART 2: t-SNE (Non-Linear Dimension Reduction) ####
################################################################################

################################################################################
#### SECTION 8: t-SNE - Non-Linear Dimension Reduction ####
################################################################################

# t-SNE (t-Distributed Stochastic Neighbor Embedding) is a non-linear method
# that's particularly good at revealing clusters and local structure.
# Unlike PCA, t-SNE can capture complex, non-linear relationships.

# Prepare the numeric data for t-SNE:
tsne_input <- chd_numeric

# Run t-SNE:
# pca = FALSE means we're not pre-processing with PCA
# check_duplicates = FALSE speeds up computation
set.seed(42)  # For reproducibility
tsne_result <- Rtsne(tsne_input, pca = FALSE, check_duplicates = FALSE)

# EXERCISE G: The code above runs t-SNE. Run it and observe the output.

# EXERCISE H: Check the t-SNE help documentation
# Type: ?Rtsne
# Important parameters to consider:
# - perplexity: affects how many neighbors each point considers (default 30)
# - max_iter: number of iterations (default 1000)
# - pca: whether to pre-process with PCA (useful for large datasets)
# Question: Should we change any parameters? Why or why not?

# EXERCISE I: Explore the t-SNE output
# Click on 'tsne_result' in the Environment panel.
# Most of the output shows the parameters we used.
# The actual 2D coordinates are in: tsne_result$Y
# - tsne_result$Y[,1] = x-coordinates
# - tsne_result$Y[,2] = y-coordinates

################################################################################
#### SECTION 9: Visualizing t-SNE Results ####
################################################################################

# Combine t-SNE coordinates with the original data:
tsne_plot_data <- bind_cols(
  depmap_data,
  tsne_x = tsne_result$Y[, 1],
  tsne_y = tsne_result$Y[, 2]
)

# EXERCISE J: Create a t-SNE plot colored by lineage
# Replace ? with lineage_1
ggplot(tsne_plot_data, aes(x = tsne_x, y = tsne_y, color = ?)) +
  geom_point(size = 2, alpha = 0.6) +
  labs(title = "t-SNE of CHD Patient Data",
       x = "t-SNE Dimension 1",
       y = "t-SNE Dimension 2")

# CHECKPOINT: Compare this t-SNE plot with the PCA plot.
# Does t-SNE show better separation of disease groups?

################################################################################
#### SECTION 10: Advanced Exercises (Optional) ####
################################################################################

# EXERCISE K: Experiment with t-SNE parameters

# Try different perplexity values (typical range: 5-50):
# tsne_result2 <- Rtsne(tsne_input, perplexity = 10, check_duplicates = FALSE)

# Try running t-SNE on PCA results (recommended for large datasets):
# tsne_result3 <- Rtsne(pca_result$x, pca = TRUE, check_duplicates = FALSE)

# BONUS: Try with the full dataset (chd_full.rds) if you have time and computing power
# chd_full <- read_rds("Data/chd_full.rds")
# For large datasets, it's best to run t-SNE on PCA results to reduce computation time

################################################################################
#### KEY TAKEAWAYS ####
################################################################################

# PCA vs t-SNE:
# - PCA: Linear, fast, interpretable, good for initial exploration
# - t-SNE: Non-linear, slower, better at revealing clusters
# 
# When to use each:
# - PCA: Understanding which variables contribute most, quick exploration
# - t-SNE: Finding hidden clusters, visualizing complex patterns
# 
# Best practice: Often use PCA first, then t-SNE for deeper investigation
