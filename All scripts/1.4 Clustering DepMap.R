#####
# 
# Clustering in R - K-Means and Hierarchical Clustering
# 
# Simon Grund Sorensen, Jakob Skou Pedersen, Søren Besenbacher, Farhad Zamani
# Aarhus University
# 
#####

################################################################################
#### PART 1: K-Means Clustering ####
################################################################################

################################################################################
#### SECTION 1: Introduction to Clustering ####
################################################################################

# Load required packages
library(tidyverse)  # For data manipulation and visualization
library(ggdendro)   # For pretty dendrograms

# WHAT IS CLUSTERING?
# Clustering groups similar data points together without prior labels.
# It's an "unsupervised" learning method - we don't tell the algorithm
# which group each point belongs to; it discovers patterns on its own.

# Two main approaches:
# 1. K-means: Partition data into K clusters (you choose K)
# 2. Hierarchical: Build a tree showing how observations group at different scales

################################################################################
#### SECTION 2: Loading Data with PCA Coordinates ####
################################################################################

# Load the data with PCA coordinates from the previous module:
depmap_pca <- read_tsv("Data/depmap_w_pca_SAVED.tsv")

# Examine the data:
glimpse(depmap_pca)

# We'll cluster patients based on their PCA coordinates rather than
# the original variables. This reduces noise and speeds up computation.

################################################################################
#### SECTION 3: K-Means Clustering ####
################################################################################

# K-means divides data into K clusters by:
#  1. Randomly assigning observations to the K cluster classes
#  2. Calculate the centroid (center) for each class
#  3. Reassign each observation to the class with the nearest centroid
#  4. Repeat steps 2 & 3 until convergence (no more changes in cluster class assignments)

# First, extract only the PCA columns usign select and starts_with():
pca_coords <- select(depmap_pca, starts_with(".fittedPC"))

# EXERCISE A: Perform k-means clustering with 4 clusters
# Replace the first ? with 'pca_coords' and the second ? with 4
kmeans_result <- kmeans(x = ?, centers = ?)

# See what's in the result:
names(kmeans_result)

# The cluster assignments are in: kmeans_result$cluster

################################################################################
#### SECTION 4: Visualizing K-Means Clusters ####
################################################################################

# Add cluster assignments to our data:
# We convert to factor so ggplot treats it as categorical (for colors)
depmap_pca$cluster_kmeans <- as.factor(kmeans_result$cluster)

# EXERCISE B: Plot the clusters in PCA space
# Replace ? with 'cluster_kmeans'
ggplot(depmap_pca, aes(x = .fittedPC1, y = .fittedPC2)) +
  geom_point(aes(color = ?), size = 2) +
  labs(title = "K-Means Clustering (K=4)",
       x = "PC1", y = "PC2")

# QUESTION: Try commenting out the as.factor() line above and re-run.
# What happens to the plot colors? Why?


################################################################################
#### SECTION 5: Choosing the Right Number of Clusters ####
################################################################################

# EXERCISE C: Experiment with different numbers of clusters

# Try K = 2:
kmeans_2 <- kmeans(pca_coords, centers = 2)
depmap_pca$cluster_2 <- as.factor(kmeans_2$cluster)

ggplot(depmap_pca, aes(x = .fittedPC1, y = .fittedPC2)) +
  geom_point(aes(color = cluster_2), size = 2) +
  labs(title = "K-Means with K=2")

# Try K = 6:
kmeans_6 <- kmeans(pca_coords, centers = 6)
depmap_pca$cluster_6 <- as.factor(kmeans_6$cluster)

ggplot(depmap_pca, aes(x = .fittedPC1, y = .fittedPC2)) +
  geom_point(aes(color = cluster_6), size = 2) +
  labs(title = "K-Means with K=6")

# EXERCISE D: Reflect on cluster numbers
# How many clusters are too many? Too few? What looks most meaningful?
# Consider: Do the clusters correspond to real biological groups?
# Write your answer as a comment below:
# 


################################################################################
#### PART 2: Hierarchical Clustering ####
################################################################################

################################################################################
#### SECTION 6: Hierarchical Clustering ####
################################################################################

# Hierarchical clustering builds a tree (dendrogram) showing how
# samples group together at different levels of similarity.
# Unlike k-means, you don't need to specify the number of clusters upfront!

# Step 1: Create a distance matrix
# This calculates how "different" each patient is from every other patient
distance_matrix <- dist(pca_coords)

# Step 2: Perform hierarchical clustering
# This builds the tree structure
hclust_result <- hclust(d = distance_matrix)

# View summary:
hclust_result

################################################################################
#### SECTION 7: Visualizing Hierarchical Clusters ####
################################################################################

# Create a dendrogram (tree diagram):
ggdendrogram(hclust_result) +
  labs(title = "Hierarchical Clustering Dendrogram")

# INTERPRETING THE DENDROGRAM:
# - The y-axis shows the distance at which clusters merge
# - Each leaf (bottom) is one cell line
# - Height of branches shows how similar clusters are
# - Cutting the tree at different heights gives different numbers of clusters

################################################################################
#### SECTION 7b: Color Leaves by Cancer Lineage ####
################################################################################

# EXERCISE: "Let's see if the dendrogram separates the cancer lineages"
# We'll color the leaf labels by `lineage_1` from `depmap_pca`.
# To ensure labels match, set rownames of the PCA matrix to the `.rownames` column.

# 1) Make sure distance/cluster labels correspond to `.rownames`
#    (Run before dist/hclust if starting fresh; safe to set here too.)

#rownames(pca_coords) <- depmap_pca$.rownames


# 2) Build dendrogram data and join lineage info for labels
library(ggdendro)
library(dplyr)

# Extract the dendrogram data so we can plot in ggplot
ddata <- ggdendro::dendro_data(hclust_result, type = "rectangle")

# Ensure compatible join key types: coerce `.rownames` in depmap_pca to character
depmap_pca <- depmap_pca %>% mutate(.rownames = as.character(.rownames))

# Add lineage to the dendrogram labels
label_df <- ddata$labels %>%
  rename(.rownames = label) %>%
  left_join(depmap_pca %>% select(.rownames, lineage_1), by = ".rownames")

# 3) Plot: segments in grey, labels colored by lineage  — For now just run the ggplot, 
# but if you have the time, go through it and undertand each line.
ggplot() +
  geom_segment(data = ddata$segments,
               aes(x = x, y = y, xend = xend, yend = yend),
               color = "grey50", linewidth = 0.4) +
  geom_point(data = label_df, aes(x = x, y = 0, color = lineage_1), size = 1) +
  scale_color_brewer(palette = "Set2") +
  labs(title = "Dendrogram (Lineage Tick Marks)", color = "Lineage") +
  theme_minimal() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid = element_blank())

################################################################################
#### SECTION 8: Cutting the Dendrogram ####
################################################################################

# We can "cut" the tree at different heights to get clusters:

# Cut at height 30 (few, large clusters):
clusters_h40 <- cutree(hclust_result, h = 40)
table(clusters_h40)  # How many cell lines in each cluster?

# Cut at height 30:
clusters_h30 <- cutree(hclust_result, h = 30)
table(clusters_h30)

# Cut at height 20 (many, small clusters):
clusters_h20 <- cutree(hclust_result, h = 20)
table(clusters_h20)

# EXERCISE: Add cutlines to the dendrogram to visualize different cuts
ggdendrogram(hclust_result) +
  geom_hline(yintercept = 40, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 30, linetype = "dashed", color = "blue") +
  geom_hline(yintercept = 20, linetype = "dashed", color = "green") +
  labs(title = "Dendrogram with Cut Heights")

# QUESTION: What is a good place to cut the tree? Why?
# Consider: Which height gives meaningful, interpretable groups?
# Bonus challenge: Could you add the same cutlines to the ggplot version?


################################################################################
#### SECTION 9: Comparing Methods ####
################################################################################

# Let's compare k-means and hierarchical clustering:
depmap_pca$cluster_hclust <- as.factor(cutree(hclust_result, k = 4))

# Plot both side by side:
library(patchwork)  # For combining plots

p1 <- ggplot(depmap_pca, aes(.fittedPC1, .fittedPC2, color = cluster_kmeans)) +
  geom_point() + labs(title = "K-Means (K=4)")

p2 <- ggplot(depmap_pca, aes(.fittedPC1, .fittedPC2, color = cluster_hclust)) +
  geom_point() + labs(title = "Hierarchical (K=4)")

p1 + p2

################################################################################
#### KEY TAKEAWAYS ####
################################################################################

# K-Means:
# - Fast and simple
# - Need to specify K upfront
# - Sensitive to initial random placement
# - Works well for spherical clusters

# Hierarchical:
# - No need to specify K upfront
# - Shows relationships at multiple scales
# - Slower for large datasets
# - Creates a tree structure you can explore

# Choosing between them:
# - Try both and compare results!
# - K-means for large datasets or when you know K
# - Hierarchical for exploring structure or when K is unknown
