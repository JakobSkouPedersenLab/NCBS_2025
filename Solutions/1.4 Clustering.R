#####
# 
# 
# Simon Grund Sorensen, Jakob Skou Pedersen, Søren Besenbacher, Aarhus University
# 
####

################################################################################
#### PART 1: K-Means Clustering ####
################################################################################

#Load packages
library(tidyverse)
library(ggdendro)

#Load the data including PCs from the last module
d = read_delim("Data/d_w_pca_SAVED.tsv")
glimpse(d)

# A) 
# Perform k-means clustering on the PCs, with four centers.
PCs = select(d, .fittedPC1, .fittedPC2, .fittedPC3, .fittedPC4, .fittedPC5)

?kmeans
kk = kmeans(x = PCs, centers = 5) #Use ?help(kmeans) if you need help on how to replace the question marks

# B)
# Plot the inferred clusters with different colors on the PC1 vs PC2 plot
# We set clusters as factors instead of numbers, 
#for the graph to look nice. What happens if you hashtag this line and runs the script?
d$Cluster = as.factor(kk$cluster) 

ggplot(d, aes(x = .fittedPC1, y = .fittedPC2))+
  geom_point(aes(color = Cluster))

# C)
# Try running kmeans and visualising with different numbers of clusters (centers). 
kk = kmeans(x = PCs, centers = 2)

d$Cluster = as.factor(kk$cluster) 

ggplot(d, aes(x = .fittedPC1, y = .fittedPC2))+
  geom_point(aes(color = Cluster, shape = ))

################################################################################
#### PART 2: Hierarchical Clustering ####
################################################################################

# D)
# Perform hierarchical clustering with base R 
# instead of k-means (hclust)

# D)
# How many clusters are too many? How many are too few?

# PART 2: Hierarchical clustering

#First we make a distance matrix of the principal components
dist_mat_PC = dist(PCs)

#Then we generate hierarchical clustering
clusters = hclust(d = dist_mat_PC)

#Use the ggdendro to make a nice ggplot of the clusters
ggdendrogram(clusters) +
  geom_hline(yintercept = 10, linetype = "dashed", color = "red")

# Interpret the Dendrogram. How many clusters are there if we cut it at y=11? y=7.5? y=5? What is a good
# place to cut the tree do you think? Try to visualize it by making a horizontal line in the plot

y = ?
  
ggdendrogram(clusters) +
  geom_?(yintercept = y, linetype = "dashed", color = "red")
