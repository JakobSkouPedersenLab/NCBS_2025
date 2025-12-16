########
# 
# Analysis of bulk RNA seq data 
#
# Simon Grund Sorensen, Jakob Skou Pedersen, Søren Besenbacher, Farhad Zamani
# Aarhus University
# 
#######

#Load packages
library(tidyverse)
library(tidymodels)

##### PART 1: Load data and do some visualisations ####
# Load counts
d <- readRDS("Data/TCGA_Formatted_data.rds")
# = filter(d, Cancertype %in% c("BRCA", "COAD")) #Only keep two cancer types for speed
#d = dplyr::select(d, -c(Cancertype)) #Remove cancer type column for now


d <- read_delim("Data/RNA_seq_counts.gtf", delim ="\t", col_names = T)

# Check it out 
glimpse(d)

# A) 
# What is in rows? What is in columns? 


# Reformat to long format and count the number of samples that had any expression of each gene
Summary_d <-
  d%>%
  pivot_longer(cols = -gene_id) %>%
  group_by(gene_id) %>%
  mutate(n_above_0 = sum(value > 0)) %>%
  distinct(gene_id, n_above_0) %>%
  ungroup()

# Comment on each line above what it does,either in writing here or discussing it with someone

# Visualize
hist(Summary_d$n_above_0)


# B) 
# How many genes have 0 expression across all samples?

# Let's only look at  genes that have an expression in at least 25% of the samples:
genes_to_include = filter(Summary_d, n_above_0 > 0.25*ncol(d))$gene_id
d = filter(d, gene_id %in% genes_to_include)


# C) 
# How many genes are we looking at now?


# To save space on the RAM, we will subset down to a random set of 200 genes
set.seed(139)
d = slice_sample(d, n = 200)

# Now, let's perform normalization and transformation of the data (logCPM)
# For each sample, we calculate the sum of reads in millions. 
# To do this, we change the form of the data to 'long-format'.
# Try to rread through each line of the code below and see if you can understand what it does.
# Feel free to add comments by adding a hashtag at the end of each line, followed by your text
d = 
  d %>%
  pivot_longer(cols = -gene_id) %>%
  group_by(name) %>%
  mutate(
    sample_summed_counts_in_million = sum(value) / 1e6
  ) %>%
  ungroup() %>%
  mutate(
    CPM = value+1/sample_summed_counts_in_million, #Add a pseudocount of one because we cannt take the log of 0.
    logCPM = log(CPM)
  )

d$sample_summed_counts_in_million = NULL #We remove this column simply to clean up a bit

# Let's quickly look at the log-transofed counts-per-million
hist(d$logCPM)

# D)
# Why is their a peak around 2 which seems to stand out from the rest?

# Let's add annotation of whether each sample came from a patient with tremor or not
# To do this, we load the metadata
meta = read_delim("Data/Metadata.tsv", delim = "\t")
glimpse(meta)

d = left_join(d, meta, by = c("name" = "ID"))

# Make a boxplot of the logCPM for the two different diagnosis groups
ggplot(d, aes(x = diagnosis, y = logCPM)) + geom_boxplot()

# E) What does this graph tell us?

# Let's plot just two samples
s1 = filter(d, name == "SRR9835803")
s2 = filter(d, name == "SRR9835947")
plot(s1$logCPM, s2$logCPM)

# F) 
# What does this graph tell us?

# G) 
# Can you make a nicer version of this graph in ggplot? 
# Make sure that axis labels are the sample names and consider adding a 
# line using geom_smooth(method = "lm")

# Let's export the formatted data so we can use it in part 2
write.table(d, "Data/RNA_seq_filtered_and_formatted.tsv", sep ="\t", col.names = T, row.names = F)

### PART 2 ####

# Load the data (if you didnt go directly from PART 1)
d <- read_delim("Data/RNA_seq_filtered_and_formatted.tsv")

# Fix the encoding of diagnosis and gender
d$diagnosis = as.factor(d$diagnosis) #It's not strictly necessary, but good to ensure that R sees categorical values as factors
d$gender = as.factor(d$gender)

# Now, we will use linear models to find 
# differentially expressed genes between essential tremor samples and controls

# First, we compare tremor and control for a single gene
tmp <- filter(d, gene_id == "ENSG00000158887.15")

# Let's see if we can visually see a difference
ggplot(tmp, aes(x = diagnosis, y = logCPM))+
  geom_boxplot()

# A) 
#Do you think the two boxplots are significantly shifted?

# Let's test it using a t-test
t_fit <- t.test(filter(tmp, diagnosis == "Control")$logCPM, 
               filter(tmp, diagnosis != "Control")$logCPM )
t_fit$p.value

# In this case we could also simply to a wilcoxon test or a t-test as we just did. 
# But, by using modelling we are ready for also considering age and gender later.

# Setup a workflow. and fit it immediately
lm_form_fit <- 
  linear_reg() %>% 
  set_engine("lm")%>% 
  fit(logCPM ~ diagnosis, data = tmp)

fit1 <- tidy(lm_form_fit)
fit1


# B)
# Do we get the same result from the linear model as from the t-test?


# C)
#Let's add age and gender to the model formula

lm_form_fit <- 
  linear_reg() %>% 
  set_engine("lm")%>% 
  fit(logCPM ~ age+gender+diagnosis, data = tmp)

fit2 = tidy(lm_form_fit)
fit2


# D) How would you interpret fit2, and why has the p-value of diagnosis changed?

# Now, let's test across all the genes. 
# For each gene, we save the estimate, std.error, and p-value
genes = unique(d$gene_id) #List of all genes
head(genes)

#Predefine model:
lm_form_fit <- 
  linear_reg() %>% 
  set_engine("lm")

# Loop over all genes and store the model coefficients for each
for(g in genes){
  tmp = filter(d, gene_id == g)
  tmp_lm = lm_form_fit %>%
    fit(logCPM ~ diagnosis, data = tmp)
  
  # Extract coefficients (could also be done using the tidy() function)
  tmp_summary = summary(tmp_lm$fit)
  coef = tmp_summary$coefficients[2]
  std_error = tmp_summary$coefficients[4]
  p_value = tmp_summary$coefficients[8]
  
  # Store in data frame
  tmp_out = data.frame(
    Gene = g,
    Coefficient = coef,
    Std_error = std_error,
    p_value = p_value
  )
  if(g == genes[1]){ #True for first gene
    out = tmp_out
  }else{ #True for all other genes
    out = bind_rows(out, tmp_out)
  }
}

fit1_all_genes = out #Save the output in a new variable
head(fit1_all_genes)

# E)
# Which gene is the most significant when comparing disease to control? 
# (HINT: Use slice_min with p_value being the variable to order by)


# F)
# Is this gene up- or down-regulated in the disease group? 


# G)
# What is the normal name of this gene? (hint: Google)


# H)
# What is the problem with the many tests we have made?

# Let's account for multiple testing by calculating the false-discovery-rate (FDR)
fit1_all_genes$q_value = p.adjust(fit1_all_genes$p_value, method = "fdr")


# I) 
# Do you think the top-hit is interesting? Why / why not?


# We now change the loop above so that we also consider the effect of age and gender. 
# Note that we had to adjust coefficient extraction.

for(g in genes){
  tmp = filter(d, gene_id == g)
  tmp_lm = lm_form_fit%>%fit(logCPM ~ age+gender+diagnosis, data = tmp)
  
  # Extract coefficients
  tmp_summary = summary(tmp_lm$fit)
  coef = tmp_summary$coefficients[4]
  std_error = tmp_summary$coefficients[8]
  p_value = tmp_summary$coefficients[16]
  
  # Store in data frame
  tmp_out = data.frame(
    Gene = g,
    Coefficient = coef,
    Std_error = std_error,
    p_value = p_value
  )
  if(g == genes[1]){ #True for first gene
    out = tmp_out
  }else{ #True for all other genes
    out = bind_rows(out, tmp_out)
  }
}

fit2_all_genes = out
head(fit2_all_genes)


#J) Is the same gene still the most significant?
