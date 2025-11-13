########
# 
# Supervised learning and tidymodels (linear regression)
# 
# Simon Grund Sorensen, Jakob Skou Pedersen, Søren Besenbacher, Aarhus University
# 
#######

################################################################################
#### PART 1: Setup and Data Preparation ####
################################################################################

## Load packages
library(tidyverse)       # for tidyverse
library(tidymodels)      # for tidymodels
library(GGally)          # Extension to our plotting library(ggplot2)
library(skimr)           # for variable summaries
library(ggpubr)          # Package to make ggplots nicer
tidymodels_prefer(quiet = TRUE) #Set tidymodels as the default whenever multiple packages have functions with the same name

### PART 0: Loading and recoding the Framingham data set ####
# load and reformat the Framingham data set. We now look at a bigger chunk of the data than
# we have done in the last few modules. But, the features should look similar to you.
# Else, you can read about the data in "Framingham Didactic Data Dictionary.pdf". 
chd_full <- read_rds("Data/chd_full.rds")

### PART 1: Explore the data and its correlation structure ####
# look at the the data. What are some of the difference compared to before the formatting?
glimpse(chd_full)

# Try also using skim instead of glimpse: This gives a different overview.
chd_full %>% 
  skimr::skim() 

# A)
# Do you prefer glimpse or skim? Why?

# B)
# Check proportions of patients that developed congenital heart disorder (chdfate == T)
chd_full %>% 
  dplyr::count(chdfate == T) %>% 
  mutate(prop = n/sum(n))


# Time to prepare for modelling.
# First we split the data, so that we save some data for testing our model.
set.seed(122)
chd_split <- initial_split(data = chd_full, prop = 0.80)


# C)
# Create data frames for the two sets:
chd_train <- training(chd_split)
chd_test <- testing(chd_split)

# Look at the outputs
chd_test
chd_train
chd_split

# D)
# Consider the above and discuss with your neighbor(s):
# i.   What does the code do?
# ii.  Why is it important to split into a train and a test set?
# iii. Why would you want a large train set and why a large test set?
# iv.  Often a validation set is also created. Why and under what circumstances 
#      is that a good idea? 

## We start by doing some exploratory plotting of the data.
#  Generally follow the paradigm: "tegne før regne"[DK] 
#  i.e: "drawing before calculating"[EN]. (It rhymes in Danish..)

# D)
# why the `-1` below? 
var_names <- names(chd_train)[-1]

# F)
# now we make an all-against-all comparison of the variables using ggpairs.
chd_train %>%
  ggpairs(columns = var_names)

# G)
# Explore and understand the output you just generated
# i.   Press zoom and adjust the window size to see as many details as possible.
#      You may select a smaller set of variables to show, by adjusting the 
#      var_names vector, if the plot does not display well on your screen.
# ii.  First focus on the diagonal plots:
#      - What do they show? 
#      - What do you think determines the plot types?
#      - What types of plots are they?
# iii. What types of plots are present in the rest of plot and when?
# iv.  Is there any indication the explanatory variables correlate with the 
#      response variable chdfate?


################################################################################
#### PART 2: Building Models ####
################################################################################ 
# Now we try to fit a linear regression model to the data. Given that chdfate 
# is binary, it is not suitable for ordinary regression analysis. We therefore 
# start by focusing on bmi and ask whether its variation can be explained by 
# the other variables.  

# We start by ordinary linear regression with the default least squares 
# inference engine:
lm_fit <- 
  linear_reg() %>% 
  fit(bmi ~ sbp + sex, data = chd_train)

# print the fit:
lm_fit

# extract a tidy representation of the inference results:
tidy(lm_fit)

# A)
# i. Look at the results and try to understand them.
#    You may refer to the lm function or to Figure 7.4 in:
#    https://argoshare.is.ed.ac.uk/healthyr_book/regression.html 

# B)
# i. Try to modify the formula in the below (the part with "~") to specify 
#    a different model dependent on other or more variables. 
# ii. Look at and interpret the results.  
lm_fit <- 
  linear_reg() %>% 
  fit(bmi ~ sbp + sex + scl, data = chd_train)
tidy(lm_fit)

# C)
#      The below formula specifies that bmi should be explained by all variables 
#      in the data apart from "id".
# i.   Execute the code.
# ii.  Which variables have significant coefficients in the model?
# iii. Which variables have the largest effect sizes.
# iV.  How do you think scaling or nornmalisation would affect effect sizes? 
# We will look more at this in the next set of exercies (1_5_1), but 
# don't dive in too deep yet.

lm_fit <- 
  linear_reg() %>% 
  fit(bmi ~ . -id , data = chd_train)
tidy(lm_fit)

## draw dot-and-whisker plots (also known as a coefficient plot)
tidy_data = tidy(lm_fit)%>%
  filter(term != "(Intercept)")

ggplot(tidy_data, aes(x = estimate, y = term))+
  geom_errorbarh(aes(xmin = estimate-std.error, xmax = estimate+std.error))+
  geom_vline(xintercept = 0, lty = 2, linewidth = 0.5, col = "gray")+
  theme_classic()
  

################################################################################
#### PART 3: Testing and Evaluation ####
################################################################################

# D)
# Predict BMI for subjects in the TEST data and evaluate the model.
# D)
# What does the dot-and-whisker plt show and why is it useful?

# E)
# Make predictions on the test data (chd_test) using the above full model:
chd_test_w_pred <- augment(lm_fit, new_data = chd_test)

# F)
# look at the chd_test_w_pred  with glimpse, skim, or similar
chd_test_w_pred %>% skim()
chd_test_w_pred
# G)
# plot the original bmi values against the predictions by filling in 
# the ?s below:
ggplot(chd_test_w_pred ) + geom_point(aes(bmi, .pred))

# plot the residuals# plot the residuals.pred
ggplot(chd_test_w_pred ) + geom_histogram(aes(.resid))

# H)
# What does the distribution of residuals tell us about the model?

# Answer: it centered around zero, which suggest normally distributed 


# I)
# use the rmse function to calculate root-mean-squared-error (look up help page)
chd_test_w_pred  %>% rmse(truth = bmi, estimate = .pred)

#Can you formulate what that means in words?
