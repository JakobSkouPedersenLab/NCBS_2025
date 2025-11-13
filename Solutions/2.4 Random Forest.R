########
# 
# Supervised learning
# Random Forest
# 
# Simon Grund Sorensen, Jakob Skou Pedersen, Søren Besenbacher, Aarhus University
# 
#######

################################################################################
#### PART 1: Setup and Building Model ####
################################################################################

## Load libraries ####
library(tidyverse)       # for tidyverse
library(tidymodels)      # for tidymodels
library(skimr)           # for variable summaries
library(ranger)          # Allows tidymodels to use random forrest
tidymodels_prefer() #Set tidymodels as the default whenever multiple packages have functions with the same name

#### PART 0: Load data #####
# load data
chd_full = read_rds("Data/chd_full.rds")

##### PART 1:  Pre-process the data ####
# We wish to train a random forest (RF) model for classification, make predictions 
# on test set, and evaluate performance.

# The Random Forrest setup in the ranger package does not use recipes as we are used to,
# so we have to pre-process the data a bit ourselves. 
# Luckily, we only need to remove incomplete rows
chd_full = chd_full[complete.cases(chd_full),]

#Split into training and test data 
set.seed(222)
chd_split <- initial_split(chd_full, prop = 3/4., strata = chdfate)

# Create data frames for the two sets:
chd_train <- training(chd_split)
chd_test <- testing(chd_split)

#### PART 2: Make a random forrest ####
# Lets setup the random forest
rf_with_seed <- 
  rand_forest(trees = 2000, mtry = tune(), mode = "classification") %>%
  set_engine("ranger", seed = 63233)

rf_fit <- rf_with_seed %>% 
  set_args(mtry = 5) %>% 
  set_engine("ranger") %>%
  fit(chdfate ~ . -id -followup, data = chd_train)

rf_fit

# Can you interpret the rf_fit? It's a bit hard, because of the abstract
# nature of the forest. Let's make some performance evaluation to help us.

################################################################################
#### PART 2: Evaluation and Comparison ####
################################################################################

# add predictions to chd_test
chd_test_w_pred_rf <- augment(rf_fit, new_data = chd_test)

# Plot a ROC curve
lr_auc <- 
  chd_test_w_pred_rf  %>% roc_curve(truth = chdfate, .pred_TRUE) 
autoplot(lr_auc)

# Calculate ROC value
chd_test_w_pred_rf  %>% roc_auc(truth = chdfate, .pred_TRUE)

classification_metrics <- metric_set(accuracy, mcc, f_meas)
classification_metrics(chd_test_w_pred_rf, truth = chdfate, estimate = .pred_class, event_level = "second")

# A)
# The AUC using cross-validated logistic regression was 0.796
# iv. How did the performance change from logistic regression to random forest?
# v.  Why may it be relevant to evaluate different modeling procedures in 
#     different situations?

