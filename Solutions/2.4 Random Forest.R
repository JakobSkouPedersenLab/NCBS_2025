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

################################################################################
#### PART 3: XGBoost Classification (Comparison) ####
################################################################################

# Load XGBoost engine via parsnip (comes with tidymodels)
# Build an XGBoost classifier comparable to RF setup

library(xgboost)

# Define XGBoost model
xgb_model <-
  boost_tree(
    trees = 2000,
    tree_depth = tune(),
    learn_rate = 0.05,
    loss_reduction = tune(),
    min_n = tune()
  ) %>%
  set_engine("xgboost") %>%
  set_mode("classification")

# Create a simple workflow using the same formula
xgb_wf <-
  workflow() %>%
  add_model(xgb_model) %>%
  add_formula(chdfate ~ . -id -followup)

# Set up cross-validation folds (stratified on outcome)
set.seed(222)
chd_folds <- vfold_cv(chd_train, v = 5, strata = chdfate)

# Define a small tuning grid to keep it fast and illustrative
xgb_grid <- grid_regular(
  tree_depth(range = c(3, 9)),
  loss_reduction(),
  min_n(range = c(2, 20)),
  levels = 4
)

# Use ROC AUC as the primary metric
xgb_metrics <- metric_set(roc_auc, accuracy, mcc, f_meas)

# Tune hyperparameters
set.seed(222)
xgb_tuned <- tune_grid(
  xgb_wf,
  resamples = chd_folds,
  grid = xgb_grid,
  metrics = xgb_metrics,
  control = control_grid(save_pred = TRUE)
)

# Select best by ROC AUC and finalize workflow
best_xgb <- select_best(xgb_tuned, metric = "roc_auc")
xgb_final_wf <- finalize_workflow(xgb_wf, best_xgb)

# Fit on full training data
xgb_fit <- fit(xgb_final_wf, data = chd_train)

# Evaluate on test data
chd_test_w_pred_xgb <- augment(xgb_fit, new_data = chd_test)

# ROC and AUC for XGBoost
xgb_roc <- chd_test_w_pred_xgb %>% roc_curve(truth = chdfate, .pred_TRUE)
autoplot(xgb_roc)

chd_test_w_pred_xgb %>% roc_auc(truth = chdfate, .pred_TRUE)

# Classification metrics side-by-side
classification_metrics(chd_test_w_pred_xgb, truth = chdfate, estimate = .pred_class, event_level = "second")

# Compare RF vs XGB ROC curves quickly by binding
bind_rows(
  xgb_roc %>% mutate(model = "XGBoost"),
  lr_auc %>% mutate(model = "RandomForest")
) %>%
  autoplot() + ggtitle("ROC: Random Forest vs XGBoost")

