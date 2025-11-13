########
# 
# Supervised learning II
# Regularization
#
# Simon Grund Sorensen, Jakob Skou Pedersen, Søren Besenbacher, Aarhus University
# 
#######

################################################################################
#### PART 1: Setup and Model Definition ####
################################################################################

## Load libraries ####
library(tidyverse)       # for tidyverse
library(tidymodels)      # for tidymodels
library(skimr)           # for variable summaries
library(vip)             # vizualise parameter importance

tidymodels_prefer() #Set tidymodels as the default whenever multiple packages have functions with the same name

#### PART 0: Load data and make recipe, as in 1_5_1
# load data
chd_full = read_rds("Data/chd_full.rds")
skim(chd_full)

#### PART 1 ####
# Now, let's preprocess our data for modelling
# First we split the data, so that we save some data for testing our model.

# Fix the random numbers by setting the seed. This enables the analysis to be reproducible when random numbers are used 
set.seed(222)

# Put 3/4 of the data into the training set, stratified by the outcome of interest
chd_split <- initial_split(chd_full, prop = 3/4., strata = chdfate)

# A)
# Create data frames for the two sets:
chd_train <- training(chd_split)
chd_test <- testing(chd_split)

# B)
# Create 10 cross-validation folds (data subsets)
set.seed(345)
folds <- vfold_cv(chd_train, v = 10)
folds

# initiate recipe
chd_rec <- 
  recipe(chdfate ~ ., data = chd_train)  %>% #Define the formula and the data
  update_role(id, new_role = "ID") %>% #Make it clear that "id" variable is not a predictor
  step_naomit(all_predictors(), all_outcomes()) %>% #Remove all rows with NA in the predictor or outcome variables
  step_dummy(all_nominal_predictors()) %>% #Convert characters and factors into binary columns
  step_zv(all_predictors()) %>% #Remove variable that have the same value across all rows
  step_center(all_numeric_predictors()) %>% #Normalize numeric data to mean of zero
  step_scale(all_numeric_predictors()) #normalize numeric data to standard deviation of one

chd_rec

################################################################################
#### PART 2: Hyperparameter Tuning ####
################################################################################

# tune penalty in lasso regression (mixture = 1)
lr_reg_mod <- 
  logistic_reg(penalty = tune(), mixture = 1) %>% 
  set_engine("glmnet")

# Make a workflow 
chd_wflow <- 
  workflow() %>%
  add_model(lr_reg_mod) %>%
  add_recipe(chd_rec)

# make grid of penalty values for the model optimization to go through
lr_reg_grid <- grid_regular(penalty(), levels = 50)
print(lr_reg_grid, n = 50)
plot(lr_reg_grid$penalty)

# fit across tuning grid of penalty values
# Note that the tuning uses the cross-validation folds to avoid overfitting.

lr_reg_res <- 
  chd_wflow %>% 
  tune_grid(folds,
            grid = lr_reg_grid,
            control = control_grid(save_pred = TRUE))
lr_reg_res

# A)
?collect_metrics()
# Use "collect_metrics()" to better understand the performance of the model
lr_reg_res %>% 
  collect_metrics()

# B)
# plot roc auc results 
# Try to read the plot-cde line by line and see if you have a basic understanding of
# what each line does
lr_plot <- 
  lr_reg_res %>% 
  collect_metrics() %>% 
  filter(.metric == "roc_auc") %>% 
  ggplot(aes(x = penalty, y = mean)) + 
  geom_point() + 
  geom_line() + 
  ylab("Area under the ROC Curve") +
  scale_x_log10(labels = scales::label_number())

lr_plot 

# C)
# What does the plot show us? (Hint, what happens when we increase the penalty)

# Under the hood, we have actually trained a whole lot of models, to evaluate the performance
# when including different variables.
# Let's select only the 15 best models
top_models <-
  lr_reg_res %>% 
  show_best(metric = "roc_auc", n = 15) %>% 
  arrange(penalty) 

top_models

# Or just the very  best model
# optimal tuning parameter
lr_reg_best <- lr_reg_res %>%
  select_best(metric = "roc_auc") #Finds the best tuning parameters to optimize model performance. 

################################################################################
#### PART 3: Final Model and Evaluation ####
################################################################################

# D) 
# What is the difference between "show_best()" and "select_best()"
# (Hint: What happens if you run "show_best()" with n = 1 ? )

# Now we have identified the optimal penalty. We can use that too
# generate  a single, optimal model:
final_chd_wflow <-
  chd_wflow%>% 
  finalize_workflow(lr_reg_best)

# Fit best model on train data
final_fit <- 
  final_chd_wflow %>%
  last_fit(chd_split)

# E)
# collect_metrics() the final fit, to evaluate performance
final_fit %>%
  collect_metrics()

# F
# plot roc usng autoplot
final_fit %>%
  collect_predictions() %>% 
  roc_curve(chdfate, .pred_TRUE) %>% 
  autoplot()

# G)
# view parameter values and vizualise importance
final_fit %>% 
  extract_fit_parsnip() %>% 
  tidy()

final_fit %>% 
  extract_fit_parsnip() %>% 
  vip(lambda = as.numeric(lr_reg_best[1,1]))

# What do you think of these feature importances? Is it meaningful or would you
#exclude any of the variables and re-run? 

#If time allows, try setting the mixture variable in the regressiom to different values between 0 and 1 and see what happens.
#Hint: Mixture = 0 is called Ridge regression, anything between 0 and 1 is called elastic net.s
