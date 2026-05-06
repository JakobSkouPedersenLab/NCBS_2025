########
# 
# Supervised Learning II - Classification and Cross-Validation
#
# Simon Grund Sorensen, Jakob Skou Pedersen, Søren Besenbacher, Farhad Zamani
# Aarhus University
#
#######

################################################################################
#### PART 1: Data Preparation ####
################################################################################

################################################################################
#### SECTION 1: Introduction to Classification ####
################################################################################

# CLASSIFICATION vs REGRESSION:
# - Regression: Predict continuous values (e.g., BMI, blood pressure)
# - Classification: Predict categories (e.g., disease yes/no, species type)
#
# Today's goal: Predict whether a patient develops heart disease (TRUE/FALSE)
# We'll use LOGISTIC REGRESSION - designed for binary outcomes.

# Load libraries
library(tidyverse)   # Data manipulation and visualization
library(tidymodels)  # Modeling framework
library(skimr)       # Data summaries

tidymodels_prefer()  # Resolve function name conflicts

################################################################################
#### SECTION 2: Load and Examine Data ####
################################################################################

chd_full <- read_rds("Data/chd_full.rds")

# Examine the data:
skim(chd_full)

################################################################################
#### SECTION 3: Train-Test Split with Stratification ####
################################################################################

# For classification, we need to ensure balanced representation of both classes.
# STRATIFICATION ensures both train and test sets have similar proportions of TRUE/FALSE.

# Set seed for reproducibility:
set.seed(222)

# Split data: 75% training, 25% testing, stratified by chdfate
# Replace ? with 'chd_full'
chd_split <- initial_split(data = ?, prop = 3/4, strata = chdfate)

# EXERCISE A: Create training and testing datasets
chd_train <- training(chd_split)
chd_test <- testing(?)  # Replace ? with 'chd_split'

# Verify the split maintains class balance (replace ? with 'chdfate'):
chd_train %>% count(chdfate) %>% mutate(prop = n/sum(n))
chd_test %>% count(?) %>% mutate(prop = n/sum(n))

################################################################################
#### SECTION 4: Data Preprocessing with Recipes ####
################################################################################

# RECIPES are a powerful way to preprocess data consistently.
# Think of it like a cooking recipe - a series of steps applied in order.

# Create a preprocessing recipe:
chd_rec <- 
  recipe(chdfate ~ ., data = chd_train) %>%          # Define formula and data
  update_role(id, new_role = "ID") %>%               # Mark 'id' as ID, not predictor
  step_naomit(all_predictors(), all_outcomes()) %>%  # Remove rows with missing values
  step_dummy(all_nominal_predictors()) %>%           # Convert categorical to binary (0/1)
  step_zv(all_predictors()) %>%                      # Remove zero-variance predictors
  step_center(all_numeric_predictors()) %>%          # Center: subtract mean (mean = 0)
  step_scale(all_numeric_predictors())               # Scale: divide by SD (SD = 1)

# WHY THESE STEPS?
# - step_dummy: Converts "Male"/"Female" to 0/1 (models need numbers)
# - step_zv: Removes constant variables (no information)
# - step_center/scale: Puts all variables on same scale (important for logistic regression!)

# EXERCISE B: Examine the recipe
chd_rec

# Are any operations unclear? Key points:
# - Recipes define steps but don't execute them yet
# - We'll "train" the recipe on training data
# - Then apply the SAME transformations to test data

# Train the recipe by learning parameters from training data:
# (e.g., learning the mean and SD for centering/scaling)
preprocessed_rec <- prep(chd_rec, training = chd_train) 
preprocessed_rec

# Apply the recipe to training data ("baking"):
chd_train_baked <- bake(preprocessed_rec, new_data = chd_train)
glimpse(chd_train_baked)

# CRITICAL: Apply the SAME recipe to test data
# (using training means/SDs, NOT test data statistics!)
chd_test_baked <- bake(preprocessed_rec, new_data = chd_test)
glimpse(chd_test_baked)

# CHECKPOINT: Notice how categorical variables like 'sex' became 'sex_Male'
# This is one-hot encoding from step_dummy()


################################################################################
#### PART 2: Model Building and Training ####
################################################################################

################################################################################
#### SECTION 5: Building a Logistic Regression Model ####
################################################################################

# LOGISTIC REGRESSION predicts probabilities for binary outcomes (0 to 1).
# Unlike linear regression, it uses a logistic function to constrain predictions.

# Define the model:
lr_model <- 
  logistic_reg() %>% 
  set_engine("glm")  # General Linear Model engine

################################################################################
#### SECTION 6: Creating a Workflow ####
################################################################################

# WORKFLOWS combine preprocessing (recipe) and modeling in one object.
# Benefits:
# - Ensures preprocessing is applied consistently
# - Makes code cleaner and less error-prone
# - Easier to deploy in production

chd_workflow <- 
  workflow() %>% 
  add_model(lr_model) %>% 
  add_recipe(?)  # We use the un-baked recipe object from above

chd_workflow

# Notice: We pass the raw recipe, not the baked data!
# The workflow will handle preprocessing automatically.

################################################################################
#### SECTION 7: Training the Model ####
################################################################################

# Fit the model using the workflow:
# This will automatically prep the recipe and bake the data
chd_fit <- 
  chd_workflow %>% 
  fit(data = chd_train)

chd_fit

# Extract coefficients in a tidy format:
model_coefs <- chd_fit %>% 
  extract_fit_parsnip() %>% 
  tidy()

model_coefs

# EXERCISE A: Identify the most significant variables
# Look at the p.value column. Which three variables have the smallest p-values?
# These are the strongest predictors of heart disease.
# Your answer:
# 1. 
# 2. 
# 3. 


################################################################################
#### SECTION 8: Making Predictions ####
################################################################################

# EXERCISE B: Get predicted probabilities on training data
fitted_train <- augment(chd_fit, chd_train)

glimpse(fitted_train)
# Notice the new columns:
# - .pred_class: Predicted category (TRUE/FALSE)
# - .pred_TRUE: Probability of heart disease
# - .pred_FALSE: Probability of no heart disease

# EXERCISE C: Visualize predicted probabilities
# Replace ? with '.pred_TRUE'
ggplot(data = fitted_train) + 
  geom_histogram(aes(x = ?), bins = 30, fill = "steelblue", alpha = 0.7) +
  labs(title = "Distribution of Predicted Probabilities (Training Data)",
       x = "Predicted Probability of Heart Disease",
       y = "Count")

# EXERCISE D: Predict on test data
# Replace ? with 'chd_test'
predict(chd_fit, new_data = ?)

# To get probabilities instead of classes, add type = "prob":
predict(chd_fit, new_data = chd_test, type = "prob")

################################################################################
#### PART 3: Model Evaluation ####
################################################################################

################################################################################
#### SECTION 9: Evaluating Predictions on Test Data ####
################################################################################

# EXERCISE E: Add predictions to test data
chd_aug <- augment(chd_fit, chd_test)
glimpse(chd_aug)

# EXERCISE F (BONUS): Create a boxplot comparing predicted probabilities
# by actual outcome
ggplot(chd_aug, aes(x = chdfate, y = .pred_TRUE, fill = chdfate)) +
  geom_boxplot() +
  labs(title = "Predicted Probabilities by Actual Outcome",
       x = "Actual Heart Disease Status",
       y = "Predicted Probability of Heart Disease")

# EXERCISE G: Assess separation ability
# Question: Can our model separate chdfate TRUE from FALSE?
# Look at the boxplot:
# - Do the boxes overlap a lot or are they well separated?
# - What would perfect separation look like?
# Your answer:
#


################################################################################
#### SECTION 10: ROC Curves and AUC ####
################################################################################

# ROC (Receiver Operating Characteristic) Curve:
# - Shows trade-off between True Positive Rate and False Positive Rate
# - X-axis: False Positive Rate (wrongly predicting disease)
# - Y-axis: True Positive Rate (correctly predicting disease)
# - Perfect classifier: goes straight up then right (area = 1.0)
# - Random classifier: diagonal line (area = 0.5)

# Plot ROC curve:
chd_aug %>% 
  roc_curve(truth = chdfate, .pred_TRUE) %>% 
  autoplot() +
  labs(title = "ROC Curve for Heart Disease Prediction")

# Calculate AUC (Area Under Curve):
auc_result <- chd_aug %>% 
  roc_auc(truth = chdfate, .pred_TRUE)

auc_result

# EXERCISE A: Understanding AUC
# If the model had NO ability to separate classes, what would AUC be?
# Answer: 0.5 (same as random guessing)
# 
# AUC interpretation:
# - 0.5: No better than random
# - 0.7-0.8: Acceptable
# - 0.8-0.9: Excellent
# - 0.9+: Outstanding

# EXERCISE B: Evaluate your model's AUC
# Look at the .estimate value. Is this a good score?
# Your assessment:
# 

# EXERCISE C: Brainstorm improvements
# How could we improve the model? Consider:
# - Adding more features
# - Feature engineering (interactions, polynomials)
# - Different model types (we'll explore this in later modules!)
# - More data
# Your ideas:
#  


################################################################################
#### PART 4: Cross-Validation ####
################################################################################

################################################################################
#### SECTION 11: Cross-Validation ####
################################################################################

# CROSS-VALIDATION provides more robust performance estimates than a single test set.
# 
# K-FOLD CROSS-VALIDATION:
# 1. Split training data into K "folds" (subsets)
# 2. Train on K-1 folds, validate on the remaining fold
# 3. Repeat K times, each fold getting a turn as validation set
# 4. Average the K performance scores
#
# Benefits:
# - Uses all data for both training and validation
# - More stable performance estimates
# - Detects overfitting

# Create 10-fold cross-validation splits:
set.seed(345)
cv_folds <- vfold_cv(chd_train, v = 10)
cv_folds

# Each row is one fold showing how many observations are in training/validation

################################################################################
#### SECTION 12: Fit Model with Cross-Validation ####
################################################################################

set.seed(456)
chd_fit_cv <- 
  chd_workflow %>% 
  fit_resamples(cv_folds)

chd_fit_cv

# Collect and view metrics from all folds:
cv_metrics <- collect_metrics(chd_fit_cv)
cv_metrics

# These are AVERAGED across all 10 folds!

################################################################################
#### SECTION 13: Comparing Cross-Validation vs Test Set ####
################################################################################

# Cross-validation performance (on training data folds):
cv_metrics

# Test set performance:
test_auc <- chd_aug %>% 
  roc_auc(truth = chdfate, .pred_TRUE)

test_accuracy <- chd_aug %>% 
  accuracy(truth = chdfate, .pred_class)

test_auc
test_accuracy

# EXERCISE A: Compare the results
# Are the cross-validation and test set metrics similar?
# Your observation:
# 

# EXERCISE B: What if they were very different?
# Possible interpretations:
# - CV much better than test: May have gotten "lucky" or "unlucky" with test split
# - Test much better than CV: Possible overfitting, or test set easier
# - Large discrepancy: Could indicate problem with data split or model instability
#
# If you see big differences, what would you do?
# - Try different random seeds
# - Use stratification (we already did!)
# - Increase test set size
# - Check for data leakage

################################################################################
#### KEY TAKEAWAYS ####
################################################################################

# 1. Classification predicts categories, not continuous values
# 2. Logistic regression outputs class probabilities (0 to 1)
# 3. Recipes ensure consistent preprocessing
# 4. Workflows combine preprocessing and modeling elegantly
# 5. ROC curves and AUC measure classification performance
# 6. Cross-validation gives more reliable performance estimates
# 7. Always compare training/CV/test performance to detect overfitting

# NEXT: We'll explore regularization to improve model generalization!
