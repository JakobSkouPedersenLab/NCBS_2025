########
# 
# Supervised Learning I - Linear Regression with Tidymodels
# 
# Simon Grund Sorensen, Jakob Skou Pedersen, Søren Besenbacher, Farhad Zamani
# Aarhus University
# 
########

################################################################################
#### PART 1: Setup and Data Preparation ####
################################################################################

################################################################################
#### SECTION 1: Introduction to Supervised Learning ####
################################################################################

# WHAT IS SUPERVISED LEARNING?
# Supervised learning uses labeled data to predict outcomes for new data.
# We "supervise" the algorithm by showing it examples with known answers.
# 
# Key concepts:
# - Features (predictors): Variables we use to make predictions (age, BMI, etc.)
# - Target (response): Variable we want to predict (e.g., BMI, disease status)
# - Training data: Data used to build the model
# - Test data: Data used to evaluate the model
#
# Today we'll use LINEAR REGRESSION to predict continuous outcomes.

################################################################################
#### SECTION 2: Loading Packages ####
################################################################################

library(tidyverse)   # For data manipulation and visualization
library(tidymodels)  # For modeling workflows
library(GGally)      # For pairwise plotting
library(skimr)       # For data summaries
library(ggpubr)      # For enhanced ggplot themes

# Set tidymodels as default when functions have the same name across packages:
tidymodels_prefer(quiet = TRUE)

################################################################################
#### SECTION 3: Loading the Full Dataset ####
################################################################################

# Load the complete Framingham dataset (larger than previous exercises):
chd_full <- read_rds("Data/chd_full.rds")

# Examine the structure:
glimpse(chd_full)

# Try skim() for a more detailed overview:
chd_full %>% 
  skimr::skim()

# EXERCISE A: Compare glimpse() and skim()
# Which do you prefer? Why? Consider:
# - Speed and readability
# - Amount of information provided
# - When you'd use each one
# Your answer:
# 


################################################################################
#### SECTION 4: Understanding the Target Variable ####
################################################################################

# Before modeling, always examine your target variable!

# EXERCISE B: Check proportions of patients with heart disease
# Replace ? with 'chdfate'
chd_full %>% 
  dplyr::count(?) %>% 
  mutate(proportion = n / sum(n))

# CHECKPOINT: Is the dataset balanced? Or do we have many more of one class?

################################################################################
#### SECTION 5: Train-Test Split ####
################################################################################

# CRITICAL CONCEPT: We NEVER test our model on the same data we trained it on!
# This would be like teaching to the test - we'd get overly optimistic results.

# Set a random seed for reproducibility:
set.seed(122)

# Split data: 80% for training, 20% for testing
# Replace ? with 'chd_full'
chd_split <- initial_split(data = ?, prop = 0.80)

# EXERCISE C: Create training and testing datasets
chd_train <- training(chd_split)
chd_test <- testing(chd_split)  # Replace ? with 'chd_split'

# Examine the split:
chd_train  # Training data
chd_test   # Testing data
chd_split  # Split information

# EXERCISE D: Discussion questions (discuss with neighbors or write answers):
# i.   What does initial_split() do? How does it decide which rows go where?
# 
# ii.  Why is it crucial to split into train and test sets?
# 
# iii. Trade-off: Why might you want a LARGE training set? A LARGE test set?
# 
# iv.  Advanced: When might you need a THIRD set (validation set)?
# 

################################################################################
#### SECTION 6: Exploratory Data Analysis (EDA) ####
################################################################################

# Danish saying: "tegne før regne" (draw before you calculate!)
# Always visualize and understand your data before modeling!

# Get all variable names except the first one (id):
var_names <- names(chd_train)[-1]

# EXERCISE E: Why do we exclude the first column with [-1]?
# Hint: Look at what the first column contains
# Answer:
# 


# EXERCISE F: Create a pairwise plot matrix
# This shows all variables plotted against each other
chd_train %>%
  ggpairs(columns = var_names)

# Note: This may take a minute to generate!

# EXERCISE G: Interpret the ggpairs output
# i. Press "Zoom" in the Plots panel for a better view.
#    If the plot is too crowded, try a smaller subset:
#    var_names_small <- c("age", "sbp", "dbp", "bmi", "scl", "chdfate")
#    Then re-run ggpairs with: columns = var_names_small

# ii. Focus on the DIAGONAL (top-left to bottom-right):
#     - What do these plots show?
#     - Why are they different plot types?
#     - What plot types do you see? (histogram, density, bar plot?)

# iii. Look at the OFF-DIAGONAL plots:
#      - What plot types appear? When?
#      - Hint: Compare numeric vs numeric, numeric vs categorical

# iv. Examine correlations with chdfate:
#     - Do any variables seem strongly related to heart disease?
#     - Which variables might be good predictors?

# Your observations:
# 


################################################################################
#### PART 2: Building Models ####
################################################################################

################################################################################
#### SECTION 7: Building Your First Linear Regression Model ####
################################################################################

# LINEAR REGRESSION BASICS:
# We're predicting a continuous outcome (BMI) from other variables.
# The formula used below is: BMI = β₀ + β₁×sbp + β₂×sex + ε
# 
# Note: chdfate is binary (TRUE/FALSE), so we'll start with BMI as our target.
# We'll get to binary classification later!

# THE TIDYMODELS WORKFLOW:
# 1. Specify the model type: linear_reg()
# 2. Fit the model to training data: fit()
# 3. Examine results: tidy()

# Example: Predict BMI from systolic blood pressure and sex
lm_fit <- 
  linear_reg() %>% 
  fit(bmi ~ sbp + sex, data = chd_train)

# View the model:
lm_fit

# Get a tidy summary of coefficients:
tidy(lm_fit)

# EXERCISE A: Interpret the results
# Look at the 'estimate' column:
# - The (Intercept) is the baseline BMI
# - Each coefficient shows how much BMI changes per unit increase in that variable
# - The p.value shows if the effect is statistically significant (< 0.05)
#
# Resources: https://argoshare.is.ed.ac.uk/healthyr_book/regression.html (Fig 7.4)
#
# Questions to consider:
# - What is the effect of sbp on BMI?
# - How does sex affect BMI?
# - Are these effects significant? 


################################################################################
#### SECTION 8: Adding More Predictors ####
################################################################################

# EXERCISE B: Build a more complex model
# Add more variables to the formula. For example, try:
# bmi ~ sbp + sex + age + scl
lm_fit2 <- 
  linear_reg() %>% 
  fit(bmi ~ sbp + sex + ?, data = chd_train)

tidy(lm_fit2)

# Compare this to the previous model. What changed?

################################################################################
#### SECTION 9: Full Model with All Predictors ####
################################################################################

# The formula "bmi ~ . - id" means:
# Predict BMI using ALL variables EXCEPT id

# EXERCISE C: Fit a model with all predictors
lm_fit_full <- 
  linear_reg() %>% 
  fit(bmi ~ . - id, data = chd_train)

tidy(lm_fit_full)

# Examine the results:
# i.   Which variables have significant coefficients (p.value < 0.05)?
# ii.  Which variables have the LARGEST effect sizes (biggest estimates)?
# iii. Which variables have the SMALLEST effect sizes?
# iv.  Advanced question: How would standardizing (scaling) variables 
#      affect the interpretation of effect sizes?
#      (We'll explore this more in the next module!)

################################################################################
#### SECTION 10: Visualizing Coefficients ####
################################################################################

# A coefficient plot (dot-and-whisker plot) helps visualize model results.

# Prepare data: remove intercept for clearer visualization
coef_data <- tidy(lm_fit_full) %>%
  filter(term != "(Intercept)")

# Create the plot:
ggplot(coef_data, aes(x = estimate, y = term)) +
  geom_errorbarh(aes(xmin = estimate - std.error, 
                     xmax = estimate + std.error)) +
  geom_point(size = 2) +
  geom_vline(xintercept = 0, lty = 2, linewidth = 0.5, col = "red") +
  labs(title = "Coefficient Plot",
       x = "Estimate (Effect Size)",
       y = "Variable") +
  theme_classic2()

# EXERCISE D: Interpret the coefficient plot
# - The dot shows the estimated effect
# - The line shows uncertainty (standard error)
# - The red dashed line is zero (no effect)
# - Coefficients that cross zero are not significant
#
# Questions:
# - Which variables have the clearest (most significant) effects?
# - Which variables are uncertain (wide error bars)?
# - Why is this visualization useful compared to just reading a table?


################################################################################
#### PART 3: Testing and Evaluation ####
################################################################################

################################################################################
#### SECTION 11: Making Predictions on Test Data ####
################################################################################

# Now for the critical step: testing on data the model has NEVER seen!

# EXERCISE E: Make predictions on test data
# Replace ? with 'chd_test'
chd_test_predictions <- augment(lm_fit_full, new_data = ?)

# EXERCISE F: Examine the predictions
# Use glimpse() or View() to see the new columns:
glimpse(chd_test_predictions)

# Notice the new columns:
# - .pred = predicted BMI values
# - .resid = residuals (actual - predicted)

################################################################################
#### SECTION 12: Evaluating Model Performance ####
################################################################################

# EXERCISE G: Plot actual vs predicted BMI
# Replace the first ? with 'bmi' (actual values)
# Replace the second ? with '.pred' (predicted values)
ggplot(chd_test_predictions, aes(x = ?, y = ?)) + 
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Actual vs Predicted BMI",
       x = "Actual BMI",
       y = "Predicted BMI") +
  coord_equal()

# The red line shows perfect predictions. Points close to it are good!

# Plot the distribution of residuals:
# Replace ? with '.resid'
ggplot(chd_test_predictions, aes(x = ?)) + 
  geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Distribution of Residuals",
       x = "Residual (Actual - Predicted)",
       y = "Count")

# EXERCISE H: Interpret the residuals
# Ideally, residuals should be:
# - Centered around zero (no systematic bias)
# - Normally distributed (bell-shaped)
# - Small magnitude (close to zero)
#
# Questions:
# - Is the distribution centered at zero?
# - Are there outliers (very large residuals)?
# - What does this tell you about model performance?

################################################################################
#### SECTION 13: Quantifying Model Performance ####
################################################################################

# EXERCISE I: Calculate Root Mean Squared Error (RMSE)
# RMSE measures the average prediction error in the same units as BMI.
# Replace first ? with 'bmi', second ? with '.pred'
chd_test_predictions %>% 
  rmse(truth = ?, estimate = ?)

# INTERPRETATION:
# RMSE tells us the typical size of prediction errors.
# For example, if RMSE = 3.5, our predictions are typically off by ±3.5 BMI points.
#
# Lower RMSE = better model performance!

# Additional metrics you can try:
# - R-squared: chd_test_predictions %>% rsq(truth = bmi, estimate = .pred)
# - Mean Absolute Error: chd_test_predictions %>% mae(truth = bmi, estimate = .pred)

################################################################################
#### KEY TAKEAWAYS ####
################################################################################

# 1. Always split data into train and test sets
# 2. Explore data before modeling (ggpairs is your friend!)
# 3. Start simple, then add complexity
# 4. Visualize coefficients to understand what matters
# 5. Evaluate on test data, not training data
# 6. Use multiple metrics (RMSE, R², visualizations) to assess performance

# NEXT STEPS:
# - Try different variable combinations
# - Consider interactions (e.g., age * sex)
# - Think about which variables make biological sense
# - Prepare for classification in the next module!
