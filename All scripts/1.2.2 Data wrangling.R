#####
# 
# Wrangling Data in R
#
# Simon Grund Sorensen, Jakob Skou Pedersen, Søren Besenbacher, Farhad Zamani
# Aarhus University
#
#####

################################################################################
#### SECTION 1: Introduction to Data Wrangling ####
################################################################################

# Load necessary packages
library(tidyverse)

# Data wrangling includes:
# - Filtering rows based on conditions
# - Selecting specific columns
# - Summarizing data (calculating means, counts, etc.)
# - Grouping data by categories
# - Creating new variables
# - Reshaping data (from wide to long format, or vice versa)

# We'll use the CHD (Coronary Heart Disease) dataset for these exercises.
# This dataset contains health information about patients.

# Load the data:
chd_data <- read_rds("Data/chd_full.rds")

# Let's first look at what variables we have:
glimpse(chd_data)

# Key variables:
# - age: patient age in years
# - sex: patient gender
# - bmi: Body Mass Index
# - sbp: Systolic Blood Pressure
# - dbp: Diastolic Blood Pressure
# - scl: Serum Cholesterol Level
# - chdfate: Whether patient has coronary heart disease (TRUE/FALSE)

################################################################################
#### PART 1: Basic Data Wrangling - Filter, Summarize, Group ####
################################################################################

################################################################################
#### SECTION 2: Filtering Data ####
################################################################################

# The filter() function lets us keep only rows that meet certain conditions.
# The pipe operator %>% passes data from one step to the next.

# Example: Find all patients with age over 60
chd_data %>%
  filter(age > 60)

# EXERCISE A: Replace the ? to find how many individuals have a BMI over 50
# (Hint: look at the variable names above. BMI is stored as 'bmi')
chd_data %>%
  filter(bmi > 50)


# You can also combine multiple conditions using & (AND) or | (OR)

# EXERCISE B: Find individuals over 50 years AND cholesterol (scl) over 300
# (Hint: use filter() with two conditions connected by &)
date_formatted =
  chd_data %>%
  filter(age > 50 & scl > 300)

################################################################################
#### SECTION 3: Summarizing Data ####
################################################################################

# The summarise() function creates summary statistics from your data.

# Example: Calculate the mean BMI
chd_data %>%
  summarise(mean_bmi = mean(bmi, na.rm = TRUE))

# EXERCISE C: Calculate the mean age in the dataset
# Replace the ? with the appropriate function and variable
chd_data %>%
  summarise(mean_age = mean(age))

################################################################################
#### SECTION 4: Grouping Data ####
################################################################################

# The group_by() function lets us calculate summaries for different groups.

# Example: Calculate mean BMI separately for each chdfate status
chd_data %>%
  group_by(chdfate) %>%
  summarise(mean_bmi = mean(bmi, na.rm = TRUE))

# EXERCISE D: Calculate the mean age for males and females separately
# Replace the first ? with the grouping variable (sex)
# Replace the second ? with a summarise statement
chd_data %>%
  group_by(sex) %>%
  summarise(mean_age = mean(age))

################################################################################
#### SECTION 5: Handling Missing Values (NA) ####
################################################################################

# Real-world data often has missing values represented as NA (Not Available).

# Let's create some missing data to demonstrate:
# WARNING: This modifies row 100 - we'll reload the data after this section!
chd_data[100, ] <- NA

# This will produce an error or NA because of the missing value:
chd_data %>%
  summarise(mean_age = mean(age))

# EXERCISE E: Look at the output above. What happened?

# EXERCISE F: Fix the code to calculate mean age while ignoring NA values
# (Hint: add the parameter na.rm = TRUE inside the mean() function)
chd_data %>%
  summarise(mean_age = mean(age, ?))

# Let's reload the clean data for the rest of the exercises:
chd_data <- read_rds("Data/chd_full.rds")

################################################################################
#### PART 2: Advanced Wrangling - Mutate, Complex Operations ####
################################################################################

################################################################################
#### SECTION 6: Creating New Variables with mutate() ####
################################################################################

# The mutate() function creates new columns based on existing ones.

# Example: Calculate age in months
chd_data %>%
  mutate(age_months = age * 12)

# EXERCISE G: Create a new variable showing the ratio of systolic to diastolic blood pressure
# Replace the ? with the calculation (sbp / dbp)
chd_data %>%
  mutate(sbp_dbp_ratio = ?)

################################################################################
#### SECTION 7: Finding Maximum and Minimum Values ####
################################################################################

# slice_max() and slice_min() help us find rows with extreme values.

# EXERCISE H: Find the patient with the largest sbp/dbp ratio
# First create the ratio, then use slice_max()
chd_data %>%
  mutate(sbp_dbp_ratio = sbp / dbp) %>%
  slice_max(?)

# EXERCISE I: Find the smallest cholesterol (scl) level among individuals with BMI over 40
# (Hint: first filter for BMI > 40, then use slice_min on scl)
chd_data %>%
  ? %>%
  ?

################################################################################
#### SECTION 8: Combining Multiple Operations ####
################################################################################

# We can chain together filter, group_by, and summarise for complex analyses.

# EXERCISE J: Among patients with heart disease (chdfate == TRUE), 
# what is the MEDIAN scl level for each gender?
# (Hint: use filter, then group_by sex, then summarise with median())
chd_data %>%
  filter(?) %>%
  group_by(?) %>%
  summarise(median_scl = ?)

################################################################################
#### SECTION 9: Creating Categorical Variables ####
################################################################################

# ifelse() lets us create binary (yes/no) variables based on conditions.

# Example: Create a variable for high BMI (over 30)
chd_data %>%
  mutate(high_bmi = ifelse(bmi > 30, TRUE, FALSE))

# EXERCISE K: Create a column called age_over_30 indicating if age > 30
# Replace the two ? marks with TRUE and FALSE
chd_data <- chd_data %>%
  mutate(age_over_30 = ifelse(age > 30, ?, ?))

################################################################################
#### SECTION 10: Multi-way Grouping ####
################################################################################

# You can group by multiple variables at once to see combinations.

# EXERCISE L: Calculate mean scl for each combination of chdfate AND age_over_30
# (Hint: put multiple variables in group_by, separated by commas)
chd_data %>%
  group_by(?, ?) %>%
  summarise(mean_scl = mean(scl, na.rm = TRUE))

# EXERCISE M: Count how many patients are in each combination of chdfate and age_over_30
# (Hint: use n() inside summarise to count rows)
chd_data %>%
  group_by(chdfate, age_over_30) %>%
  summarise(count = ?)

################################################################################
#### SECTION 11: Combining Data Wrangling with Visualization ####
################################################################################

# One of the most powerful aspects of the tidyverse is combining wrangling and plotting!

# EXERCISE N: Calculate mean scl for each Age and chdfate combination
summary_data <- chd_data %>%
  group_by(?, ?) %>%
  summarise(mean_scl = mean(scl, na.rm = TRUE))

# EXERCISE O: Make a scatter plot with mean_scl on x-axis, age on y-axis, 
# colored by chdfate
# (Hint: use the summary_data from above)
summary_data %>%
  ggplot(aes(x = ?, y = ?, color = ?)) +
  geom_point()

# EXERCISE P: Extend the plot to show point SIZE based on number of observations
# First, modify the summary to include a count
summary_data <- chd_data %>%
  group_by(age, chdfate) %>%
  summarise(
    mean_scl = mean(scl, na.rm = TRUE),
    n_patients = ?  # Replace ? with n() to count observations
  )

# Then add size to the plot:
summary_data %>%
  ggplot(aes(x = mean_scl, y = age, color = chdfate, size = ?)) +
  geom_point()

# EXERCISE Q: Remove points based on fewer than 10 observations
# (Hint: add a filter() step before plotting)
summary_data %>%
  filter(? >= 10) %>%
  ggplot(aes(x = mean_scl, y = age, color = chdfate, size = n_patients)) +
  geom_point()

# CONGRATULATIONS! You've learned the core data wrangling operations:
# filter, select, mutate, group_by, summarise, and combining with ggplot!