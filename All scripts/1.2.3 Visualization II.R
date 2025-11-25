#####
# 
# Visualising Data in R II - Advanced ggplot2
#
# Simon Grund Sorensen, Jakob Skou Pedersen, Søren Besenbacher, Farhad Zamani
# Aarhus University
# 
#####

################################################################################
#### SECTION 1: Introduction to the Framingham Heart Study Dataset ####
################################################################################

# Load required packages
library(tidyverse)

# We'll use data from the Framingham Heart Study - a famous long-term 
# cardiovascular study. The full dataset has 4,699 observations.
# Read more: https://www.causeweb.org/tshs/framingham-didactic/

# Load a 500-patient subsample:
chd_data <- read_rds("Data/chd_500.rds")

# Examine the structure:
glimpse(chd_data)

# KEY VARIABLES (see "Framingham Didactic Data Dictionary_2.pdf" for details):
# - sbp: Systolic blood pressure (mm Hg)
# - dbp: Diastolic blood pressure (mm Hg) 
# - scl: Serum cholesterol (mg/100 mL)
# - bmi: Body Mass Index
# - age: Patient age (years)
# - sex: Patient sex (male/female)
# - month: Month of enrollment
# - chdfate: Whether patient developed coronary heart disease (TRUE/FALSE)

################################################################################
#### PART 1: Exploring Plot Types ####
################################################################################

################################################################################
#### SECTION 2: Scatter Plots and Regression Lines ####
################################################################################

# One powerful feature of ggplot is that you can save the base plot and reuse it!

# Example: Create a base plot showing systolic vs diastolic blood pressure
base_plot <- ggplot(chd_data, aes(sbp, dbp))
base_plot + geom_point()

# EXERCISE A: Create a plot comparing BMI (x-axis) against DBP (y-axis)
# Replace the ?'s with variable names
my_plot <- ggplot(chd_data, aes(?, ?)) 
my_plot + geom_point()

# EXERCISE B: Add a linear regression line to show the trend
# Replace the ? with the parameter name 'method'
my_plot + 
  geom_point() + 
  geom_smooth(? = "lm")

# CHECKPOINT: What does the regression line tell you about the relationship between BMI and blood pressure?

################################################################################
#### SECTION 3: Bar Plots for Categorical Data ####
################################################################################

# Bar plots are great for showing counts of categorical variables.

# Example: Count patients by sex
ggplot(chd_data) + 
  geom_bar(aes(sex))

# EXERCISE C: Make a bar plot showing patients enrolled per month
# Replace ? with geom_bar
ggplot(chd_data) + 
  ?(aes(month))

# We can also show TWO categorical variables at once using stacked bars.

# EXERCISE D: Make a stacked bar plot of enrollment by month, filled by sex
# Replace the first ? with 'sex' and the second ? with "stack"
ggplot(chd_data) + 
  geom_bar(aes(month, fill = ?), position = ?)

################################################################################
#### SECTION 4: Box Plots, Violin Plots, and Dot Plots ####
################################################################################

# These plots show the distribution of a continuous variable across categories.

# Example: Distribution of age by sex
ggplot(chd_data, aes(sex, age)) + 
  geom_boxplot()

# EXERCISE E: Make a box plot of sex vs BMI
# Replace ? with geom_boxplot
ggplot(chd_data, aes(sex, bmi)) + 
  geom_?()

# EXERCISE F: Change it to a violin plot (shows the full distribution shape)
# (Hint: replace geom_boxplot with geom_violin)
ggplot(chd_data, aes(sex, bmi)) + 
  ?

# EXERCISE G: Make a dot plot where each dot represents a patient
# Replace the first ? with 'geom_dotplot' and second ? with a number like 1
ggplot(chd_data, aes(sex, bmi)) + 
  geom_?(binaxis = "y", stackdir = "center", binwidth = ?)

# CHECKPOINT: Which plot type (box, violin, or dot) do you find most informative? Why?

################################################################################
#### PART 2: Customizing Plots - Colors, Scales, Labels ####
################################################################################

################################################################################
#### SECTION 5: Customizing Colors ####
################################################################################

# ggplot has many built-in color palettes to make your plots more attractive.


# First, let's create a base plot we'll reuse:
bp_plot <- ggplot(chd_data, aes(sbp, dbp))

# CONTINUOUS COLOR SCALES (for numeric variables like BMI):

# Example: Default color scale
bp_plot + geom_point(aes(color = bmi))

# EXERCISE A: Try a better color scale for continuous data
# Add one of these after the geom_point:
# - scale_color_viridis_c() (colorblind-friendly!)
# - scale_colour_gradientn(colours = rainbow(4))
bp_plot + 
  geom_point(aes(color = bmi)) + 
  ?

# DISCRETE COLOR SCALES (for categorical variables):

# EXERCISE B: Try the Brewer color palette for categorical data
# Add: scale_colour_brewer(palette = "Paired")
bp_plot + 
  geom_point(aes(color = sex)) +
  ?

# Try it with month too:
bp_plot + 
  geom_point(aes(color = month)) +
  scale_colour_brewer(palette = "Paired")

################################################################################
#### SECTION 6: Adjusting Axes and Labels ####
################################################################################

# EXERCISE C: Reverse the x-axis to show high values on the left
# Add scale_x_reverse() to the plot
ggplot(chd_data, aes(sbp, dbp)) +
  geom_point(aes(color = month)) + 
  ?

# EXERCISE D: Add informative labels and a title
# Use: xlab("Label"), ylab("Label"), ggtitle("Title")
ggplot(chd_data, aes(sbp, dbp)) + 
  geom_point(aes(color = month)) + 
  xlab("?") +
  ylab("?") +
    ggtitle("?")

################################################################################
#### PART 3: Faceting and Advanced Layouts ####
################################################################################

################################################################################
#### SECTION 7: Faceting - Multiple Subplots ####
################################################################################

# Faceting creates multiple small plots, one for each category.
# This is incredibly powerful for comparing patterns across groups!

# Example: Separate plot for each sex
ggplot(chd_data, aes(sbp, dbp)) + 
  geom_point() + 
  facet_wrap(~sex)

# EXERCISE A: Create separate plots for each month
# Replace ? with 'month'
ggplot(chd_data, aes(sbp, dbp, color = month)) + 
  geom_point() + 
  facet_wrap(~?)

# EXERCISE B: Arrange the facets in 4 rows
# Add the parameter nrow = 4 inside facet_wrap
ggplot(chd_data, aes(sbp, dbp, color = month)) + 
  geom_point() + 
  facet_wrap(~month, ?)

# EXERCISE C: Instead of rows, specify 3 columns using ncol
# (Hint: use ncol = 3)
ggplot(chd_data, aes(sbp, dbp, color = month)) + 
  geom_point() + 
  facet_wrap(~month, ?)

# FACET_GRID: For two categorical variables, use facet_grid

# EXERCISE D: Create a grid with sex as rows and chdfate as columns
# Format: facet_grid(rows ~ columns)
my_facet_plot <- ggplot(chd_data, aes(sbp, dbp, color = month)) + 
  geom_point() + 
  facet_grid(? ~ ?)

# Print the plot:
my_facet_plot

# CHECKPOINT: How does facet_wrap differ from facet_grid?
# - facet_wrap: wraps panels into multiple rows/columns (for one variable)
# - facet_grid: creates a matrix layout (for two variables)

################################################################################
#### PART 4: Advanced Topics and Exporting ####
################################################################################

################################################################################
#### SECTION 8: Advanced Exercise - Proportions ####
################################################################################

# EXERCISE E (BONUS): Create a plot showing the PROPORTION of heart disease 
# by sex, not just counts.
# Hint: Use geom_bar with position = "fill"
# See R4DS section 2.5.2 for guidance

ggplot(chd_data, aes(x = sex, fill = chdfate)) + 
  geom_bar(position = "fill") +
  ylab("Proportion") +
  ggtitle("Proportion of Heart Disease by Sex")


################################################################################
#### SECTION 9: Saving Your Plots ####
################################################################################

# Once you've created a beautiful plot, you'll want to save it!
# Use ggsave() to export plots in various formats.

# Example: Save the last plot you created
# ggsave("Output/my_plot.pdf", width = 20, height = 15, units = "cm")

# EXERCISE A: Create a plot and save it as PDF
my_facet_plot  # This displays the plot
ggsave("Output/chd_facet_plot.pdf", width = 20, height = 20, units = "cm")

# EXERCISE B: Save the same plot as PNG (better for presentations/websites)
# Just change the file extension to .png
ggsave("Output/chd_facet_plot.png", width = 20, height = 20, units = "cm")

# EXERCISE C: Download the plots to your laptop
# (If working on a server: look in the Output folder and download the files)

# TIP: You can control the resolution of PNG files with the 'dpi' parameter:
# ggsave("plot.png", dpi = 300)  # High quality for publications

################################################################################
#### SECTION 10: Bonus - Plotting Maps (Optional) ####
################################################################################

# R can also create geographic visualizations!

# EXERCISE D: Plot a map of Denmark
# First install the package (if not already installed):
# install.packages("devtools")
# devtools::install_github("sebastianbarfort/mapDK")

library(mapDK)
mapDK()

# EXERCISE E: Think about adding data to this map
# What variables would you need in a dataframe to add data points?
# Consider:
# - Geographic coordinates (latitude, longitude)
# - Region/municipality names
# - The data values you want to visualize

# EXERCISE F (ADVANCED BONUS): Plot a map of India
# Resources:
# - https://www.rpubs.com/spoonerf/countrymapggplot2
# - Search: "ggplot2 map India"
# Note: India is much larger than Denmark, so this may require more computing power!

################################################################################
#### CONGRATULATIONS! ####
################################################################################

# You've learned:
# - Different plot types (scatter, bar, box, violin, dot)
# - Customizing colors with palettes
# - Adding labels and titles
# - Faceting for small multiples
# - Saving plots in different formats
# - Basic geographic visualization

# Keep practicing with your own data!