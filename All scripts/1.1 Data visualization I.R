#####
#
# R-programming and Visualising Data in R I
#
# Simon Grund Sorensen, Jakob Skou Pedersen, Søren Besenbacher, Farhad Zamani
# Aarhus University
#
#####

################################################################################
#### PART 1: Introduction to R Programming ####
################################################################################

################################################################################
#### SECTION 1: Introduction to R Programming ####
################################################################################

# Any line starting with a hash-tag (#) is a comment - just text!
# Comments help you document your code and leave notes for yourself.
# Try to make a new line below and write a hashtag followed by a message.
# This is my first comment in R!


# Now, let's try to run some code! To run code, move your
# cursor to the line you want to run and press CTRL+Enter (Windows)
# or CMD + Enter (Mac). Try it on the line below:

print("Hello World")

# Can you see the output in the Console below? Good! Let's try
# some more commands. Just run each line and see what happens.
print("Welcome to NCBS 2025, Bioinformatics and Genomics Workshop!")

# In R, we can store information in variables using the <- operator.
# Let's create a variable called 'my_number' and assign it the value 42:
my_number <- 42

# Now print it to see the value:
print(my_number)

# We can also do math with variables:
my_number + 8

# A data.frame is like a spreadsheet - it has rows and columns.
# Let's create a simple data.frame with two columns (x and y):
my_data <- data.frame(x = 1:5, y = 11:15)

# Look at the data by just typing its name:
my_data

# We can access individual columns using the $ symbol:
my_data$x
my_data$y

# Let's make a simple plot of our data:
plot(my_data)

# We can also get summary statistics for each column:
summary(my_data$x)
summary(my_data$y)

# CHECKPOINT: Do you understand the difference between a variable and a data.frame?
# A variable can hold a single value, while a data.frame holds multiple rows and columns.

################################################################################
#### PART 2: Introduction to ggplot2 ####
################################################################################

# First, we load the tidyverse - an essential collection of packages for data analysis.

################################################################################
#### SECTION 2: Introduction to ggplot2 ####
################################################################################
# First, we load the tidyverse - an essential collection of packages for data analysis.
# The tidyverse includes ggplot2, which we'll use for data visualization.
# Run this line by pressing CMD+Enter (Mac) or CTRL+Enter (Windows):
library(tidyverse)

# Let's load the 'mpg' dataset that is built into R. This dataset contains
# information about different car models and their fuel efficiency.
# Later you will learn to load your own data files.
mpg_data <- mpg

# EXPLORE THE DATA: Click on 'mpg_data' in the Environment panel (upper right corner).
# A data viewer will open showing the data like a spreadsheet.
# Notice the columns: manufacturer, model, displ (engine size), year, cyl (cylinders),
# hwy (highway mpg), cty (city mpg), and more.
# Close the data viewer and return to this script.

# Let's also look at a quick summary of the data:
glimpse(mpg_data)

################################################################################
#### SECTION 3: Building Your First ggplot ####
################################################################################

# ggplot builds plots in layers. Let's start with an empty plot:
ggplot(data = mpg_data)

# An empty plot isn't very useful! We need to add a "geom" (geometric object).
# Let's start simple with a histogram showing ONE variable.

# A histogram shows the distribution of 'highway miles per gallon' (hwy):
ggplot(data = mpg_data) +
  geom_histogram(aes(x = hwy))

# CHECKPOINT: A histogram needs only ONE variable (x).
# The y-axis (count) is calculated automatically.

################################################################################
#### SECTION 4: Scatter Plots - Comparing Two Variables ####
################################################################################

# Now let's compare TWO variables using a scatter plot.
# We'll look at engine size (displ) vs highway fuel efficiency (hwy).

# First, here's a complete example:
ggplot(data = mpg_data) +
  geom_point(mapping = aes(x = displ, y = hwy))

# EXERCISE 1: Now you try! Replace the '?' to plot 'hwy' on the y-axis
# (Hint: look at the example above)
ggplot(data = mpg_data) +
  geom_point(mapping = aes(x = displ, y = hwy))

################################################################################
#### SECTION 5: Adding More Aesthetics (Color, Size, Shape) ####
################################################################################

# We can map additional variables to visual properties like size and color.
# This is done INSIDE the aes() function.

# EXERCISE 2: Replace the ? to make point SIZE reflect the number of cylinders (cyl)
ggplot(data = mpg_data) + 
  geom_point(mapping = aes(x = displ, y = hwy, size = cyl))

# EXERCISE 3: Replace the ? to make point COLOR reflect the number of cylinders (cyl)
ggplot(data = mpg_data) + 
  geom_point(mapping = aes(x = displ, y = hwy, color = cyl))

# CHECKPOINT: What's the difference between these two?
# aes(color = cyl)          # Maps cylinder data to colors (INSIDE aes)
# geom_point(color = "red") # Sets ALL points to red (OUTSIDE aes)

# Let's see what the aes() function actually returns:
aes(x = displ, y = hwy)

################################################################################
#### SECTION 6: Adding Smooth Lines and Trend Lines ####
################################################################################

# We can add a smooth trend line using geom_smooth:
ggplot(data = mpg_data) +
  geom_smooth(mapping = aes(x = displ, y = hwy))

# EXERCISE 4: Make the line dotted by setting linetype parameter.
# Try values: 1 (solid), 2 (dashed), 3 (dotted), or use "dashed", "dotted"
ggplot(data = mpg_data) +
  geom_smooth(mapping = aes(x = displ, y = hwy), linetype = 3)

################################################################################
#### SECTION 7: Back to Histograms - Adjusting Bins ####
################################################################################

# EXERCISE 5: Change the geom below from geom_point to geom_histogram
ggplot(data = mpg_data) + 
  geom_histogram(mapping = aes(x = hwy))

# Histograms group data into bins. We can control the bin width or number of bins.
# EXERCISE 6: Set the number of bins to 5
# (Hint: the parameter is called 'bins')
ggplot(data = mpg_data) + 
  geom_histogram(mapping = aes(x = hwy), bins = 5)

################################################################################
#### SECTION 8: Styling Points - Fixed Colors and Shapes ####
################################################################################

# When we set aesthetics OUTSIDE aes(), they apply to ALL points equally.
# EXERCISE 7: Set color to "red" and shape to "+" (note: these go OUTSIDE aes)
ggplot(data = mpg_data, aes(x = cty, y = displ)) +
  geom_point(color = "red" , shape = "+")

################################################################################
#### SECTION 9: Layering Multiple Geoms ####
################################################################################

# One of ggplot's superpowers is layering multiple geoms!
# EXERCISE 8: Add a geom_smooth layer to the plot above
# (Hint: use + to add another layer)
ggplot(data = mpg_data, aes(x = cty, y = displ)) +
  geom_point(color = "red", shape = "+") +
  geom_smooth()
  

# FINAL STEP: Press Export > Save as PDF to save your visualization!

################################################################################
#### BONUS: Want More Practice? ####
################################################################################

# If time allows, continue with Visualization II or the exercises in R for data science chapter 2: 
#https://r4ds.hadley.nz/data-visualize#introduction 
