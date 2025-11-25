#####
# 
# Loading Data in R
#
# Simon Grund Sorensen, Jakob Skou Pedersen, Søren Besenbacher, Farhad Zamani
# Aarhus University
# 
#####

################################################################################
#### SECTION 1: Loading Required Packages ####
################################################################################

# Load the tidyverse package which includes many data import/export functions
library(tidyverse)

################################################################################
#### PART 1: Loading Data from Different File Formats ####
################################################################################

################################################################################
#### SECTION 2: Understanding Different File Formats ####
################################################################################

# In the previous exercise we used the built-in 'mpg' dataset.
# In real-world data analysis, you'll need to load data from files.
# R can read many different file formats. Let's explore the most common ones.

# NOTE: Each time we load data into 'd', we overwrite the previous version.
# This is just for demonstration - normally you'd use different variable names!

################################################################################
#### Loading R Data Files (.rds) ####
################################################################################

# .rds files are R's native format - they preserve data types perfectly
# and load quickly. This is the best format for saving R objects.
chd_data <- read_rds("Data/chd_500.rds")

# Let's look at the structure of what we just loaded:
glimpse(chd_data)

################################################################################
#### Loading CSV Files (Comma Separated Values) ####
################################################################################

# CSV files are very common - they can be opened in Excel or any text editor.
# Values are separated by commas.
chd_data <- read_csv("Data/chd_500.csv")

# Notice that read_csv() shows you information about column types!
# This helps you verify the data loaded correctly.

################################################################################
#### Loading TSV Files (Tab Separated Values) ####
################################################################################

# TSV files use tabs instead of commas to separate values.
# They're useful when your data contains commas within fields.
chd_data <- read_delim("Data/chd_500.tsv", delim = "\t")

################################################################################
#### Loading Excel Files (.xlsx) ####
################################################################################

# Excel files require the 'readxl' package
library(readxl)
chd_data <- read_excel("Data/chd_500.xlsx")

# NOTE: If your Excel workbook has multiple sheets, you can specify which one:
# chd_data <- read_excel("Data/chd_500.xlsx", sheet = "Sheet1")
# or by sheet number: sheet = 2

################################################################################
#### SECTION 3: Inspecting and Understanding Your Data ####
################################################################################

# EXERCISE 1: Load the data using any of the methods above (if you haven't already)
# The data is now stored in the variable 'chd_data'

# Let's examine the data structure:
glimpse(chd_data)

# You can also view the first few rows:
head(chd_data)

# Or open it in a spreadsheet-like viewer:
# View(chd_data)  # Uncomment and run this line to open the viewer

# QUESTION A: Look at the output of glimpse() above.
# What are the data types of each column? (chr = character/text, dbl = decimal number, int = integer)
# Write your observations as comments below:
# 


# QUESTION B: Did R detect the correct data types automatically?
# For example, should any columns be factors instead of characters?
# Should any numeric columns be integers instead of doubles?
# Write your thoughts below:
# 


################################################################################
#### SECTION 4: Converting Data Types ####
################################################################################

# Sometimes R doesn't guess the correct data type. We can convert them!

# EXERCISE 2: Convert the "chdfate" column to logical (TRUE/FALSE) instead of character
# Replace both ? marks with: chdfate
chd_data$? <- as.logical(chd_data$?)

# Verify the conversion worked:
glimpse(chd_data$chdfate)

# Common conversion functions:
# as.logical()  - converts to TRUE/FALSE
# as.numeric()  - converts to numbers
# as.character() - converts to text
# as.factor()   - converts to categorical variable

################################################################################
#### SECTION 5: Exporting Data ####
################################################################################

# After cleaning or modifying data, you often want to save it.

# EXERCISE 3: Export your formatted data as a .tsv file
# Replace the first ? with your data variable name: chd_data
# Replace the second ? with the separator character: "\t" (tab)
write.table(x = ?, file = "Data/chd_500_formatted.tsv", sep = ?, col.names = TRUE, row.names = FALSE)

# Other export options:
# write_csv(chd_data, "Data/chd_500_formatted.csv")  # Export as CSV
# write_rds(chd_data, "Data/chd_500_formatted.rds")  # Export as RDS (recommended!)

# CHECKPOINT: Do you understand the difference between the file formats?
# - .rds: R's format, preserves all data types perfectly
# - .csv: Universal format, can be opened in Excel
# - .tsv: Like CSV but uses tabs, better when data contains commas