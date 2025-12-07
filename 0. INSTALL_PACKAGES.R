install.packages("BiocManager", repos = "https://cloud.r-project.org")

# Disable Bioconductor upgrade prompts and package updates
Sys.setenv(R_BIOC_UPGRADE = "FALSE")
BiocManager::install("M3C", update = FALSE, ask = FALSE)

packages_to_install = c(
  "tidyverse", "broom", "Rtsne", "tidymodels", "skimr", "ggdendro",
  "GGally", "dotwhisker", "vip", "ranger", "data.table", "glmnet",
  "gapminder", "riskCommunicator", "pheatmap", "ape", "cluster", "ggpubr",
  "keras", "patchwork", "keras3", "xgboost"
)

new_packages = packages_to_install[!(packages_to_install %in% installed.packages()[,"Package"])]

if(length(new_packages) > 0){
  for (pkg in new_packages) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

# Some packages need to be installed directly from github. Don't panic if this doesnt work at first,
## we can come help you (especially windows computers sometimes have isssues install devtools, and need
## and updated software called Rtools first which can be downloaded at https://cran.r-project.org/bin/windows/Rtools/)
# install.packages("devtools")
# devtools::install_github("sebastianbarfort/mapDK")
#devtools::install_github("https://github.com/eddelbuettel/bh") #Has to be done like this with slower internet
#BiocManager::install("BiocParallel")
#BiocManager::install("DESeq2")
