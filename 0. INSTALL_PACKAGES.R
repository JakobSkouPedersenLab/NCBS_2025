install.packages("BiocManager", repos = "https://cloud.r-project.org")
BiocManager::install("M3C")

packages_to_install = c(
  "tidyverse", "broom", "Rtsne", "tidymodels", "skimr", "ggdendro",
  "GGally", "dotwhisker", "vip", "ranger", "data.table", "glmnet",
  "gapminder", "riskCommunicator", "pheatmap", "ape", "cluster", "ggpubr",
  "keras", "patchwork"
)

new_packages = packages_to_install[!(packages_to_install %in% installed.packages()[,"Package"])]

if(length(new_packages) > 0){
  for (pkg in new_packages) {
    install.packages(pkg, ask = FALSE, Ncpus = 4L)
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
