if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("TCGAbiolinks")
install.packages("ggpubr")
install.packages("dplyr")
install.packages("purrr")
