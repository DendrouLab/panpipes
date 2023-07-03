if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
    "tidyverse", 
    "tidyselect",
    "Seurat",
    "ComplexHeatmap",
    "xtable",
    "forcats",
    "purrr",
    "tidyr",
    "tibble",
    "stringr",
    "scales",
    "reticulate",
    "readr",
    "optparse",
    "magrittr",
    "foreign",
    "dplyr",
    "clustree",
    "ggraph",
    "ggplot2",
    "RColorBrewer",
    "Matrix"))

remotes::install_github("grimbough/rhdf5")
remotes::install_github("ilia-kats/MuData")
remotes::install_github("pmbio/MuDataSeurat")