chooseCRANmirror(ind=38)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
install.packages("systemfonts", dependencies = TRUE)

BiocManager::install(c(
    "tidyverse", 
    "tidyselect",
    "Matrix",
    "Seurat",
    "ggforce",
    "ggraph",
    "ComplexHeatmap",
    "xtable",
    "forcats",
    "purrr",
    "tidyr",
    "scales",
    "reticulate",
    "readr",
    "optparse",
    "magrittr",
    "foreign",
    "clustree",
    "ggplot2",
    "RColorBrewer"))

remotes::install_github("grimbough/rhdf5")
remotes::install_github("ilia-kats/MuData")
remotes::install_github("pmbio/MuDataSeurat")