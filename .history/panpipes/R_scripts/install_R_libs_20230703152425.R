chooseCRANmirror(ind=38)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
install.packages("systemfonts", dependencies = TRUE)

BiocManager::install(c(
    "optparse",
    "tidyverse",
    "tidyselect",
    "Matrix",
    "foreign",
    "ggforce",
    "ggraph",
    "ComplexHeatmap",
    "xtable",
    "scales",
    "reticulate",
    "readr",
    "magrittr",
    "clustree",
    "RColorBrewer",
    "Seurat"))

remotes::install_github("grimbough/rhdf5")
remotes::install_github("ilia-kats/MuData")
remotes::install_github("pmbio/MuDataSeurat")