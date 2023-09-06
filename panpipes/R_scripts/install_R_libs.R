install.packages("systemfonts", dependencies = TRUE)
cran_packages <- c("optparse", "foreign","clustree","remotes", "table")
for (pp in cran_packages){
  message("============ Installing ", pp, " ============")
  install.packages(pp)
}

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("ComplexHeatmap", "ggraph", 
    "tidyverse", "ggforce","Seurat"))


remotes::install_github("grimbough/rhdf5")
# remotes::install_github("ilia-kats/MuData")
# remotes::install_github("pmbio/MuDataSeurat")