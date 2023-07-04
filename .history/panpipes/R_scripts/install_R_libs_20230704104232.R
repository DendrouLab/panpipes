chooseCRANmirror(ind=38)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
install.packages("systemfonts", dependencies = TRUE)

cran_packages <- c("ggforce","ggraph",
                    "foreign","xtable","clustree","remotes")
for (pp in cran_packages){
    message("============ Installing ", pp, " ============")
    install.packages(pp)
}

BiocManager::install(c("ComplexHeatmap"))
remotes::install_github("grimbough/rhdf5")

install.packages("Seurat")
remotes::install_github("ilia-kats/MuData")
remotes::install_github("pmbio/MuDataSeurat")