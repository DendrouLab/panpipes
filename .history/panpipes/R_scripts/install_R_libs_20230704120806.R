chooseCRANmirror(ind=38)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
install.packages("systemfonts", dependencies = TRUE)
BiocManager::install(c("ComplexHeatmap"))

cran_packages <- c("foreign","clustree","remotes","Seurat")
for (pp in cran_packages){
    message("============ Installing ", pp, " ============")
    install.packages(pp)
}



remotes::install_github("grimbough/rhdf5")
remotes::install_github("ilia-kats/MuData")
remotes::install_github("pmbio/MuDataSeurat")