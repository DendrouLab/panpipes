#options(repos=structure(c(CRAN="YOUR FAVORITE MIRROR")))
options(repos=structure(c(CRAN="https://cran.mirror.garr.it/CRAN/")))
#standard stuff
install.packages(c("tidyverse","remotes",
                "magrittr","Seurat",
                "optparse","yaml","stringi","Matrix",
                "reticulate","stringr","scales","clustree"))
#github stuff
remotes::install_github("PMBio/MuDataSeurat")
remotes::install_github("mojaveazure/seurat-disk")
remotes::install_github("jokergoo/ComplexHeatmap")
