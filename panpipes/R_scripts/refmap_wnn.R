suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)  
  library(optparse)
  library(yaml)
})

option_list <- list(
    make_option(c("--reference_path"), default = "",
                help="path to reference h5mu or seurat object"),
    make_option(c("--query_path"), default = "",
                help="path to query h5mu or seurat object"),
    make_option(c("--wnn_mapping_normalization"), default = "",
                help="wich normalization to apply/was applied to query"),
    make_option(c("--normalize"), default=FALSE,
                help="normalize query?")
    )

opt <- parse_args(OptionParser(option_list=option_list))
message("Running with options:")
print(opt)

check_and_read<- function(reference_path){
if (endsWith( reference_path, ".rds")){
    reference<- readRDS(reference_path)
    refname <- gsub(".rds","",basename(reference_path))
}

if (endsWith( reference_path, ".h5seurat")){
    library(SeuratDisk)
    reference <- LoadH5Seurat(reference_path)
    refname <- gsub(".h5seurat","",basename(reference_path))
}

if (endsWith( reference_path, ".h5mu")){
    library(MuDataSeurat)
    reference = MuDataSeurat::ReadH5MU(reference_path)
    refname <- gsub(".h5mu","",basename(reference_path))
}
if (endsWith( reference_path, ".h5ad")){
    library(SeuratDisk)
    Convert(reference_path, dest = "h5seurat", overwrite = TRUE)
    reference <- LoadH5Seurat(gsub("h5ad","h5seurat",reference_path)
    refname <- gsub(".h5ad","",basename(reference_path))
}



if (!exists("reference")){
    stop("you did not provide a valid object. 
            we support: .h5mu, Seurat saved as .rds or .h5seurat")
}
return(reference)
}

reference<- check_and_read(opt$reference_path)
query <- check_and_read(opt$query_path)


if (!is.null(opt$wnn_mapping_normalization)){
    if(opt$wnn_mapping_normalization == "SCT"){
        if (opt$normalize){
            message("will normalize query data with SCTransform")
            normchoice= "SCT"
            query <- SCTransform(query)
        }else{
            message("your query RNA is pre-normalized with SCTransform")
        }
    }
    if(opt$wnn_mapping_normalization == "LogNormalize"){
        if (opt$normalize){
            message("will normalize query data with LogNormalize")
            normchoice= "LogNormalize"
            query<- NormalizeData(query)
        }else{
            message("your query RNA is pre-normalized with LogNormalize")
        }
    }
}

#force reference reduction to go through spca
if(!"spca"  %in% names(reference@reductions)){
    message("updating reference by computing spca guided by wnn graph")
    reference <- RunSPCA(reference, assay = 'RNA', graph = 'wsnn') 
    #is this how it is called in mudata? check for compatibility
}

anchors <- FindTransferAnchors(
  reference = reference,
  query = query,
  normalization.method = normchoice,
  reference.reduction = "spca",
  dims = 1:50
)

message("Found anchors, now mapping query")
if (file.exists("pipeline.yml")){
    params<- read_yaml("pipeline.yml")
    refdata<-params["refdata"]
}else{
    stop("I don't have a valid refdata param, don't know what to transfer")
    # refdata = list(
    # celltype.l1 = "celltype.l1",
    # celltype.l2 = "celltype.l2",
    # predicted_ADT = "ADT")
}

message("will attempt to predict the following columns/continuous values matrices:")
print(refdata)

query <- MapQuery(
  anchorset = anchors,
  query = query,
  reference = reference,
  refdata = refdata,
  reference.reduction = "spca",
  reduction.model = "wnn.umap"
)

# merge reference and query and save the recomputed Umap instead of the original one
reference$is_reference <- "Reference"
query$is_reference <- "Query"
refquery <- merge(reference, query)
refquery[["spca"]] <- merge(reference[["spca"]], query[["ref.spca"]])
g<-DimPlot(refquery, group.by = 'is_reference', shuffle = TRUE)



ggsave(g, filename = paste0("figures/wnn_reference",refname,"_umap.png", height = 6, width = 8, dpi=300)

refquery <- RunUMAP(refquery, reduction = 'spca', dims = 1:50)
g<-DimPlot(refquery, group.by = 'is_reference', shuffle = TRUE)

ggsave(g, filename = "figures/wnn_recomputed",refname,"_umap.png", height = 6, width = 8, dpi=300)


umap <- FetchData(refquery, vars = c('UMAP_1','UMAP_2'))
write.csv(umap, file= "refmap/umap_",refname,"_wnn.csv", sep=",", col.names=T, row.names=T,quote=F)



saveRDS(refquery, file = "refmap/query_to_reference",refname,"_wnn.rds")
#convert to mudata

