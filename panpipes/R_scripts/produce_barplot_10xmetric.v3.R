print(.libPaths())
suppressPackageStartupMessages({library(tidyverse)
library(optparse)
library(Seurat)
library(scales)
library(ggrepel) })
options(stringsAsFactors = F)
options(bitmapType='cairo')

option_list <-  list(
  make_option("--csvpaths", default=NULL,
              help="The txt file containing the paths to metrics.cvs"),
  make_option("--outdir", default=NULL,
              help="The output directory"),
  make_option("--figdir", default=NULL,
              help="The figure directory"),
  make_option("--project", default=NULL,
              help="The project number, or a key word to describe the project"),
  make_option("--height", default=12,
              help="The heigth of the png to save, in cm."),
  make_option("--width", default=12,
              help="The width of the png to save, in cm."),
  make_option("--kneeplot", default=FALSE,
              help="Produce Knee plot? default to false cause it takes time. if set to true specify 10x Chemistry Version")

)

myFun <- function(n = 5000) {
  a <- do.call(paste0, replicate(5, sample(LETTERS, n, TRUE), FALSE))
  paste0(a, sprintf("%04d", sample(9999, n, TRUE)), sample(LETTERS, n, TRUE))
}



opt <- parse_args(OptionParser(option_list=option_list))

if(is.null(opt$outdir)) { opt$outdir <- paste0(getwd(),"/")}
if(!endsWith(opt$outdir,"/")) {opt$outdir <- paste0(opt$outdir,"/")}

if(is.null(opt$figdir)) { opt$figdir <- paste0(getwd(),"/figures/")}
if(!endsWith(opt$figdir,"/")) {opt$figdir <- paste0(opt$figdir,"/")}

if(!dir.exists(opt$figdir)) dir.create(opt$figdir)
# translate strings into bools (because cgat pipeline wil pipe strings)
if(opt$kneeplot == "True"|opt$kneeplot =="TRUE"){
  opt$kneeplot = TRUE
}else{
  opt$kneeplot = FALSE
}

print("Running with options:")
print(opt)
project <- opt$project

if(!is.null(opt$csvpaths)) {
  caf <- read.table(opt$csvpaths, header=T, sep="\t")
  if(!"gex_path" %in% colnames(caf)){
    caf$gex_path=caf$path
    caf <- caf %>% filter(filetype=="cellranger")
    if(nrow(caf)==0){
      message("no cellranger inputs")
      exit()
    }

  }
  # paths <- gsub("filtered_feature_bc_matrix","metrics_summary.csv",caf$gex_path)
  paths <- file.path(caf$gex_path, "metrics_summary.csv")
  
  # paths <- gsub("\\/$","", paths)

  #for(n in paths){ print(n); dd <- read.csv(n, header=T, stringsAsFactors = F) ; mat[n,colnames(dd)] <- gsub("%|,","", dd[1,])}
  #for(n in paths){ print(n); mat[n,] <-gsub("%|,","",read.csv(paste0(n,"metrics_summary.csv"), header=T, stringsAsFactors = F)[1,])}
  listfile <- list()
  for(n in paths){ listfile[[n]] <- read.csv(n, header=T, stringsAsFactors = F) }
  mat<-plyr::ldply(listfile, rbind)
  mat <- as.data.frame(sapply(mat,function(x) gsub("%|,","", x)))
  mat <- mat[,-1]

  if(!is.null(caf$sample_id)){
    rownames(mat)<-caf$sample_id
  }else{
    nm <- strsplit(paths, "/", fixed=T)
    f <- length(nm[[1]]) - 2
    nm <- sapply(nm,"[[",f)
    rownames(mat)<-nm
  }


}

tag.mat <- mat[,grepl("Antibody", colnames(mat))]
tag.mat <- na.omit(tag.mat)
mat <- mat[,!grepl("Antibody", colnames(mat))]
mat <- na.omit(mat)
swr = function(string, nwrap=20) {paste(strwrap(string, width=nwrap), collapse="\n")}
swr = Vectorize(swr)

if(ncol(mat)>0){
  reshape2::melt(as.matrix(mat)) %>%
    mutate(value=as.numeric(as.character(value))) %>% mutate(variable=gsub("."," ",Var2, fixed=T)) %>% mutate(variable=swr(variable)) %>%
    ggplot(aes(Var1,value, fill=Var1)) + facet_wrap(~variable, scales="free") +theme_gray()+geom_bar(stat="identity", position="dodge", color="grey20") +
    theme_bw() +
    theme(axis.text.x = element_text(angle=35,hjust=1,size=6),legend.position = "none",
          strip.text.x = element_text(size = 9, face='bold')) + xlab("") + ylab("") +
    scale_fill_manual(values = hue_pal(l = 80,c(0, 360))(nrow(mat)))-> g
  message("doing first plot")
  ggsave(g,height = opt$height ,width = opt$width, units = "cm", file=paste0(opt$figdir,"10xMetrics_",project,"_plot.pdf") )
  write.table(mat, quote=F, row.names=T,col.names=T,file=paste0(opt$outdir,"10xMetrics_",project,".txt"), sep='\t')
message("done plot")
  mat %>% rownames_to_column("sample_id") %>%
  select(Sequencing.Saturation,Number.of.Reads, Mean.Reads.per.Cell, Estimated.Number.of.Cells, sample_id) %>%
    mutate(Sequencing.Saturation=as.numeric(as.character(Sequencing.Saturation)),
           Number.of.Reads=as.numeric(as.character(Number.of.Reads)),
           Mean.Reads.per.Cell=as.numeric(as.character(Mean.Reads.per.Cell)),
           Estimated.Number.of.Cells=as.numeric(as.character(Estimated.Number.of.Cells))) %>%
    ggplot(aes(Sequencing.Saturation,Number.of.Reads,size=Mean.Reads.per.Cell, color=Estimated.Number.of.Cells,label=sample_id)) +
    geom_point() + theme_bw() + scale_color_viridis_c() + geom_text_repel(size=3, color="black") -> g2

  ggsave(g2,height = 10 ,width = 16, units = "cm", file=paste0(opt$figdir,"10xMetrics_",project,"sequencing_saturation_plot.pdf"))

}
if(ncol(tag.mat)>0){
  reshape2::melt(as.matrix(tag.mat)) %>%
    mutate(value=as.numeric(as.character(value))) %>% mutate(variable=gsub("."," ",Var2, fixed=T)) %>% mutate(variable=swr(variable)) %>%
    ggplot(aes(Var1,value, fill=Var1)) + facet_wrap(~variable, scales="free") +geom_bar(stat="identity", position="dodge") +
    theme_bw()+
    theme(axis.text.x = element_text(angle=35,hjust=1,size=9),legend.position = "none",
          strip.text.x = element_text(size = 9, face='bold')) + xlab("") + ylab("") +
    scale_fill_manual(values = hue_pal(l = 80,c(0, 360))(nrow(tag.mat)))-> g

  ggsave(g,height = opt$height ,width = opt$width, units = "cm", file=paste0(opt$figdir,"10xMetrics_Antibodytag_",project,"_plot.pdf"))
  write.table(tag.mat, quote=F, row.names=T,col.names=T,file=paste0(opt$outdir,"10xMetrics_Antibodytag_",project,".txt"), sep='\t')
  colnames(tag.mat) <- gsub("Antibody..","",colnames(tag.mat), fixed = T)
  tag.mat %>% 
    rownames_to_column("sample_id") %>%
    select(Sequencing.Saturation,Number.of.Reads, Mean.Reads.per.Cell, Antibody.Reads.Usable.per.Cell, sample_id) %>%
    mutate(Sequencing.Saturation=as.numeric(as.character(Sequencing.Saturation)),
           Number.of.Reads=as.numeric(as.character(Number.of.Reads)),
           Mean.Reads.per.Cell=as.numeric(as.character(Mean.Reads.per.Cell)),
           Antibody.Reads.Usable.per.Cell=as.numeric(as.character(Antibody.Reads.Usable.per.Cell))) %>%
    ggplot(aes(Sequencing.Saturation,Number.of.Reads,size=Antibody.Reads.Usable.per.Cell, color=Mean.Reads.per.Cell, label=sample_id)) +
    geom_point() + theme_bw() + scale_color_viridis_c() + geom_text_repel(size=3, color="black") -> g2

  ggsave(g2,height = 10 ,width = 15, units = "cm", file=paste0(opt$figdir,"10xMetrics_Antibodytag",project,"sequencing_saturation_plot.pdf"))

}

if(nrow(mat)==0){mat <- tag.mat}
#modify this to deal with AB reads and GEX reads, need to split the knee plotting in 2 if you want both GEX and ADT knees in 2 different graphs
if(opt$kneeplot){
  print("plotting kneeplot")
  gsub("metrics_summary.csv","raw_feature_bc_matrix" , paths) ->normpath
  raw<-list()
  
  for( a in 1:length(normpath)){
    print(normpath[a])
    mm<-Read10X(normpath[a])
    if(class(mm)=="list") { mm<-do.call("rbind",mm)} #modify this to deal with AB reads and GEX reads, need to split the knee plotting in 2 if you want both GEX and ADT in 2 different graphs
    colnames(mm) <- paste0(colnames(mm),"_",rownames(mat)[a])
    xx <- CreateSeuratObject(mm,assay = "RNA",  min.cells = 0, min.features = 0)
    xx@meta.data %>%
      rownames_to_column("barcode") %>% 
      arrange(desc(nCount_RNA)) %>% 
      mutate(nobarcode = rownames(mat)[a], 
      barcodes=rank(-nCount_RNA,ties.method = "first")) -> rawrank
    raw[[a]] <- rawrank
  }
  raw.data <-do.call("rbind",raw)
  gsub("raw_feature_bc_matrix", "filtered_feature_bc_matrix", normpath) -> normpath
  
  raw<-list()
  for( a in 1:length(normpath)){
    print(normpath[a])
    mm<-read.table(paste0(normpath[a],"/barcodes.tsv.gz"), header=F, sep="\t") 
    #mm<-Read10X(normpath[a])
    #if(class(mm)=="list") { mm<-do.call("rbind",mm)}
    #colnames(mm) <- paste0(colnames(mm),"_",rownames(mat)[a])
    df <- data.frame( barcode = paste(unlist(mm),rownames(mat)[a], sep="_"), nobarcode =rep(rownames(mat)[a]) )
    raw[[a]] <-df 
  }
  filtered<- do.call("rbind", raw)
  
  raw.data$passfilter <- raw.data$barcode %in% filtered$barcode

  raw.data %>%
  {ggplot(.,aes(y=nCount_RNA, x=barcodes )) + geom_line(size=2,aes(col = nobarcode  ))+
    scale_y_log10(breaks=c(1,10,100,1000,10000, 100000)) +
    scale_x_log10(breaks = c(1,10,100,1000,10000,100000)) + #+ geom_point()
    theme_bw() + geom_hline(yintercept = 1000, linetype="dashed") + 
    theme(legend.key.size = unit(0.2,"cm"), axis.text=element_text(size=13,face="bold", color="black")) +
    scale_color_manual(values = hue_pal(l = 80,c(0, 360))(nrow(mat)))} ->gb
  fl<-gb

  for (bce in unique(raw.data$nobarcode)) {
    print(bce)
    gb+ geom_line(data=dplyr::filter(raw.data, passfilter=="TRUE", nobarcode==bce), aes(y=nCount_RNA, x=barcodes), color="blue", alpha=0.4, size=3.5) +
      scale_y_log10(breaks=c(1,10,100,1000,10000, 100000)) +
      scale_x_log10(breaks = c(1,10,100,1000,10000,100000)) + #+ geom_point()
      theme_bw() + geom_hline(yintercept = 1000, linetype="dashed") + 
      theme(legend.key.size = unit(0.2,"cm"), axis.text=element_text(size=13,face="bold", color="black"))  ->gb }
  gb + fl$layers ->gb

  ggsave(gb,height = 15 ,width = 17, dpi=300, file=paste0(opt$figdir,"10xMetrics_",project,"_gex_knee_plot.png"), type="cairo-png")


}

# if ADT exists then repeat for ADT

message("finished pipeline")
