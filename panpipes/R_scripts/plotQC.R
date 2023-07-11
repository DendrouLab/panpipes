# plotting utilities for qc 
suppressPackageStartupMessages({
  library(tidyverse)
  library(optparse)
})

do_scatter_plot <- function(df, x, y, facet=NULL, hue=NULL){
  g <- df %>% 
    ggplot(aes_string(x=x,y=y, color=hue)) + 
    geom_point(size=0.5) 
  if (!is.null(facet)){
    g <- g + facet_wrap(as.formula(paste0("~", sc)), ncol=6) 
  }
  
  if(is.numeric(df[hue])){
    g <- g + scale_color_viridis_c()
  }
  g <- g+ theme_bw()+
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
          strip.text.x = element_text(size=6))
  return(g)
}

do_violin_plot <- function(df, qc, group){
  g <- df %>% drop_na() %>%
    ggplot(aes_string(x="sample_id",y=qc)) +
    geom_violin(aes_string(fill=group)) + 
    {if(qc=="doublet_scores") geom_hline(yintercept=0.25, color="grey50", linetype="dashed") }+
    {if("pct_counts" %in% qc) geom_hline(yintercept=c(5, 10, 20 ,70), color="grey50", linetype="dashed")}+
    {if("pct_counts" %in% qc) coord_cartesian(ylim=c(0,100))}+
    theme_bw()+
    theme(axis.text.x=element_text(size=8, angle=90),
          axis.text.y=element_text(size=8))
  if(length(unique(df[[group]])) > 10){
    message(paste0(sc, "has too many categories, removing legend"))
    g <- g + theme(legend.position="none")
  }
  
  return(g)
}

do_bar_plot <- function(df, qc, group){
  g <- df %>%
    ggplot(aes_string(x=group, fill=qc)) +
    geom_bar(aes_string()) + 
    theme_bw()+
    theme(axis.text.x=element_text(size=8, angle=90),
          axis.text.y=element_text(size=8))
#   if(length(unique(df[[group]])) > 10){
#     message(paste0(sc, "has too many categories, removing legend"))
#     g <- g + theme(legend.position="none")
#   }
  
  return(g)
}

options(stringsAsFactors = F)
options(bitmaptype="cairo")

option_list <- list(
  make_option(c("--cell_metadata"), default=NULL,
              help="the path to the anndata object"),
  make_option(c("--groupingvar"), default="sample_id,tissue,patient,channel",
              help="names of grouping variables"),
  make_option(c("--rna_qc_metrics"), default=NULL, 
              help="the qc_metrics to plot"),
  make_option(c("--prot_qc_metrics"), default=NULL, 
              help="the qc_metrics to plot"),
  make_option(c("--rep_qc_metrics"), default=NULL, 
              help="the qc_metrics to plot"),
  make_option(c("--atac_qc_metrics"), default=NULL, 
              help="the qc_metrics to plot"),
  make_option(c("--outdir"), default="./figures/",
              help="the name of the output folder"),
  make_option(c("--prefilter"), default=TRUE,
              help="am i parsing data before or after filtering?"),
  make_option(c("--sampleprefix"), default="",
              help="the prefix to prepend to save summary filtering plots"),
  make_option(c("--scanpy_or_muon"), default="scanpy", 
              help="was the input file written out from the obs of scanpy or muon")
)


message("Plot QC data")
opt <- parse_args(OptionParser(option_list=option_list))

if(is.null(opt$outdir)) { opt$outdir <- paste0(getwd(),"/")}
if(!grepl("\\/$", opt$outdir)){opt$outdir <- paste(opt$outdir, "/", sep = "")}
if(!file.exists(opt$outdir)){dir.create(opt$outdir)}

print("Running with options:")
print(opt)

run<- opt$outdir
opt[which(opt=="NULL")] <- NULL
opt[which(opt=="None")] <- NULL
opt$prefilter <- as.logical(opt$prefilter)

# load_data --------------------------------------------------------------------

data_plot = read.delim(opt$cell_metadata)

# define source facet for all plots
if (!is.null(opt$groupingvar)){
  source_facet <- strsplit(opt$groupingvar,",")[[1]]
  check_excl <- source_facet[!source_facet %in% colnames(data_plot)] 
  keep_source <- NULL
  for (cc in check_excl){
    for (mod in c("rna","prot","atac","rep")){
      id <- paste(mod,cc,sep=".")
      if(id %in% colnames(data_plot)){
        keep_source <- c(keep_source,id)
      }
    }
  }

  source_facet <- source_facet[source_facet %in% colnames(data_plot)]
  source_facet <- unique(c(source_facet, keep_source) )
  if(length(source_facet)>0){
    print("Facet plotting with")
    # add sample_id as a minimum requirement if it's not there already
    source_facet = unique(c("sample_id", source_facet))
    print(source_facet)
  }
}else{
  stop("i don't have the minimum variable _sampleid_ to facet on, will stop here")
}
# RNA plots --------------------------------------------------------------------


if(opt$scanpy_or_muon=="scanpy"){
  rna_data_plot <- data_plot
}else{
  rna_data_plot <- data_plot[,grep("^rna\\.",colnames(data_plot))]
  colnames(rna_data_plot) <- gsub("^rna\\.", "", colnames(rna_data_plot))
}

outpath = file.path(run, "rna")
if (!dir.exists(outpath)) dir.create(outpath)


# this is not forced anymore
# I'm removing the illusion of choice because QC metrics only applies to the violin plots anyway.
# qcmetrics=c("total_counts",
#             "log1p_total_counts",
#             "n_genes_by_counts",
#             "log1p_n_genes_by_counts",
#             "doublet_scores", 
#             "pct_counts_mt", 
#             "pct_counts_rp", 
#             "pct_counts_ig",
#             "pct_counts_hb")

# check these qc metrics are in the file
if (!is.null(opt$rna_qc_metrics)) {
  qcmetrics <- strsplit(opt$rna_qc_metrics,",")[[1]]
}
qcmetrics <- qcmetrics[qcmetrics %in% colnames(rna_data_plot)]
uniq_sample_id <- nrow(unique(rna_data_plot["sample_id"]))
rna_source_facet <- gsub("^rna\\.", "",grep("^rna.", source_facet, value = TRUE))
rna_source_facet <- unique(c(rna_source_facet, source_facet[!grepl("^rna.", source_facet)]))
rna_source_facet <- rna_source_facet[rna_source_facet %in% colnames(rna_data_plot)]
for (qc in qcmetrics){
  print(qc)
  for (sc in rna_source_facet){ #add gsub temp here
    print(sc)
    g <- do_violin_plot(rna_data_plot, qc, sc)
    if (uniq_sample_id  > 50){width=12}else{width=6}
    ggsave(g, filename=file.path(outpath, paste0("violin_", sc, "_rna-", qc,".png")), type="cairo", width= width, height=6)
    
  }
}

for (sc in rna_source_facet){
  uniq_source <- nrow(unique(rna_data_plot[sc]))
  if(uniq_source >6){
    ncols=6
    nrows=ceiling(uniq_source/6)
  }else{
    ncols=uniq_source
    nrows=1
  }
  # plot once per source facet
  if (all(c("total_counts","n_genes_by_counts")%in% colnames(rna_data_plot))){
    g <- do_scatter_plot(rna_data_plot,x="total_counts",y="n_genes_by_counts", facet=sc)
    ggsave(g, filename=file.path(outpath, paste0("scatter_", sc, "_rna-nUMI_v_rna-genes.png")), type="cairo",
           width= 2*ncols, height=2*nrows, dpi=120)
  }
  
  if (all(c("log1p_total_counts","log1p_n_genes_by_counts")%in% colnames(rna_data_plot))){
    
    g <- do_scatter_plot(rna_data_plot,x="log1p_total_counts",y="log1p_n_genes_by_counts", facet=sc)
    ggsave(g, filename=file.path(outpath, paste0("scatter_", sc, "_rna-log1p_nUMI_v_rna-log1p_genes.png")), type="cairo",
           width= 2*ncols, height=2*nrows, dpi=120)
  }
  if (all(c("total_counts","pct_counts_mt")%in% colnames(rna_data_plot))){
    g <- do_scatter_plot(rna_data_plot,x="total_counts",y="pct_counts_mt", facet=sc)
    ggsave(g, filename=file.path(outpath, paste0("scatter_", sc, "_rna-nUMI_v_rna-pct_mt.png")), type="cairo",
           width= 2*ncols, height=2*nrows, dpi=120)
  }
  if (all(c("n_genes_by_counts","doublet_scores","total_counts")%in% colnames(rna_data_plot))){
    g <- do_scatter_plot(rna_data_plot,x="n_genes_by_counts",y="doublet_scores",
                         hue="total_counts", facet=sc)
    ggsave(g, filename=file.path(outpath, paste0("scatter_", sc, "_rna-genes_rna-doublet_scores_rna-numi.png")),    
           type="cairo", width= 2*ncols, height=2*nrows, dpi=120)
  }
  
}



# Protein plots ----------------------------------------------------------------


if(!is.null(opt$prot_qc_metrics)){
  qcmetrics <- strsplit(opt$prot_qc_metrics,",")[[1]]
  prot_data_plot <- data_plot[,grep("^prot\\.",colnames(data_plot))]
  colnames(prot_data_plot) <- gsub("^prot\\.", "", colnames(prot_data_plot))
  prot_source_facet <- gsub("^prot\\.", "",grep("^prot.", source_facet, value = TRUE))
  prot_source_facet <- unique(c(prot_source_facet, source_facet[!grepl("^prot.", source_facet)]))
  prot_source_facet <- prot_source_facet[prot_source_facet %in% colnames(prot_data_plot)]

  outpath = file.path(run, "prot")
  if (!dir.exists(outpath)) { 
    dir.create(outpath)
    }
  uniq_sample_id <- nrow(unique(prot_data_plot["sample_id"]))

  # check these qc metrics are in the file
  qcmetrics <- qcmetrics[qcmetrics %in% colnames(prot_data_plot)]
  for (qc in qcmetrics){
    print(qc)
    for (sc in prot_source_facet){
      g <- do_violin_plot(prot_data_plot, qc, sc)
      if (uniq_sample_id  > 50){width=12}else{width=6}
      ggsave(g, filename=file.path(outpath, paste0("violin_", sc, "_prot-", qc,".png")), type="cairo", width= width, height=6)
    }
  }
  
  # do the following plots:
  # - prot:total_counts vs prot_nadt_by_counts
  # - prot:total_counts vs prot:isotype_counts
  # - prot:log1p_total_counts vs prot:log1p_isotype_counts
  # - prot:total_counts vs prot:pct_isotype_counts
  
  
  for (sc in prot_source_facet){
    uniq_source <- nrow(unique(prot_data_plot[sc]))
    if(uniq_source >6){
      ncols=6
      nrows=ceiling(uniq_source/6)
    }else{
      ncols=uniq_source
      nrows=1
    }
    # plot once per source facet
    message("protein scatter plots")  
    if (all(c("total_counts","n_adt_by_counts")%in% colnames(prot_data_plot))){
      g <- do_scatter_plot(prot_data_plot,x="total_counts",y="n_adt_by_counts", facet=sc)
      ggsave(g, filename=file.path(outpath, paste0("scatter_", sc, "_prot-nUMI_v_prot-adt.png")), type="cairo",
             width= 2*ncols, height=2*nrows, dpi=120)
    }
    if (all(c("log1p_total_counts","log1p_n_adt_by_counts")%in% colnames(prot_data_plot))){
      g <- do_scatter_plot(prot_data_plot,x="log1p_total_counts",y="log1p_n_adt_by_counts", facet=sc)
      ggsave(g, filename=file.path(outpath, paste0("scatter_", sc, "_prot-log1p_nUMI_v_prot-log1p_adt.png")), type="cairo",
             width= 2*ncols, height=2*nrows, dpi=120)
    }
    if (all(c("total_counts","pct_counts_isotype")%in% colnames(prot_data_plot))){
      g <- do_scatter_plot(prot_data_plot,x="total_counts",y="pct_counts_isotype", facet=sc)
      ggsave(g, filename=file.path(outpath, paste0("scatter_", sc, "_prot-nUMI_v_prot-pct_isotype.png")), type="cairo",
             width= 2*ncols, height=2*nrows, dpi=120)
    }
    if (all(c("total_counts","total_counts_isotype")%in% colnames(prot_data_plot))){
      g <- do_scatter_plot(prot_data_plot,x="total_counts",y="total_counts_isotype", facet=sc)
      ggsave(g, filename=file.path(outpath, paste0("scatter_", sc, "_prot-nUMI_v_prot-total_counts_isotype.png")), type="cairo",
             width= 2*ncols, height=2*nrows, dpi=120)
    }
    if (all(c("log1p_total_counts","log1p_total_counts_isotype", "isotype_exclude")%in% colnames(prot_data_plot))){
      g <- do_scatter_plot(prot_data_plot,x="log1p_total_counts",y="log1p_total_counts_isotype", hue="isotype_exclude", facet=sc)
      ggsave(g, filename=file.path(outpath, paste0("scatter_", sc, "_prot-log1p_nUMI_v_prot-log1p_total_counts_isotype.png")), type="cairo",
             width= 2*ncols, height=2*nrows, dpi=120)
    }
  }

}
# Atac plots ----------------------------------------------------------------


if(!is.null(opt$atac_qc_metrics)){
  message("Atac plots")
  atac_data_plot <- data_plot[,grep("^atac\\.",colnames(data_plot))]
  colnames(atac_data_plot) <- gsub("^atac\\.", "", colnames(atac_data_plot))
  atac_source_facet <- gsub("^atac\\.", "",grep("^atac.", source_facet, value = TRUE))
  atac_source_facet <- unique(c(atac_source_facet, source_facet[!grepl("^atac.", source_facet)]))
  atac_source_facet <- atac_source_facet[atac_source_facet %in% colnames(atac_data_plot)]


  outpath = file.path(run, "atac")
  if (!dir.exists(outpath)) dir.create(outpath)
  
  uniq_sample_id <- nrow(unique(atac_data_plot["sample_id"]))

  # check these qc metrics are in the file
  qcmetrics <- qcmetrics[qcmetrics %in% colnames(atac_data_plot)]
  for (qc in qcmetrics){
    print(qc)
    for (sc in atac_source_facet){
      g <- do_violin_plot(atac_data_plot, qc, sc)
      if (uniq_sample_id  > 50){width=12}else{width=6}
      ggsave(g, filename=file.path(outpath, paste0("violin_", sc, "_atac-", qc,".png")), type="cairo", width= width, height=6)
    }
  }
}
  
# Rep plots ----------------------------------------------------------------


if (!is.null(opt$rep_qc_metrics)) {
  message("Repertoire plots")
  qcmetrics <- strsplit(opt$rep_qc_metrics,",")[[1]]
  qcmetrics <- gsub("rep:", "", qcmetrics)
  rep_data_plot <- data_plot[,grep("^rep\\.",colnames(data_plot))]
  colnames(rep_data_plot) <- gsub("^rep\\.", "", colnames(rep_data_plot))
  rep_data_plot = rep_data_plot %>% filter(sample_id!="")
  
  rep_source_facet <- gsub("^rep\\.", "",grep("^rep.", source_facet, value = TRUE))
  rep_source_facet <- unique(c(rep_source_facet, source_facet[!grepl("^rep.", source_facet)]))
  rep_source_facet <- rep_source_facet[rep_source_facet %in% colnames(rep_data_plot)]

  
  outpath = file.path(run, "rep")
  if (!dir.exists(outpath)) dir.create(outpath)
  
  
  uniq_sample_id <- nrow(unique(rep_data_plot["sample_id"]))

  # check these qc metrics are in the file
  qcmetrics <- qcmetrics[qcmetrics %in% colnames(rep_data_plot)]
  for (qc in qcmetrics){
    for (sc in rep_source_facet){
      g <- do_bar_plot(rep_data_plot, qc, sc)
      if (uniq_sample_id  > 50){width=12}else{width=6}
        ggsave(g, filename=file.path(outpath, paste0("bar_", sc, "_rep-", qc,".png")), type="cairo", width= width, height=6)

      if (!(qc %in% c('has_ir', "receptor_type"))){
        g <- do_bar_plot(rep_data_plot, qc, sc) + facet_grid(~receptor_type)
        ggsave(g, filename=file.path(outpath, paste0("bar_facet_", sc, "_rep-", qc,".png")), type="cairo", width= width*4, height=6)

      }
    }
  }
}


if(!is.null(opt$prot_qc_metrics)){

# rna vs prot plots ----------------------------------------------------------------
  

  outpath = file.path(run, "rna_v_prot")
  if (!dir.exists(outpath)) dir.create(outpath)    
  

  for (sc in source_facet){
    
    message(sc)
    uniq_source <- nrow(unique(data_plot[sc]))
    if(uniq_source >6){
      ncols=6
      nrows=ceiling(uniq_source/6)
    }else{
      ncols=uniq_source
      nrows=1
    }
    message("rna v protein scatter plots")  
    if (all(c("rna.total_counts","prot.total_counts")%in% colnames(data_plot))){
      g <- do_scatter_plot(data_plot,x="rna.total_counts",y="prot.total_counts", facet=sc)
      ggsave(g, filename=file.path(outpath, paste0("scatter_", sc, "-nUMI_v_rna-nUMI.png")), type="cairo",
              width= 2*ncols, height=2*nrows, dpi=120)
    }
    if (all(c("rna.log1p_total_counts","prot.log1p_total_counts")%in% colnames(data_plot))){
      g <- do_scatter_plot(data_plot,x="rna.log1p_total_counts",y="prot.log1p_total_counts", facet=sc)
      ggsave(g, filename=file.path(outpath, paste0("scatter_", sc, "-log1p_nUMI_v_rna-log1p_nUMI.png")), type="cairo",
              width= 2*ncols, height=2*nrows, dpi=120)
    }
    if (all(c("rna.total_counts","prot.total_counts_isotype")%in% colnames(data_plot))){
      g <-  do_scatter_plot(data_plot,x="rna.total_counts",y="prot.total_counts_isotype", facet=sc)
      ggsave(g, filename=file.path(outpath, paste0("scatter_", sc, "-nUMI_v_prot-counts_isotype.png")), type="cairo",
              width= 2*ncols, height=2*nrows, dpi=120)
    }    
    if (all(c("rna.log1p_total_counts","prot.log1p_total_counts_isotype") %in% colnames(data_plot))){
      g <-  do_scatter_plot(data_plot,x="rna.log1p_total_counts",y="prot.log1p_total_counts_isotype", facet=sc)
      ggsave(g, filename=file.path(outpath, paste0("scatter_", sc, "-log1p_nUMI_v_prot-log1p_counts_isotype.png")), type="cairo",
              width= 2*ncols, height=2*nrows, dpi=120)
    }
    if (all(c("rna.doublet_scores","prot.log1p_total_counts") %in% colnames(data_plot))){
      g <-  do_scatter_plot(data_plot,x="rna.doublet_scores",y="prot.log1p_total_counts", facet=sc)
      ggsave(g, filename=file.path(outpath, paste0("scatter_", sc, "-doublet_scores_v_prot-log1p_nUMI.png")), type="cairo",
              width= 2*ncols, height=2*nrows, dpi=120)
      
    }
    
    
  }
  
}


sprefix <- opt$sampleprefix
print(sprefix)
message ("saving some counts tables for references")


if(opt$prefilter){
  if(all(c("pct_counts_mt", "pct_counts_hb", "n_genes_by_counts", "doublet_scores") %in% colnames(rna_data_plot))){
    f1 <- rna_data_plot %>% 
      dplyr::filter(pct_counts_mt<=20 & pct_counts_hb<=70 &n_genes_by_counts>=100 & doublet_scores<=0.25) %>%
      group_by_at(.vars=c(rna_source_facet)) %>%
      group_by_at(.vars=c(rna_source_facet)) %>%
      summarise(cell.count= n()) %>%
      group_by_at(.vars="sample_id") %>% 
      rename(cell.count_f1=cell.count) 
    
    f2 <- rna_data_plot %>% 
      dplyr::filter(pct_counts_mt<=10 & pct_counts_hb<=50 &n_genes_by_counts>=100 & doublet_scores<=0.25) %>%
      group_by_at(.vars=c(rna_source_facet)) %>%
      summarise(cell.count= n()) %>% 
      group_by_at(.vars="sample_id") %>% 
      rename(cell.count_f2=cell.count) 
    
    f3 <- rna_data_plot %>% 
      dplyr::filter(pct_counts_mt<=5 & pct_counts_hb<=50 &n_genes_by_counts>=100 & n_genes_by_counts<=3000) %>%
      group_by_at(.vars=c(rna_source_facet)) %>%
      summarise(cell.count= n()) %>% 
      group_by_at(.vars="sample_id") %>% 
      rename(cell.count_f3=cell.count) 
    
    baseline <- rna_data_plot %>% 
      group_by_at(.vars=c(rna_source_facet)) %>%
      summarise(cell.count= n()) %>% 
      group_by_at(.vars="sample_id") %>% 
      rename(baseline.counts=cell.count)
    
    info <- merge(merge(merge(f1,f2,by=rna_source_facet,all=TRUE),f3, by=rna_source_facet, all=TRUE), baseline, by=rna_source_facet, all=T) %>%
      mutate(percent_retain_f1 = 100*cell.count_f1/baseline.counts,
             percent_retain_f2 = 100*cell.count_f2/baseline.counts,
             percent_retain_f3 = 100*cell.count_f3/baseline.counts) 
    
    
    
    write.table(info, file=paste0( sprefix,"_threshold_filter.tsv"), col.names=T, row.names=F, sep="\t", quote=F)
    
    explain <- data.frame(qcmetric=c("pct_counts_mt_max","pct_counts_hb_max","n_genes_by_counts_min","doublet_scores_max","n_genes_by_counts_max"),
                          f1=c(20,70,100,0.25,NA),
                          f2=c(10,50,100,0.25,NA),
                          f3=c(5,50,100,NA,3000)) 
    
    write.table(explain, file=paste0(sprefix,"_threshold_filter_explained.tsv"), col.names=T, row.names=F, sep="\t", quote=F)
    
    lab <- c("%Mt <=10 & %HB<70 & minGenes>=100 & scrublet<=0.25",
             "%Mt <=5 & %HB<50 & minGenes>=100 & scrublet<=0.25",
             "%Mt <=5 & %HB<50 & minGenes>=100 & maxGenes<=3000")
    names(lab) <- c("percent_retain_f1","percent_retain_f2","percent_retain_f3") 
    
    g <- info %>% 
      pivot_longer(cols=starts_with("percent"), names_to="filterclass", values_to="percent")%>%
      ggplot(aes(sample_id,percent, fill=filterclass)) +
      geom_bar(stat="identity", position="dodge", color="black") +
      facet_wrap(~filterclass, ncol=1, labeller = labeller(filterclass=lab)) +
      theme_bw()+
      theme(axis.text.x=element_text(size=8,angle=45, hjust=1, vjust=1),
            axis.text.y=element_text(size=13)) +
      scale_fill_manual(values = c("red", "yellow","blue"), 
                        breaks=c("percent_retain_f1","percent_retain_f2","percent_retain_f3"),
                        limits=c("percent_retain_f1","percent_retain_f2","percent_retain_f3")) + 
      coord_cartesian(ylim=c(0,100)) 
    
    ggsave(g, file = paste0(run,"barplot_cellcounts_thresholds_filter.png"), type="cairo", width=9, height=9)
  }
  
}else{
  message("producing files with final counts for cells after filtering")
  
  baseline <- rna_data_plot %>% 
    group_by_at(.vars=c(rna_source_facet)) %>% 
    summarise(cell.count= n()) %>% 
    group_by_at(.vars="sample_id") 
  
  write.table(baseline, file=paste0( sprefix,"_filtered_data.tsv"), col.names=T, row.names=F, sep="\t", quote=F)
  
  g <- baseline %>% 
    ggplot(aes_string(x="sample_id", y="cell.count")) +
    geom_bar(stat="identity", position="dodge", color="black") +
    theme_bw()+
    theme(axis.text.x=element_text(size=8,angle=45, hjust=1.05, vjust=0.95),
          axis.text.y=element_text(size=13)) 
  ggsave(g, file = file.path(run, "rna","barplot_cellcounts_filtered_data.png"), type="cairo", width=9, height=6)
  
}


message("done")
