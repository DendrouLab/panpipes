import scanpy as sc
from muon import read, MuData
from anndata import AnnData
import argparse
import pandas as pd
import numpy as np
from scipy.sparse import issparse
from panpipes.funcs.processing import check_for_bool
from panpipes.funcs.scmethods import find_all_markers_pseudo_seurat


import sys
import logging
L = logging.getLogger()
L.setLevel(logging.INFO)
log_handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s: %(levelname)s - %(message)s')
log_handler.setFormatter(formatter)
L.addHandler(log_handler)


parser = argparse.ArgumentParser()
parser.add_argument("--infile",
                    default="anndata_log1p.h5ad",
                    help="file name, format: .h5ad")
parser.add_argument('--modality',
                    default=None,
                    help='')
parser.add_argument('--layer',
                    default="logged_counts",
                    help='')
parser.add_argument("--output_file_prefix",
                    default="./",
                    help="output directory/markers")
parser.add_argument("--cluster_file",
                    default=None,
                    help="file name, format: .h5ad")
parser.add_argument("--testuse", default='wilcoxon', help="wilcoxon, t-test or logreg")
parser.add_argument("--pseudo_seurat", default=False, type=check_for_bool,
                    help="apply seurat like filtering before marker test")
parser.add_argument("--minpct", type=str, default='0.1', help="minimum fraction of cells expressing gene")
parser.add_argument("--mindiffpct", type=str, default='-inf',
                    help="minimum fraction difference between cluster and other cells. Setting not recommended.")
parser.add_argument("--threshuse",  type=str, default='0.25',
                    help="testing limited to genes with this (log scale) difference in mean expression level.")
parser.add_argument("--mincells", type=str, default='3',
                    help="minimum number of cells required (applies to cluster and to the other cells)")


args, opt = parser.parse_known_args()


# script -------------------------------

def check_log1p_dict(adata):
    # temp code to get around this bug
    # https://github.com/scverse/muon/issues/40
    if 'log1p' in adata.uns.keys():
        if not adata.uns['log1p']:
            L.debug( "log1p dict is Empty, refilling dict")
            adata.uns['log1p'] = {'base': None}


def load_and_merge_clusters(adata, cluster_file):
    # load clusterss
    df = pd.read_csv(cluster_file, sep="\t", index_col=0)
    adata_shape = adata.shape
    # add clusters to adata.obs
    adata.obs = pd.merge(adata.obs, df, how="left", left_index=True, right_index=True)
    # check check we have not lost cells in merge
    if adata.shape != adata_shape:
        L.info("some cells lost in merge, not all cells have a cluster?")
    L.debug(adata.shape)


def get_header():
    # write out the correct header before we start, need to werite it out before, because we use append mode for the multimodal find markers
    # and we don't wantto end up with a new header half
    if sc.__version__ < "1.7.0":
        markers_columns = ['cluster', 'scores', 'gene', 'avg_logFC', 'pvals', 'p.adj.bonferroni']
    elif sc.__version__ >= "1.7.0":
        markers_columns = ['cluster', 'gene', 'scores', 'avg_logFC', 'pvals', 'p.adj.bonferroni']
    return markers_columns


def filter_zero_count_genes(adata, layer=None):
    if layer is None:
        if issparse(adata.X):
            bool_list = (adata.X.sum(axis=0) > 0).tolist()[0]
        else:
            bool_list = (adata.X.sum(axis=0) > 0).tolist()
    else:
        if issparse(adata.layers[layer]):
            bool_list = (adata.layers[layer].sum(axis=0) > 0).tolist()[0]
        else:
            bool_list = (adata.layers[layer].sum(axis=0) > 0).tolist()
    adata = adata[:, bool_list]
    L.debug(adata.shape)


def run_clustering(adata, 
                   cluster_col=None, 
                   cluster_groups='all', 
                   layer=None, 
                   methoduse=None, 
                   pseudo_seurat=False):
    if type(cluster_groups) is list:
        # make sure the clusters are strings else rank_gene_groups crashes
        cluster_groups = [str(x) for x in cluster_groups]
    adata.obs[cluster_col] = adata.obs[cluster_col].astype('str').astype('category')
    if pseudo_seurat:
        L.info('running rank gene groups with pseudo-seurat method')
        markers, filter_stats = find_all_markers_pseudo_seurat(adata,   
                                        groups=cluster_groups, 
                                        groupby=cluster_col, 
                                        method=methoduse,
                                        n_genes=float("inf"), 
                                        corr_method="bonferroni",
                                        layer=layer,
                                        arg_minpct=0.1,
                                        arg_mindiffpct=-float("inf"), 
                                        arg_logfcdiff=0.25)
        markers = markers.reset_index().drop(columns='level_1').rename(columns={'level_0':'gene'})
    else:
        # find markers
        L.info('running rank gene groups with standard scanpy implementation')
        sc.tl.rank_genes_groups(adata, 
                                groups=cluster_groups, 
                                groupby=cluster_col, 
                                method=methoduse,
                                n_genes=float("inf"), 
                                corr_method="bonferroni",
                                layer=layer)
        markers = sc.get.rank_genes_groups_df(adata, group=None)
        L.info("Rank gene groups complete")
        filter_stats = None
    return markers, filter_stats


def main(adata, 
         cluster_file,
         mod, 
         mincells=None,
         layer=None,
         testuse=None,
         pseudo_seurat=False,
         output_file_prefix="markers") :
    # check the X slot actually contains data
    if adata.X.shape[1] == 0:
        L.info("exiting because adata.X contains no values")
        sys.exit("exiting because adata.X contains no values")
    # prevent log1p base error.
    check_log1p_dict(adata)
    # load clusters
    load_and_merge_clusters(adata, cluster_file)
    # filter out genes with zero counts
    filter_zero_count_genes(adata, layer=layer)
    # filter out clusters with less than min cells
    if mincells is not None:
        min_cell_df = adata.obs['clusters'].value_counts() 
        if any(min_cell_df < mincells):
            exclude = min_cell_df[min_cell_df < mincells].index
            L.info('exluding clusters for haveing too few cells: %s' % str(list(exclude)))
            clust_vals = list(set(min_cell_df[min_cell_df >= mincells].index))
        else:
            clust_vals = list(set(adata.obs["clusters"]))
    all_markers, all_filter_stats = run_clustering(adata,
                                                   cluster_groups=clust_vals,
                                                   cluster_col="clusters",
                                        layer=layer, 
                                        methoduse=testuse, 
                                        pseudo_seurat=pseudo_seurat)
    # make sure all files have consitent headers
    all_markers.columns = get_header()
    all_markers['mod'] = mod
    all_markers.to_csv(output_file_prefix + ".txt", index=False, sep='\t')
    # subset to significant markers
    signif_markers = all_markers[(all_markers["p.adj.bonferroni"] < 0.05) &( all_markers.avg_logFC > 0)]
    if signif_markers.shape[0] !=0:
        signif_markers.to_csv(output_file_prefix + mod + "_signif.txt", index=False, sep='\t')
        # write out to excel
        excel_file_top = output_file_prefix + mod + "_signif.xlsx"
        with pd.ExcelWriter(excel_file_top) as writer:
            for xx in signif_markers.cluster.unique():
                df_sub = signif_markers[signif_markers.cluster == xx]
                if df_sub.shape[0] !=0:
                    df_sub.to_excel(writer, sheet_name="cluster" + str(xx), index=False)
    #  when relevent write out filter stats
    if pseudo_seurat is True:
        all_filter_stats = all_filter_stats.reset_index().rename(columns={'index':'cluster'})
        all_filter_stats.to_csv(output_file_prefix + "_filter_stats_" + mod + ".txt", index=False, sep="\t")



L.info("Running with options")
L.info(args)
L.info('\n')

# read data
mdata = read(args.infile)
    

if type(mdata) is AnnData:
    adata = mdata
    # main function only does rank_gene_groups on X, so 
elif type(mdata) is MuData and args.modality is not None:
    adata = mdata[args.modality]    
else:
    sys.exit('if inputting a mudata object, you need to specify a modality')
    


main(adata, 
     mod=args.modality,
     cluster_file=args.cluster_file,
    mincells=int(args.mincells),
    layer=args.layer,
    testuse=args.testuse,
    pseudo_seurat=args.pseudo_seurat,
    output_file_prefix=args.output_file_prefix)

L.info("done")

