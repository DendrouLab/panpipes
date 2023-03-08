import scanpy as sc
from muon import read, MuData
from anndata import AnnData
import argparse
import pandas as pd
import re
import numpy as np
from scipy.sparse import issparse
from itertools import compress

from panpipes.funcs.processing import check_for_bool
from panpipes.funcs.scmethods import pseudo_seurat
from panpipes.funcs.io import read_yaml

import sys
import logging
L = logging.getLogger()
L.setLevel(logging.DEBUG)
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
parser.add_argument("--pseudo_seurat", default=False, help="apply seurat like filtering before marker test")
parser.add_argument("--minpct", type=str, default='0.1', help="minimum fraction of cells expressing gene")
parser.add_argument("--mindiffpct", type=str, default='-inf',
                    help="minimum fraction difference between cluster and other cells. Setting not recommended.")
parser.add_argument("--threshuse",  type=str, default='0.25',
                    help="testing limited to genes with this (log scale) difference in mean expression level.")
parser.add_argument("--mincells", type=str, default='3',
                    help="minimum number of cells required (applies to cluster and to the other cells)")
parser.add_argument("--use_dense", type=str, default="False")

args, opt = parser.parse_known_args()
# args = argparse.Namespace(infile='prot/ncomps3_nneigh30_meteuclidean_dir/neighbors.h5ad', modality='prot', layer='rna', output_file_prefix='prot/ncomps3_nneigh30_meteuclidean_dir/algleiden_res1/markers', cluster_file='prot/ncomps3_nneigh30_meteuclidean_dir/algleiden_res1/clusters.txt.gz', testuse='wilcoxon', pseudo_seurat='False', minpct='0.1', mindiffpct='-inf', threshuse='0.25', mincells='10', use_dense='False')


# script -------------------------------

def main(adata, args, mod="", marker_columns=[]):
    # check the X slot actually contains data
    if adata.X.shape[1] == 0:
        L.info("skipping %s because adata.X contains no values" % mod)
        return None
    
    # temp code to get around this bug
    # https://github.com/scverse/muon/issues/40
    if 'log1p' in adata.uns.keys():
        if not adata.uns['log1p']:
            L.debug( "log1p dict is Empty, refilling dict")
            adata.uns['log1p'] = {'base': None}

    # load clusterss
    df = pd.read_csv(args.cluster_file, sep="\t", index_col=0)
    adata_shape = adata.shape

    # add clusters to adata.obs
    adata.obs = pd.merge(adata.obs, df, how="left", left_index=True, right_index=True)
    # check check we have not lost cells in merge
    if adata.shape != adata_shape:
        L.info("some cells lost in merge, not all cells have a cluster?")
    L.debug(adata.shape)
    if issparse(adata.X):
        bool_list = (adata.X.sum(axis=0) > 0).tolist()[0]
    else:
        bool_list = (adata.X.sum(axis=0) > 0).tolist()

    adata = adata[:, bool_list]
    L.debug(adata.shape)
    run_pseudo_seurat = check_for_bool(args.pseudo_seurat)
    use_dense = check_for_bool(args.use_dense)

    clust_vals = set(adata.obs["clusters"])
    markers_dict = {}
    filter_stats_dict = {}

    for cv in clust_vals:
        L.info(cv)
        # 1. set idents
        adata.obs["idents"] = ["1" if cc == cv else "0" for cc in adata.obs["clusters"]]
        adata.obs["idents"] = adata.obs["idents"].astype("category")
        # check we have enough cells
        if args.mincells == "None":
            min_cells = 3
        else:
            min_cells = int(args.mincells)
        n_cluster = sum([cell == '1' for cell in adata.obs['idents']])
        n_other = sum([cell == '0' for cell in adata.obs['idents']])
        if n_cluster <= int(args.mincells) | n_other <= int(args.mincells):
            sys.exit("not enough cells in cluster, min cells: %i" % int(args.mincells))
        L.info("Cluster %s number of cells: %i" % (cv, n_cluster))
        L.info("Other number of cells: %i\n" % n_other)
        # adata = adata.raw.to_adata().copy()
        if run_pseudo_seurat is True:
            filter_stats = pseudo_seurat(adata, use_dense=use_dense)
            L.info("number of genes remaining after filtering:  %i\n" % filter_stats['background'].sum())
            adata_rg = adata[:, filter_stats['background'].tolist()].copy()

            sc.tl.rank_genes_groups(adata_rg, layer=None,
            groups=["1"], groupby="idents", reference="0",
            method="wilcoxon", n_genes=float("inf"), corr_method="bonferroni")

            markers = sc.get.rank_genes_groups_df(adata_rg, group="1",)
            # remove adata from mem
            adata_rg = None
        else: 
            sc.tl.rank_genes_groups(adata, groups=["1"], groupby="idents", reference="0", 
                                    method="t-test_overestim_var", 
                                    n_genes=float("inf"), 
                                    corr_method="bonferroni", layer=None)
            markers = sc.get.rank_genes_groups_df(adata, group="1")
        # filter positive only and adjusted pval < 0.05
        # markers = markers[(markers.logfoldchanges > 0) & (markers.pvals_adj < 0.05)]

        markers.head()
        L.info("Rank gene groups complete")
        markers_dict[cv] = markers
        if pseudo_seurat is True:
            filter_stats_dict[cv] = filter_stats
    
    
    all_markers = pd.concat(markers_dict, axis=0, keys=markers_dict.keys())
    all_markers.columns = marker_columns[1:6]
    all_markers = all_markers.reset_index().rename(columns={'level_0':'cluster'}).drop(columns="level_1")
    all_markers['mod'] = mod
    all_markers.to_csv(args.output_file_prefix + ".txt", index=False, sep='\t', header=None, mode="a")

    signif_markers = all_markers[(all_markers["p.adj.bonferroni"] < 0.05) &( all_markers.avg_logFC > 0)]
    if signif_markers.shape[0] !=0:
        signif_markers.to_csv(args.output_file_prefix + mod + "_signif.txt", index=False, sep='\t')
        L.debug("test_debug")
        excel_file_top = args.output_file_prefix + mod + "_signif.xlsx"
        with pd.ExcelWriter(excel_file_top) as writer:
            for xx in signif_markers.cluster.unique():
                L.debug(xx)
                df_sub = signif_markers[signif_markers.cluster == xx]
                if df_sub.shape[0] !=0:
                    df_sub.to_excel(writer, sheet_name="cluster" + str(xx), index=False)

    if pseudo_seurat is True:
        all_filter_stats = pd.concat( filter_stats_dict, axis=0, keys= filter_stats_dict.keys(), ignore_index=True)
        all_filter_stats = all_filter_stats.reset_index().rename(columns={'index':'cluster'})
        all_filter_stats.to_csv(args.output_file_prefix + "_filter_stats_" + mod + ".txt", index=False, sep="\t")



L.info("Running with options")
L.info(args)
L.info('\n')

# get layer
try:
    layers = read_yaml(args.layer)
except AttributeError:
    layers = args.layer

# read data
mdata = read(args.infile)
    
if type(mdata) is MuData and args.modality is None and type(layers) is not dict:
    sys.exit('if inputting a mudata object, layers must be in a dictionary')
    
# write out the correct header before we start, need to werite it out before, because we use append mode for the multimodal find markers
# and we don't wantto end up with a new header half
if sc.__version__ < "1.7.0":
    markers_columns = ['cluster', 'scores', 'gene', 'avg_logFC', 'pvals', 'p.adj.bonferroni', 'mod']
elif sc.__version__ >= "1.7.0":
    markers_columns = ['cluster', 'gene', 'scores', 'avg_logFC', 'pvals', 'p.adj.bonferroni', 'mod']
pd.DataFrame(columns=markers_columns).to_csv(args.output_file_prefix + ".txt", index=False, sep='\t', mode="w")


if type(mdata) is AnnData:
    adata = mdata
    if layers != "X" :
        adata.X = adata.layers[layers[args.modality]]
    main(adata, args,marker_columns=markers_columns, mod=args.modality)
elif args.modality is not None:
    adata = mdata[args.modality]
    if type(layers) is dict:
        ll = layers[args.modality]
    else:
        ll = layers
    adata.X = adata.layers[ll]
    main(adata, args, marker_columns=markers_columns, mod=args.modality)
else:
    # we have multimodal object, and some kind of multimodal clustering output
    for mod in mdata.mod.keys():
        if type(layers) is dict:
            ll = layers[args.modality]
        else:
            ll = layers
        mdata[mod].X = mdata[mod].layers[ll]
        main(mdata[mod], args, mod=mod, marker_columns=markers_columns)


