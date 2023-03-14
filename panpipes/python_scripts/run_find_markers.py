import scanpy as sc
import argparse
import pandas as pd
import re
import numpy as np
from scipy.sparse import issparse
from panpipes.funcs.processing import check_for_bool
from panpipes.funcs.io import read_anndata
from panpipes.funcs.scmethods import pseudo_seurat

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
                    default="rnapbmc3k_nn30_drharmony_metcosine_neighbors.h5ad",
                    help="file name, format: .h5ad")
parser.add_argument('--use_muon',
                    default=False,
                    help='')
parser.add_argument('--modality',
                    default="rna",
                    help='')
parser.add_argument("--outputfile",
                    default="rnacluster1_markers.txt",
                    help="file name, format: .txt")
parser.add_argument("--cluster_file",
                    default="rnapbmc3k_nn30_res0.2_algleiden_cluster.txt",
                    help="file name, format: .h5ad")
parser.add_argument("--cluster_value", type=str,
                    default='1', 
                    help="the ID of the cluster that we want to look for markers in.")
parser.add_argument("--testuse", 
                    default='wilcoxon', 
                    help="wilcoxon, t-test or logreg")
parser.add_argument("--pseudo_seurat",
                    default=False, type=check_for_bool,
                    help="apply seurat like filtering before marker test")
parser.add_argument("--minpct", 
                    default=0.1,
                    help="minimum fraction of cells expressing gene")
parser.add_argument("--mindiffpct", type=str, 
                    default='-inf',
                    help="minimum fraction difference between cluster and other cells. Setting not recommended.")
parser.add_argument("--threshuse",  
                    default=0.25,
                    help="testing limited to genes with this (log scale) difference in mean expression level.")
parser.add_argument("--mincells", 
                    default=3,
                    help="minimum number of cells required (applies to cluster and to the other cells)")


args, opt = parser.parse_known_args()
use_muon = args.use_muon


# script -------------------------------

# check for options

L.info("Running with options")
L.info(args)
L.info('\n')
# read data
adata = read_anndata(args.infile, use_muon=use_muon, modality=args.modality)

# adata = sc.read_h5ad("../../data/scanpy_pipe_test/pbmc3k_nn30_drharmony_metcosine_neighbors.h5ad")

# load clusters
df = pd.read_csv(args.cluster_file, sep="\t", index_col=0)
# df = pd.read_csv("../../data/scanpy_pipe_test/pbmc3k_res0.6_nn30_algleiden_cluster.txt", sep="\t", index_col=0 )
#
adata_shape = adata.shape


# add clusters to adata.obs
adata.obs = pd.merge(adata.obs, df, how="left", left_index=True, right_index=True )
# check check we have not lost cells in merge
if adata.shape != adata_shape:
    L.info("some cells lost in merge, not all cells have a cluster?")



# how many unique clusters?

# 1. set idents
adata.obs["idents"] = ["1" if cc==int(args.cluster_value) else "0" for cc in adata.obs["clusters"]]
adata.obs["idents"] = adata.obs["idents"].astype("category")

# check we have enough cells
if args.mincells == "None":
    min_cells = 3
else:
    min_cells = int(args.mincells)

n_cluster = sum([cell is '1' for cell in adata.obs['idents']])
n_other = sum([cell is '0' for cell in adata.obs['idents']])
if n_cluster <= int(args.mincells) | n_other <= int(args.mincells):
    sys.exit("not enough cells in cluster, min cells: %i" % int(args.mincells))

L.info("Cluster %s number of cells: %i" % (args.cluster_value, n_cluster))
L.info("Other number of cells: %i\n" % (n_other))

# there is some weird difference between dense and sparse array?
if issparse(adata.X):
    bool_list = (adata.X.sum(axis=0) > 0).tolist()[0]
else:
    bool_list = (adata.X.sum(axis=0) > 0).tolist()
adata = adata[:, bool_list]

if isinstance(args.pseudo_seurat, str):
    if args.pseudo_seurat == 'True':
        run_pseudo_seurat = True
    elif args.pseudo_seurat == 'False':
        run_pseudo_seurat = False
    else:
        L.info(args.pseudo_seurat)
        sys.exit("spelling_error: please specify --pseudoseurat as True or False")
elif type(args.pseudo_seurat) == 'bool':
    run_pseudo_seurat = args.pseudo_seurat
else:
    sys.exit("type error %s: please specify --pseudo_seurat as True or False" % type(args.pseudo_seurat))


if run_pseudo_seurat:
    L.info("running pseudo seurat with params: \nminpct: %s \nmindiffpct: %s \nlogFCthreshold: %s\n\n" % (args.minpct, args.mindiffpct, args.threshuse))

    filter_stats = pseudo_seurat(adata,
                                 arg_minpct=float(args.minpct),
                                 arg_mindiffpct=float(args.mindiffpct),
                                 arg_logfcdiff=float(args.threshuse))
    L.info("number of genes remaining after filtering:  %i\n" % filter_stats['background'].sum())
    # adata = adata.raw.to_adata().copy()
    adata = adata[:, filter_stats['background'].tolist()]
#     adata.raw = adata
    L.info(adata)
    L.info('\n')
else:
    # get backgrounds
    pass
    # write out background after rank_gene_groups


# now we can actually run find markers,
# run on a subset of genes
sc.tl.rank_genes_groups(adata, groups=["1"], groupby="idents", reference="0",
                        method="wilcoxon", n_genes=float("inf"), corr_method="bonferroni", layer="logged_counts")

markers = sc.get.rank_genes_groups_df(adata, group="1")

# merge markers with filter_stats_df, so that we can see all the scroes for gsea
#

# filter positive only and adjusted pval < 0.05
# markers = markers[(markers.logfoldchanges > 0) & (markers.pvals_adj < 0.05)]

if sc.__version__ < "1.7.0":
    markers.columns = ['scores', 'gene', 'avg_logFC', 'pvals', 'p.adj.bonferroni']
elif sc.__version__ == "1.7.0":
    markers.columns = ['gene', 'scores', 'avg_logFC', 'pvals', 'p.adj.bonferroni']

markers.head()
L.info("Number of Markers for cluster %s: %i" % (args.cluster_value, markers.shape[0]))
L.info("Rank gene groups complete")


L.info("saving markers to:" + args.outputfile)
markers.to_csv(args.outputfile, sep="\t", index=False)

if run_pseudo_seurat:
    fname = re.sub('_markers.txt', '_background.txt', args.outputfile)
    L.info("saving background to:" + fname)
    filter_stats.to_csv(fname, sep="\t", index=False)


# quick do some plotting!

L.info("Done")

