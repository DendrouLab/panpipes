"""
Fetching data for plotting
CRG 2020-06-19
"""

import scanpy as sc
from anndata import AnnData
import muon as mu
import pandas as pd
import argparse

sc.settings.autoshow = False

import matplotlib
matplotlib.use('agg')

import sys
import logging
L = logging.getLogger()
L.setLevel(logging.INFO)
log_handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s: %(levelname)s - %(message)s')
log_handler.setFormatter(formatter)
L.addHandler(log_handler)

# parse argumetns
parser = argparse.ArgumentParser()
parser.add_argument("--infile",
                    default="mdata.h5mu, full object, not filtered by HVG",
                    help="file name, format: .h5mu")
parser.add_argument('--modality',
                    default=None,
                    help='')
parser.add_argument("--layer",
                    default="logged_counts",
                    help="layer containing data to be plotted, default=X")
parser.add_argument("--group_col",
                    default="cluster",
                    help="column in obs containing groups for plots,")
parser.add_argument("--marker_file",
                    default="",
                    help="markers from sc.rank_genes_group() as txt(.gz)")
parser.add_argument("--figure_prefix", default="figures/markers_",
                    help="figures path")
parser.add_argument("--n", default=None,
                    help="number of genes to plot per cluster, default=None will plot all genes in your marker_file")



args, opt = parser.parse_known_args()
L.info(args)

sc.settings.figdir = args.figure_prefix

# script

def calc_dendrogram(adata, group_col):
    if len(adata.obs[group_col].cat.categories) > 5:
        incl_dendrogram = True
        if "X_pca" not in adata.obsm.keys():
            sc.pp.pca(adata)
        L.info("calculating dendrogram")
        try:
            sc.tl.dendrogram(adata, groupby=group_col, use_rep="X_pca")
        except ValueError:
            L.info("cannot calculate dendrogram")
            incl_dendrogram = False
    else:
        incl_dendrogram = False
    return incl_dendrogram


def do_plots(adata, mod, group_col, mf, n=10, layer=None):
    L.debug("check layers")
    # get markers for plotting
    mf = mf[mf['avg_logFC'] > 0]
    df = mf.groupby(group_col).apply(lambda x: x.nlargest(n, ['scores'])).reset_index(drop=True)
    marker_list={str(k): list(v) for k,v in df.groupby(group_col)["gene"]}
    # add cluseter col to obs
    # check whether a dendrogram is computed/
    incl_dendrogram = calc_dendrogram(adata, group_col)
    L.info("start plotting")
    sc.pl.stacked_violin(adata,
                marker_list,
                groupby=group_col,
                layer=layer,
                save = '_top_markers'+ mod +'.png',
                dendrogram=incl_dendrogram,
                # figsize=(24, 5)
                )
    sc.pl.matrixplot(adata,
                    marker_list,
                    groupby=group_col,
                    save=  '_top_markers'+ mod +'.png',
                    dendrogram=incl_dendrogram,
                    figsize=(24, 5))
    sc.pl.dotplot(adata,
                marker_list,
                groupby=group_col,
                save=  '_top_markers'+ mod +'.png',
                dendrogram=incl_dendrogram,
                figsize=(24, 5))
    sc.pl.heatmap(adata,
                marker_list,
                groupby=group_col,
                save=  '_top_markers'+ mod +'.png',
                dendrogram=incl_dendrogram,
                # figsize=(24, 5)
                )


# read data
mdata = mu.read(args.infile)

if type(mdata) is AnnData:
    adata = mdata
    # main function only does rank_gene_groups on X, so 
elif type(mdata) is mu.MuData and args.modality is not None:
    adata = mdata[args.modality]    
else:
    sys.exit('if inputting a mudata object, you need to specify a modality')
    

L.info("load marker file")
mf = pd.read_csv(args.marker_file, sep='\t' )
mf[args.group_col] = mf['cluster'].astype('category')

# get layer
adata.obs[args.group_col] = mdata.obs[args.group_col]
do_plots(adata,
         mod=args.modality, 
         group_col=args.group_col,
         mf=mf,
         layer=args.layer, 
         n=int(args.n))

L.info("Done")

