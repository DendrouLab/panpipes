"""
Fetching data for plotting
CRG 2020-06-19
"""

import scanpy as sc
from anndata import AnnData
import muon as mu
import pandas as pd
import argparse
import os
import re
sc.settings.autoshow = False
from panpipes.funcs.io import read_yaml
import matplotlib
matplotlib.use('agg')


# from panpipes.funcs.processing import check_for_bool
# from panpipes.funcs.io import read_anndata, write_anndata

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


def main(adata, mod, layer, group_col, mf, n):
    L.debug("check layers")
    if layer != 'X':
        adata.X = adata.layers[layer]
    # run pca if not done already.( for dendorgram calc)
    if("X_pca" in adata.obsm.keys()):
        sc.pp.pca(adata)
    # read file
    df = mf.groupby(group_col).apply(lambda x: x.nlargest(n, ['scores'])).reset_index(drop=True)
    marker_list={str(k): list(v) for k,v in df.groupby(group_col)["gene"]}
    # add cluseter col to obs
    adata.obs[group_col] = adata.obs[group_col].astype('str').astype('category')
    L.debug("do plots")
    sc.pl.stacked_violin(adata,
                marker_list,
                groupby=group_col,
                save=  '_top_markers'+ mod +'.png',
                figsize=(24, 5))
    sc.pl.matrixplot(adata,
                    marker_list,
                    groupby=group_col,
                    save=  '_top_markers'+ mod +'.png',
                    figsize=(24, 5))
    sc.pl.dotplot(adata,
                marker_list,
                groupby=group_col,
                save=  '_top_markers'+ mod +'.png',
                figsize=(24, 5))
    sc.pl.heatmap(adata,
                marker_list,
                groupby=group_col,
                save=  '_top_markers'+ mod +'.png',
                figsize=(24, 5))


L.debug("load data")
mdata = mu.read(args.infile)


L.debug("load marker file")
mf = pd.read_csv(args.marker_file, sep='\t' )
mf[args.group_col] = mf['cluster'].astype('str').astype('category')
# get layer
try:
    layers = read_yaml(args.layer)
except AttributeError:
    layers = args.layer

if type(mdata) is mu.MuData and args.modality is None and type(layers) is not dict:
    sys.exit('if inputting a mudata object, layers must be in a dictionary')
    

if type(mdata) is AnnData:
    adata = mdata
    main(adata,mod=args.modality, layer=args.layer, n=args.n)
elif args.modality is not None:
    adata = mdata[args.modality]
    if type(layers) is dict:
        ll = layers[args.modality]
    else:
        ll = layers
    main(adata,
         mod=args.modality, 
         layer = ll, 
         mf=mf, 
         group_col = args.group_col,
         n=int(args.n))
else:
    # we have multimodal object, and some kind of multimodal clustering output
    for mod in mf['mod'].unique():
        print(mod)
        mf_sub = mf[mf['mod'] == mod]
        mdata[mod].obs[args.group_col] = mdata.obs.loc[mdata[mod].obs_names,args.group_col].astype('category')
        mdata.update_obs()
        main(mdata[mod], mod=mod, layer = layers[mod], group_col = args.group_col, mf=mf_sub, n=int(args.n))


