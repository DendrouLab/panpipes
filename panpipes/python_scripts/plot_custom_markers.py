"""
Fetching data for plotting
CRG 2020-06-19
"""

import scanpy as sc
import muon as mu
from anndata import AnnData

import pandas as pd
import argparse
from itertools import product
sc.settings.autoshow = False
import os
import matplotlib
matplotlib.use('agg')
import re
from panpipes.funcs.io import  read_yaml

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
                    default="pbmc3k.h5ad",
                    help="file name, format: .h5ad")
parser.add_argument("--modalities",
                    default="rna",
                    help="modality keys containing data to be plotted, comma spearated list e.g. rna,prot default 'rna'")
parser.add_argument("--layers",
                    default="X",
                    help="layer containing data to be plotted, default=X")
parser.add_argument("--group_cols", default=["sample_id"],
                help="grouping variables", nargs='+',)
parser.add_argument("--marker_file",
                    default="Marker_Lists-myeloid.csv",
                    help="comma separateed list of marker_list csvs, one column of gene_ids")
parser.add_argument("--base_figure_dir", default="figures/",
                    help="figures path")


args, opt = parser.parse_known_args()
L.info("running with:")
L.info(args)
sc.settings.figdir = args.base_figure_dir

# ---- script

def main(adata, mod, df, grouping_var, pfx, layer_choice=None):
    for gc in df['group'].unique():
        # define file name
        if layer_choice is None or layer_choice =="X":
            layer_choice = None
            layer_string = ""
        else:
            layer_string = layer_choice
        # get features
        fetches = df[df['group'] == gc]['feature']
        plot_features = [gg for gg in fetches if gg in adata.var_names]
        plot_features = list(set(plot_features))
        # do PCA if it is missing, (prerequisite for dendrogram)
        if "X_pca" not in adata.obsm.keys():
            sc.pp.pca(adata)
        use_dendrogram=True
        if len(adata.obs[grouping_var].unique()) > 2:
            try:
                sc.tl.dendrogram(adata, grouping_var, use_rep="X_pca",linkage_method="average")
                use_dendrogram=True
            except ValueError: 
                use_dendrogram=False

        sc.settings.figdir  = os.path.join(args.base_figure_dir, mod ,  re.sub(":", "_", grouping_var))
        fname_prefix = "_".join([layer_string, pfx, gc ])
        fname_prefix = re.sub(":", "_", fname_prefix)
        sc.pl.dotplot(adata,
                        var_names=plot_features,
                        groupby=grouping_var,
                        layer=layer_choice,
                        dendrogram=use_dendrogram,
                        save=fname_prefix + '.png',
                        figsize=(24, 5))
        sc.pl.matrixplot(adata,
                        var_names=plot_features,
                        groupby=grouping_var,
                        dendrogram=use_dendrogram,
                        layer=layer_choice,
                        save=fname_prefix + '.png',
                        figsize=(24, 5))

L.debug("load data")
mdata = mu.read(args.infile)
modalities = args.modalities.split(',')

L.debug("load marker file")
df = pd.read_csv(args.marker_file )

pfx = re.sub(".csv", "csv", os.path.basename(args.marker_file))
# write out the features that are not found in adata var (due to filtering, or incorrect name)
not_found = [gg for gg in df['feature'] if gg not in mdata.var_names]
not_found_file = re.sub(".csv", "_features_not_found.txt", os.path.basename(args.marker_file))
if len(not_found) > 0 :
    with open(not_found_file, 'w') as f:
        for item in not_found:
            f.write("%s\n" % item)
            
# get layer
try:
    layers = read_yaml(args.layers)
except AttributeError:
    # this assumes that we have tried to parse a dict and nstead found a string
    # there is probably a better solution
    layers = args.layers

group_vars = args.group_cols

if type(mdata) is AnnData:
    adata = mdata
    for gv in group_vars:
        adata.obs[gv] = adata.obs[gv].astype('category')
        main(adata, layer_choice=args.layers, group = gv, pfx = pfx, df = df)
else:
    # we have multimodal object
    for mod in modalities:
        print(mod)
        df_sub = df[df['mod'] == mod]
        for gv in group_vars:
            mdata[mod].obs[gv] = mdata.obs.loc[mdata[mod].obs_names,gv].astype('category')
        mdata.update_obs()
        try:
            ll = layers[mod]
        except KeyError:
            ll = [None]
        if len(group_vars) > 0 and ll is not None:
            for gv, layer in product(group_vars, ll):
                print(gv, layer)
                main(adata=mdata[mod], 
                    mod=mod,
                    layer_choice = layer,
                    grouping_var=gv, 
                    pfx = pfx,
                    df=df_sub)

L.info("Done")

