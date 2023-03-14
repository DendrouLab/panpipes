"""
Plotting umaps of a small list of custom genes
CRG 2020-06-19
"""
import scanpy as sc
import muon as mu
from anndata import AnnData
import pandas as pd
import argparse
from panpipes.funcs.io import read_yaml
from itertools import product

sc.settings.autoshow = False
import os
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
                    default="pbmc3k.h5mu",
                    help="file name, format: .h5mu/.h5ad")
parser.add_argument("--modalities",
                    default="rna",
                    help="modality keys containing data to be plotted, comma spearated list e.g. rna,prot default 'rna'")
parser.add_argument("--layers",
                    default="X",
                    help="layer containing data to be plotted, default=X")
parser.add_argument("--basis_dict",
                    default="X_umap",
                    help="basis in obsm to plot")
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

def main(adata, mod, layer_choice, df, basis):
    for gc in df['group'].unique():
        # define file name
        if layer_choice is None or layer_choice =="X":
            layer_choice = None
            layer_string = ""
        else:
            layer_string = layer_choice
        fname_prefix = "_".join(["_" + mod, layer_string, gc])
        # get features
        fetches = df[df['group'] == gc]['feature']
        plot_features = [gg for gg in fetches if gg in adata.var_names]
        plot_features = list(set(plot_features))
        L.info("plotting")
        L.info(plot_features)
        sc.settings.figdir  = os.path.join(args.base_figure_dir, mod)
        mu.pl.embedding(adata, basis=basis, layer=layer_choice, color=plot_features, save = fname_prefix + ".png")


L.debug("load data")
mdata = mu.read(args.infile)
modalities = args.modalities.split(',')

L.debug("load marker file")
df = pd.read_csv(args.marker_file )

# get layer
try:
    layers = read_yaml(args.layers)
except AttributeError:
    # this assumes that we have tried to parse a dict and nstead found a string
    # there is probably a better solution
    layers = args.layers

# get bases
try:
    basis_dict = read_yaml(args.basis_dict)
except AttributeError:
    # this assumes that we have tried to parse a dict and nstead found a string
    # there is probably a better solution
    basis_dict = args.basis_dict


if type(mdata) is AnnData:
    adata = mdata
    basis= basis_dict
    main(adata, layer_choice=args.layers,  df = df, basis=basis)
else:
    # we have multimodal object
    for mod in modalities:
        print(mod)
        df_sub = df[df['mod'] == mod]
        mdata.update_obs()
        try:
            ll = layers[mod]
        except KeyError:
            ll = [None]
        if mod in basis_dict.keys():
            bb= basis_dict[mod]
        else:
            bb = []
        if len(bb) > 0 :
            for basis, layer in product(bb, ll):
                print(basis,layer)
                main(adata=mdata[mod], 
                    mod=mod,
                    layer_choice = layer,
                    df=df_sub,
                    basis=basis)

L.info("Done")

