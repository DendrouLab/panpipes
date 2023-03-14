"""
Scanpy inbuilt plots for metadata variables
CRG 2020-06-19
"""

import scanpy as sc
import muon as mu
from anndata import AnnData
import pandas as pd
import argparse
from panpipes.funcs.io import read_yaml, dictionary_stripper
from itertools import product, chain

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
                    help="file name, format: .h5mu")
parser.add_argument("--basis_dict",
                    default="X_umap",
                    help="basis in obsm to plot")
parser.add_argument("--categorical_variables",
                    default=None,
                    help="categorical values within mdata.obs")
parser.add_argument("--continuous_variables",
                    default=None,
                    help="continuous values within mdata.obs")
parser.add_argument("--base_figure_dir", default="figures/",
                    help="figures path")
parser.add_argument("--fig_suffix", default="variables.png",
                    help="figures path")
args, opt = parser.parse_known_args()

L.info("running with:")
L.info(args)
# args = argparse.Namespace(infile='../run_clustering_mm/mdata_clustered.h5mu',
#  basis_dict="{'prot': ['X_umap', 'X_pca']}", 
#  base_figure_dir='./')

# ---------
def main(adata, mod, plot_features,  basis, fig_suffix):
    # define file name
    fname_prefix = "_".join(["_" + mod,  fig_suffix])
    # get features
    L.info("plotting")
    L.info(plot_features)
    sc.settings.figdir  = os.path.join(args.base_figure_dir, mod)
    plot_features = [pf for pf in plot_features if pf in mdata.obs.columns]
    pointsize = 120000 / adata.shape[0]
    mu.pl.embedding(adata, basis=basis, color=plot_features, size=pointsize, save = fname_prefix + ".png")


L.debug("load data")
mdata = mu.read(args.infile)


# get bases
basis_dict = dictionary_stripper(read_yaml(args.basis_dict))

if args.categorical_variables is not None:
    cat_vars = dictionary_stripper(read_yaml(args.categorical_variables))
    
    uniq_discrete = list(set(chain(*cat_vars.values())))
    # make sure they are categories
    mdata.obs[uniq_discrete] = mdata.obs[uniq_discrete].apply(lambda x: x. astype('category'))
else:
    cat_vars = {}


if args.continuous_variables is not None:
    try:
        cont_vars = dictionary_stripper(read_yaml(args.continuous_variables))
    except AttributeError:
        # this assumes that we have tried to parse a dict and nstead found a string
        # there is probably a better solution
        cont_vars = args.continuous_variable_dict
else:
    cont_vars = {}




# we have multimodal object
params_dict = {}
for mod in ['rna', 'prot', 'atac', 'multimodal']:
    params_dict['mod' ] = mod
    try:
        basis_list = basis_dict[mod]
    except KeyError:
        L.info("no basis given for mod: %s, no plot's produced" % mod)
        continue
    params_dict['plot_features'] = []
    if mod in cat_vars.keys():
        params_dict['plot_features'] = cat_vars[mod]
    if mod in cont_vars.keys():
        params_dict['plot_features'] = params_dict['plot_features'] +  cont_vars[mod]
    params_dict['fig_suffix'] = args.fig_suffix
    L.info("plotting with params")
    L.info(params_dict)
    if len(basis_list) > 0:
        for basis in basis_list:
            main(adata=mdata, basis = mod + ":" + basis, **params_dict)



L.info("Done")

