"""
Scanpy inbuilt plots for custom markers
CRG 2020-06-19
"""

import scanpy as sc
from muon import read, MuData
sc.settings.autoshow = False
import pandas as pd
import argparse
from anndata import AnnData
import matplotlib
matplotlib.use('agg')
import os 
import re

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
                    default="mdata.h5mu",
                    help="file name, format: .h5mu")
parser.add_argument("--modalities",
                    default="rna",
                    help="list of modalities to search for UMAPs in")
args, opt = parser.parse_known_args()

L.info("running with:")
L.info(args)
# args = argparse.Namespace(infile='mdata_clustered.h5mu', figdir=None)
# ---------

def main(adata,figdir):
    # get all possible umap coords
    pattern="X_umap(.*)"
    obsm_keys = [x for x in adata.obsm.keys() if re.search(pattern, x)]
    L.info("Umap keys founds %s" % obsm_keys)
    # get all possible clustersclusters
    pattern=re.compile(r'^leiden|^louvain') 
    cluster_keys  = [x for x in adata.obs.columns if re.search(pattern, x)]
    L.info("Cluster keys founds %s" % cluster_keys)
    if len(obsm_keys) == 0  or len(cluster_keys) == 0:
        return
    
    # make sure clusters are categories
    for ck in cluster_keys:
        adata.obs[ck] = adata.obs[ck].astype('category')
    # plot all the umaps
    for ok in obsm_keys:
        fig = sc.pl.embedding(adata, basis = ok,color=cluster_keys,
         show=False, return_fig=True, legend_loc='on data')
        for ax in fig.axes:
            ax.set(xlabel="UMAP_1", ylabel="UMAP_2")
        fig.suptitle(ok, y=1.0)
        fig.savefig(os.path.join(figdir, ok +  "_clusters.png"))

def plot_spatial(adata,figdir):
    # get all possible umap coords
    pattern="spatial(.*)"
    obsm_keys = [x for x in adata.obsm.keys() if re.search(pattern, x)]
    L.info("Umap keys founds %s" % obsm_keys)
    # get all possible clustersclusters
    pattern=re.compile(r'^leiden|^louvain') 
    cluster_keys  = [x for x in adata.obs.columns if re.search(pattern, x)]
    L.info("Cluster keys founds %s" % cluster_keys)
    if len(obsm_keys) == 0  or len(cluster_keys) == 0:
        return
    
    # make sure clusters are categories
    for ck in cluster_keys:
        adata.obs[ck] = adata.obs[ck].astype('category')
    # plot all the umaps
    for ok in obsm_keys:
        fig = sc.pl.embedding(adata, basis = ok,color=cluster_keys,
         show=False, return_fig=True, legend_loc='on data')
        for ax in fig.axes:
            ax.set(xlabel="spatial1", ylabel="spatial2")
        fig.suptitle(ok, y=1.0)
        fig.savefig(os.path.join(figdir, ok +  "_clusters.png"))


L.debug("load data")
mdata = read(args.infile)

mods = args.modalities.split(',')
# detemin initial figure directory based on object type

# do plotting
if 'multimodal' in mods:
    if os.path.exists("multimodal/figures") is False:
        os.makedirs("multimodal/figures")
    main(mdata, figdir="multimodal/figures")


# we also need to plot per modality
if type(mdata) is MuData:
    for mod in mdata.mod.keys():
        if mod in mods:
            L.info("plotting for modality: %s" % mod)
            figdir  = os.path.join(mod, "figures")
            if os.path.exists(figdir) is False:
                os.makedirs(figdir)
            if mod == "spatial": # added separate function for spatial
                plot_spatial(mdata[mod], figdir)
            else:
                main(mdata[mod], figdir)





L.info('done')