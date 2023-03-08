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
parser.add_argument("--figdir",
                    default=None,
                    help="figures directory")
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
    pattern=re.compile(r'^leiden|^louvain_res') 
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
            ax.set(xlabel="1", ylabel="2")
        fig.suptitle(ok, y=1.0)
        fig.savefig(os.path.join(figdir, ok +  "_clusters.png"))

L.debug("load data")
mdata = read(args.infile)

# detemin initial figure directory based on object type
if args.figdir is None: 
    if type(mdata) is MuData :
        figdir  = "multimodal/figures"
        L.info("plotting multimodal umaps")
    else:
        figdir  = "figures"
        L.info("plotting umaps from anndata")
else:
    figdir = args.figdir
# do plotting
main(mdata, figdir)


# we also need to plot per modality
if type(mdata) is MuData:
    for mod in mdata.mod.keys():
        L.info("plotting for modality: %s" % mod)
        figdir  = os.path.join(mod, "figures")
        main(mdata[mod], figdir)
