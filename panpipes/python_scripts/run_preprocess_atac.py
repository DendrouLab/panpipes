'''
scanpy Normalize script ATAC
'''
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import os
import argparse
import seaborn as sns
from muon import atac as ac
import muon as mu
import sys
import logging
L = logging.getLogger()
L.setLevel(logging.INFO)
log_handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s: %(levelname)s - %(message)s')
log_handler.setFormatter(formatter)
L.addHandler(log_handler)
L.debug("testing logger works")


from panpipes.funcs.scmethods import X_is_raw
from panpipes.funcs.scmethods import findTopFeatures_pseudo_signac


sc.settings.verbosity = 3

parser = argparse.ArgumentParser()

parser.add_argument("--input_mudata",
                    default="mudata_unfilt.h5mu",
                    help="")
parser.add_argument("--output_mudata",
                    default="mudata_unfilt.h5mu",
                    help="")
parser.add_argument("--figdir",
                    default="./figures/",
                    help="path to save the figures to")
parser.add_argument("--binarize",
                    default=False,
                    help="whether to binarize the peak counts matrix before any other normalization, like TFIDF")
parser.add_argument("--normalize",
                    default="TFIDF",
                    help="what operation to use to normalize data. Options are: log1p, TFIDF, None")
parser.add_argument("--TFIDF_flavour",
                    default="signac",
                    help="how to log the data. Default to signac(ArchR) default (1) from https://satijalab.org/signac/reference/runtfidf")
parser.add_argument("--min_mean",
                    default=0.05,
                    help="min mean expression for HVF selection")
parser.add_argument("--max_mean",
                    default=1.5,
                    help="max mean expression for HVF selection")
parser.add_argument("--min_disp",
                    default=0.5,
                    help="dispersion min for HVF selection")
parser.add_argument("--dimred",
                    default="PCA",
                    help="which dimensionality red to use for atac")
parser.add_argument("--dim_remove",
                    default=None,
                    help="which dimensionality red components to remove")
parser.add_argument("--feature_selection_flavour",
                    default=None,
                    help="which HVF selection to perform, 'scanpy' or 'signac'")
parser.add_argument("--min_cutoff",
                    default=None,
                    help="cutoff for Signac's HVF selection")


args, opt = parser.parse_known_args()

L.info("running with args:")
L.debug(args)

figdir = args.figdir

if not os.path.exists(figdir):
    os.mkdir(figdir)

sc.settings.figdir = figdir
sc.set_figure_params(scanpy=True, fontsize=14, dpi=300, facecolor='white', figsize=(5,5))

mdata = mu.read(args.input_mudata)
atac = mdata.mod['atac']

# save raw counts
if X_is_raw(atac):
    atac.layers["raw_counts"] = atac.X.copy()
elif "raw_counts" in atac.layers :
    L.info("raw_counts layer already exists")
else:
    L.error("X is not raw data and raw_counts layer not found")
    sys.exit("X is not raw data and raw_counts layer not found")


# NORMALIZE

if args.binarize is True:
    L.info("binarizing peak count matrix")
    ac.pp.binarize(atac)    

if args.normalize is not None:
    if args.normalize == "log1p":
        if args.binarize:
            L.warning("Careful, you have decided to binarize data but also to normalize per cell and log1p. Not sure this is meaningful")
        sc.pp.normalize_per_cell(atac,counts_per_cell_after=1e4)
        sc.pp.log1p(atac)
        atac.layers["lognorm"] = atac.X.copy()
    elif args.normalize == "TFIDF":
        if args.TFIDF_flavour=="signac":
            ac.pp.tfidf(atac, scale_factor=1e4, log_tfidf=True, log_tf=False, log_idf=False)
            atac.layers["signac_norm"] = atac.X.copy()
        elif args.TFIDF_flavour=="logTF":
            ac.pp.tfidf(atac, scale_factor=1e4, log_tfidf=False, log_tf=True, log_idf=False)    
            atac.layers["logTF_norm"] = atac.X.copy()
        elif args.TFIDF_flavour=="logIDF":
            ac.pp.tfidf(atac, scale_factor=1e4, log_tfidf=False, log_tf=False, log_idf=True)
            atac.layers["logIDF_norm"] = atac.X.copy()
else: 
    L.error("Require a normalization strategy")
    sys.exit("Exiting because no normalization was specified. If None was intended, check your pipeline.yml file")


#highly variable feature selection

if args.feature_selection_flavour == "scanpy":
	sc.pp.highly_variable_genes(atac, min_mean=float(args.min_mean), max_mean=float(args.max_mean), min_disp=float(args.min_disp))
elif args.feature_selection_flavour == "signac":
	findTopFeatures_pseudo_signac(atac, args.min_cutoff)
else:
    L.warning("No highly variable feature selection was performed!")

if "highly_variable" in atac.var: 
    L.warning( "You have %s Highly Variable Features", np.sum(atac.var.highly_variable))




# The combined steps of TF-IDF followed by SVD are known as latent semantic indexing (LSI), 
# and were first introduced for the analysis of scATAC-seq data by Cusanovich et al. 2015.

if args.dimred == "PCA":
    sc.pp.scale(atac)
    atac.layers["scaled_counts"] = atac.X.copy()
    sc.tl.pca(atac, n_comps=min(50,atac.var.shape[0]-1), svd_solver='arpack', random_state=0) 
if args.dimred == "LSI":
    ac.tl.lsi(atac)

if args.dim_remove is not None:
    dimrem=int(args.dim_remove)
    if args.dimred == "LSI":
        atac.obsm['X_lsi'] = atac.obsm['X_lsi'][:,dimrem:]
        atac.varm["LSI"] = atac.varm["LSI"][:,dimrem:]
        atac.uns["lsi"]["stdev"] = atac.uns["lsi"]["stdev"][dimrem:]
    if args.dimred == "PCA":
        atac.obsm['X_pca'] = atac.obsm['X_pca'][:,dimrem:]
        atac.varm["PCs"] = atac.varm["PCs"][:,dimrem:]
        
mdata.update()
mdata.write(args.output_mudata)

L.info("Done")

