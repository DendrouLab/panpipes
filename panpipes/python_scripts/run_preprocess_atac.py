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
from panpipes.funcs.scmethods import X_is_raw
from panpipes.funcs.scmethods import findTopFeatures_pseudo_signac
from panpipes.funcs.scmethods import lsi, calc_tech_corr, extract_lsi
from panpipes.funcs.plotting import plot_lsi_corr


L = logging.getLogger()
L.setLevel(logging.INFO)
log_handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s: %(levelname)s - %(message)s')
log_handler.setFormatter(formatter)
L.addHandler(log_handler)
L.debug("testing logger works")

sc.settings.verbosity = 3

parser = argparse.ArgumentParser()

parser.add_argument("--input_mudata",
                    default="mudata_unfilt.h5mu",
                    help="")
parser.add_argument("--output_mudata",
                    default="mudata_unfilt.h5mu",
                    help="")
parser.add_argument("--figdir",
                    default="./figures/atac",
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
parser.add_argument('--n_top_features', default=None,
                    help = "if using scanpy and setting this, selecting these many HVF")
parser.add_argument("--dimred",
                    default="PCA",
                    help="which dimensionality red to use for atac")
parser.add_argument("--n_comps",
                    default=50,
                    help="how many components to compute")
parser.add_argument("--solver",
                    default="arpack",
                    help="what pca solver to use")
parser.add_argument("--dim_remove",
                    default=None,
                    help="which dimensionality red components to remove")
parser.add_argument("--feature_selection_flavour",
                    default=None,
                    help="which HVF selection to perform, 'scanpy' or 'signac'")
parser.add_argument("--min_cutoff",
                    default=None,
                    help="cutoff for Signac's HVF selection")
parser.add_argument("--filter_by_hvf", default=False) 
parser.add_argument("--color_by", default="batch") 


args, opt = parser.parse_known_args()

L.info("running with args:")
L.info(args)

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
if X_is_raw(atac):
    if args.binarize is True:
        L.info("binarizing peak count matrix")
        ac.pp.binarize(atac)    

    if args.normalize is not None:
        if args.normalize == "log1p":
            if args.binarize:
                L.warning("Careful, you have decided to binarize data but also to normalize per cell and log1p. Not sure this is meaningful")
            sc.pp.normalize_per_cell(atac,counts_per_cell_after=1e4)
            sc.pp.log1p(atac)
            atac.layers["logged_counts"] = atac.X.copy()
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
        L.error("We Require a normalization strategy for ATAC")
        sys.exit("Exiting because no normalization was specified. If None was intended, check your pipeline.yml file")


#highly variable feature selection
if "highly_variable" in atac.var: 
    L.warning( "You have %s Highly Variable Features", np.sum(atac.var.highly_variable))
else:

    if args.feature_selection_flavour == "scanpy":
        if args.n_top_features is None:
            sc.pp.highly_variable_genes(atac, min_mean=float(args.min_mean), 
            max_mean=float(args.max_mean), min_disp=float(args.min_disp))
        else:
            sc.pp.highly_variable_genes(atac, n_top_genes=int(args.n_top_features))
    elif args.feature_selection_flavour == "signac":
        findTopFeatures_pseudo_signac(atac, args.min_cutoff)
    else:
        L.warning("No highly variable feature selection was performed!")

if "highly_variable" in atac.var: 
    L.warning( "You have %s Highly Variable Features", np.sum(atac.var.highly_variable))

# filter by hvgs
filter_by_hvgs = args.filter_by_hvf

if filter_by_hvgs:
    L.info("filtering object to only include highly variable features")
    mdata.mod["atac"] = atac[:, atac.var.highly_variable]
    if isinstance(mdata, mu.MuData):
        mdata.update()
    atac = mdata["atac"]
    genes = atac.var
    genes['gene_name'] = atac.var.index
    genes.to_csv("filtered_variable_features.tsv", sep="\t", index=True)

# The combined steps of TF-IDF followed by SVD are known as latent semantic indexing (LSI), 
# and were first introduced for the analysis of scATAC-seq data by Cusanovich et al. 2015.

if args.dimred == "LSI" and args.normalize == "TFIDF":
    lsi(adata=atac, num_components=int(args.n_comps))
else:
    L.info("Applying LSI on logged_counts is not recommended. Changing dimred to PCA")
    args.dimred = "PCA"
if args.dimred == "PCA":
    sc.pp.scale(atac)
    atac.layers["scaled_counts"] = atac.X.copy()
    sc.tl.pca(atac, n_comps=int(args.n_comps), svd_solver=args.solver, 
    random_state=0)

col_variables = args.color_by.split(",")
col_variables = [a.strip() for a in col_variables]
col_use = [var for var in col_variables if var in atac.obs.columns]

#some plotting before removing the components

if args.dimred =='PCA':
    sc.pl.pca_variance_ratio(atac, log=True, n_pcs=int(args.n_comps), save=".png")
    sc.pl.pca(atac, color=col_use, save = "_vars.png")
    sc.pl.pca_loadings(atac, components="1,2,3,4,5,6", save = ".png")
    sc.pl.pca_overview(atac, save = ".png")
if args.dimred =='LSI':
    sc.pl.embedding(atac, color=col_use,basis="X_lsi", save = ".png")
    correlation_df = calc_tech_corr(extract_lsi(atac), 
        tech_covariates = ['n_genes_by_counts', 'total_counts'])

    correlation_df.to_csv(os.path.join(args.figdir,"tech_covariates_corr_LSI.tsv"), sep='\t')

    plot_lsi_corr(correlation_df, tech_covariates=['n_genes_by_counts', 'total_counts'],
        filename=os.path.join(args.figdir,"LSI_corr_plot.png"))


if args.dim_remove is not None:
    L.info("removing component from dimred")
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

