'''
scanpy Normalize script ATAC
'''
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

# import scipy.io
# import matplotlib.pyplot as plt
# import numpy as np
# import pandas as pd
import scanpy as sc
# import anndata as ad
import os
import argparse
# import seaborn as sns
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

from panpipes.funcs.io import read_anndata, write_anndata, write_obs
from panpipes.funcs.processing import check_for_bool
from panpipes.funcs.scmethods import identify_isotype_outliers
from panpipes.funcs.plotting import adjust_x_axis

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
                    help="wheter to binarize the peak counts matrix before any other normalization, like TFIDF")
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

args, opt = parser.parse_known_args()

L.info("running with args:")
L.debug(args)

figdir = args.figdir

if not os.path.exists(figdir):
    os.mkdir(figdir)

sc.settings.figdir = figdir
sc.set_figure_params(scanpy=True, fontsize=14, dpi=300, facecolor='white', figsize=(5,5))

if args.is_paired:
    args.use_muon= True
    L.info("I'm working on a multiome experiment")
else:
    L.info("qc'ing an atac standalone assay")

mdata = mu.read(args.input_mudata)
atac = mdata.mod['atac']


# NORMALIZE

if args.binarize:
    L.info("binarizing peak count matrix")
    ac.pp.binarize(atac)    

if args.normalize is not None:
    if args.normalize == "log1p":
        if args.binarize:
            L.warning("Careful, you have decided to binarize data but also to normalize per cell and log1p. Not sure this is meaningful")
        sc.pp.normalize_per_cell(atac,counts_per_cell_after=1e4)
        sc.pp.log1p(atac)
    elif args.normalize == "TFIDF":
        if args.TFIDF_flavour=="signac":
            ac.pp.tfidf(atac, scale_factor=1e4, log_tfidf=True, log_tf=False, log_idf=False)
        elif args.TFIDF_flavour=="logTF":
            ac.pp.tfidf(atac, scale_factor=1e4, log_tfidf=False, log_tf=True, log_idf=False)    
        elif args.TFIDF_flavour=="logIDF":
            ac.pp.tfidf(atac, scale_factor=1e4, log_tfidf=False, log_tf=False, log_idf=True)
else: 
    L.error("Require a normalization strategy")
    sys.exit("Exiting because no normalization was specified. If None was intended, check your pipeline.yml file")


#highly variable feat selection
sc.pp.highly_variable_genes(atac, min_mean=float(args.min_mean), max_mean=float(args.max_mean), min_disp=float(args.min_disp))

mdata.update()
mdata.write(args.output_mudata)