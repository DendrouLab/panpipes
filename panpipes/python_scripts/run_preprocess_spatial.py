'''
Preprocess spatial transcriptomics data
'''
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import numpy as np
import scanpy as sc
import muon as mu
import scanpy.experimental as sce

import os
import argparse
import sys
import logging
import re

L = logging.getLogger()
L.setLevel(logging.INFO)
log_handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s: %(levelname)s - %(message)s')
log_handler.setFormatter(formatter)
L.addHandler(log_handler)
L.debug("testing logger works")


from panpipes.funcs.scmethods import X_is_raw


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

parser.add_argument("--norm_hvg_flavour",
                    default=None,
                    help="how to normalize the data and perform HVG selection, either 'squidpy' or 'seurat'")
parser.add_argument("--n_top_genes",
                    default=2000,
                    help="how many variable genes to select")
parser.add_argument("--filter_by_hvg",
                    default=False,
                    help="if True, subset the data to highly-variable genes after finding them")
parser.add_argument("--hvg_batch_key",
                    default=None,
                    help="if specified, highly-variable genes are selected within each batch separately and merged")
parser.add_argument("--squidpy_hvg_flavour",
                    default='seurat',
                    help="flavor for identifying highly variable genes. 'seurat', 'cellranger' or 'seurat_v3'")
parser.add_argument("--min_mean",
                    default=0.05,
                    help="min mean expression for HVG selection")
parser.add_argument("--max_mean",
                    default=1.5,
                    help="max mean expression for HVG selection")
parser.add_argument("--min_disp",
                    default=0.5,
                    help="dispersion min for HVG selection")
parser.add_argument("--theta",
                    default=100,
                    help="the negative binomial overdispersion parameter theta for Pearson residuals")
parser.add_argument("--clip",
                    default=None,
                    help="determines if and how residuals are clipped. Further information: https://scanpy.readthedocs.io/en/stable/generated/scanpy.experimental.pp.normalize_pearson_residuals.html")
parser.add_argument("--n_pcs",
                    default=50,
                    help="how many PCs to compute")


args, opt = parser.parse_known_args()

L.info(args)

figdir = args.figdir

if not os.path.exists(figdir):
    os.mkdir(figdir)

sc.settings.figdir = figdir
sc.set_figure_params(scanpy=True, fontsize=14, dpi=300, facecolor='white', figsize=(5,5))

L.info("Reading in MuData from '%s'" % args.input_mudata)
mdata = mu.read(args.input_mudata)
spatial = mdata.mod['spatial']

input_data = os.path.basename(args.input_mudata)
pattern = r"_filtered.h5(.*)"
match = re.search(pattern, input_data)
sprefix = input_data[:match.start()]


# check if raw data is available
#maybe layer of raw data as parameter
L.info("Checking if raw data is available")
if X_is_raw(spatial):
    L.info("Saving raw counts from .X to .layers['raw_counts']")
    spatial.layers['raw_counts'] = spatial.X.copy()
elif "raw_counts" in spatial.layers :
    L.info(".layers['raw_counts'] already exists and copying it to .X")
    spatial.X = spatial.layers['raw_counts'].copy()
else:
    sys.exit("X is not raw data and 'raw_counts' layer not found")

L.info("Using raw counts for HVG selection & normalization")

# Normalization + HVG selection based on flavour
if args.norm_hvg_flavour == "squidpy":
    if args.squidpy_hvg_flavour == "seurat_v3":
        L.info("Running HVG selection with flavor seurat_v3")
        sc.pp.highly_variable_genes(spatial, flavor="seurat_v3", n_top_genes=int(args.n_top_genes), subset=args.filter_by_hvg,
                                    batch_key=args.hvg_batch_key)
        L.info("Log-normalizing data")
        sc.pp.normalize_total(spatial)
        sc.pp.log1p(spatial)
    else:
        L.info("Log-normalizing data")
        sc.pp.normalize_total(spatial)
        sc.pp.log1p(spatial)
        L.info("Running HVG selection with flavor %s" % args.squidpy_hvg_flavour)
        sc.pp.highly_variable_genes(spatial, flavor=args.squidpy_hvg_flavour,
                                    min_mean=float(args.min_mean),
                                    max_mean=float(args.max_mean),
                                    min_disp=float(args.min_disp), subset=args.filter_by_hvg, batch_key=args.hvg_batch_key)
    L.info("Saving log-normalized counts to .layers['lognorm']")
    spatial.layers["lognorm"] = spatial.X.copy()
    # plot HVGs:
    sc.pl.highly_variable_genes(spatial, show=False, save="_genes_highlyvar" + "."+ sprefix+ ".png")

elif args.norm_hvg_flavour == "seurat":
    if args.clip is None:
        clip = args.clip
    elif args.clip == "None":
        clip = None
    elif args.clip == "np.Inf":
        clip = np.Inf
    else:
        clip = float(args.clip)
    L.info("Running Pearson Residuals HVG selection")
    sce.pp.highly_variable_genes(spatial, theta=float(args.theta), clip=clip, n_top_genes=int(args.n_top_genes),
                                 batch_key=args.hvg_batch_key, flavor='pearson_residuals',
                                 layer="raw_counts", subset=args.filter_by_hvg)
    L.info("Running Pearson Residuals normalization")
    sce.pp.normalize_pearson_residuals(spatial, theta=float(args.theta), clip=clip, layer="raw_counts")
    L.info("Saving log-normalized counts to .layers['norm_pearson_resid']")
    spatial.layers["norm_pearson_resid"] = spatial.X.copy()
else:
    # error or warning?
    L.warning("No normalization and HVG selection was performed! To perform, please specify the 'norm_hvg_flavour' as either 'squidpy' or 'seurat'")


if "highly_variable" in spatial.var:
    L.info("You have %s Highly Variable Features", np.sum(spatial.var.highly_variable))



#PCA
L.info("Running PCA")
sc.pp.pca(spatial, n_comps=int(args.n_pcs), svd_solver='arpack', random_state=0)
L.info("Plotting PCA")
sc.pl.pca(spatial, save = "_vars" + "."+ sprefix+".png")
sc.pl.pca_variance_ratio(spatial, log=True, n_pcs=int(args.n_pcs), save= "."+ sprefix+".png")


        
mdata.update()
L.info("Saving updated MuData to '%s'" % args.output_mudata)
mdata.write(args.output_mudata)

L.info("Done")

