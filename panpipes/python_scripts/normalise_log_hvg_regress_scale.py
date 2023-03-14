## Processing anndata based on n_genes, percent mito and min cells/
## originally written by Tom Thomas (https://github.com/tomthomas3000/TAURUS)
## adapted for this pipeline by Charlotte Rich-Griffin 2020-09-30


# import numpy as np
import pandas as pd
import scanpy as sc
# import scipy.io
# import matplotlib.pyplot as plt
import os
# import anndata
import argparse
from muon import MuData
from anndata import AnnData
import re


from panpipes.funcs.processing import check_for_bool
from panpipes.funcs.io import read_anndata, write_anndata

import sys
import logging
L = logging.getLogger()
L.setLevel(logging.INFO)
log_handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s: %(levelname)s - %(message)s')
log_handler.setFormatter(formatter)
L.addHandler(log_handler)

# parse arguments
parser = argparse.ArgumentParser()


parser.add_argument('--input_mudata',
                    default='data/anndata_filt.h5mu',
                    help='')
parser.add_argument('--output_logged_mudata',
                    default=None,
                    help='if not specified, then the input file will be overwritten')
parser.add_argument('--output_scaled_anndata',
                    default=None,
                    help='if not specified then the output will be written to output_logged_anndata path')
parser.add_argument('--use_muon',
                    default=False, type=check_for_bool,
                    help='')  
parser.add_argument('--fig_dir', 
                    default="./figures",
                    help="save plots here")
parser.add_argument('--exclude_file', default=None,
                    help='')

# highly variable genes options
parser.add_argument('--flavor', default='seurat')
parser.add_argument('--n_top_genes', default=None)
parser.add_argument('--min_mean', default=0.0125)
parser.add_argument('--max_mean', default=3)
parser.add_argument('--min_disp', default=0.5)
parser.add_argument("--filter_by_hvg", default=False, type=check_for_bool)
# regress out options
parser.add_argument('--regress_out', default=None)
# scale options
parser.add_argument('--scale_data', default=True, type=check_for_bool)
parser.add_argument('--scale_max_value', default=None)
# pca options
parser.add_argument("--n_pcs", default=50)
parser.add_argument("--color_by", default="batch") 

parser.set_defaults(verbose=True)
args, opt = parser.parse_known_args()
L.info(args)

use_muon = args.use_muon

# sc.settings.verbosity = 3
# sc.logging.print_header()

figdir = args.fig_dir
if not os.path.exists(figdir):
    os.mkdir(figdir)

sc.settings.figdir = figdir
sc.set_figure_params(scanpy=True, fontsize=14, dpi=300, facecolor='white', figsize=(5,5))

mdata = mu.read(args.input_anndata)

# if we have actually loaded an anndata object, make into a mudata
if isinstance(mdata, AnnData):
    mdata =  MuData({'rna': mdata})

adata = mdata['rna']


adata.layers['raw_counts'] = adata.X.copy()
# Normalise to depth 10k, store raw data, assess and drop highly variable genes, regress mitochondria and count

# sc.pp.highly variabel genes Expects logarithmized data, except when flavor='seurat_v3' in which count data is expected.
# change the order accordingly
L.info("run hvgs")
if args.flavor == "seurat_v3":
    if args.n_top_genes is None:
        raise ValueError("if seurat_v3 is used you must give a n_top_genes value")
        # sc.pp.highly_variable_genes(adata, flavor="seurat_v3",)
    else:
        sc.pp.highly_variable_genes(adata, flavor="seurat_v3",
                                    n_top_genes=int(args.n_top_genes))
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
else:
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, flavor=args.flavor,
                                min_mean=float(args.min_mean), max_mean=float(args.max_mean),
                                min_disp=float(args.min_disp))

sc.pl.highly_variable_genes(adata,show=False, save ="_genes_highlyvar.png")


# exclude hvgs if there is an exclude file
if args.exclude_file is not None:
    L.info("exclude genes from hvg")
    excl = pd.read_csv(args.exclude_file, delimiter='\t')
    L.info(excl.status.value_counts())
    L.info("number of hvgs prior to filtering")
    L.info(adata.var.highly_variable.sum())
    #adata.var['hvg_exclude'] = [True if gg in excl.gene_id.tolist() else False for gg in adata.var.index]
    adata.var['hvg_exclude'] = [True if gg in excl.gene_name.tolist() else False for gg in adata.var.index]
    adata.var.loc[adata.var.hvg_exclude, 'highly_variable']= False
    L.info("number of hvgs after filtering")
    L.info(adata.var.highly_variable.sum())
    sc.pl.highly_variable_genes(adata,show=False, save ="_exclude_genes_highlyvar.png")

# filter by hvgs
filter_by_hvgs = args.filter_by_hvg

if filter_by_hvgs is True:
    adata = adata[:, adata.var.highly_variable]
    if isinstance(adata, MuData):
        adata.update()


genes = adata.var
genes['gene_name'] = adata.var.index
genes.to_csv("filtered_genes.tsv", sep="\t", index=True)
# save object (has to occur before scaling) 
if args.output_logged_anndata is None:
    args.output_logged_anndata = args.input_anndata

if args.output_scaled_anndata == args.output_logged_anndata:
    pass
else:
    write_anndata(mdata, args.output_logged_anndata, use_muon=use_muon, modality="all")


#------
# OPINION: I don't like using raw very much. I'd kinda prefer to use a layer
# adata.raw = adata # this means that log normalised counts are saved in raw
adata.layers['logged_counts'] = adata.X.copy()
# regress out
if args.regress_out is not None:
    regress_opts = args.regress_out.split(",")
    sc.pp.regress_out(adata, regress_opts)

if args.scale_data is True:
    if args.scale_max_value is not None:
        sc.pp.scale(adata)
else:
    sc.pp.scale(adata, max_value=args.scale_max_value)

# run pca
L.info("running pca")
sc.tl.pca(adata, n_comps=250, svd_solver='arpack', random_state=0) #given args above this should work
# extract pca coordinates for plotting (in R??)
pca_coords = pd.DataFrame(adata.obsm['X_pca'])
# add in the rownames 
pca_coords.index = adata.obs_names
# save coordinates to file
# (note this saves values values up to 6 significant figures, because why save 20 for a plot
pca_fname = re.sub(".h5ad|.h5mu", "_pca.txt.gz", args.output_scaled_anndata)
pca_coords.to_csv(pca_fname, sep='\t')



# do some plots!
sc.pl.pca_variance_ratio(adata, log=True, n_pcs=int(args.n_pcs), save=".png")

col_variables = args.color_by.split(",")
# for cv in col_variables:
#     sc.pl.pca(adata, color=cv, save="_" + cv + ".png")

sc.pl.pca(adata, color=col_variables, save = "_vars.png")

sc.pl.pca_loadings(adata, components="1,2,3,4,5,6", save = ".png")
sc.pl.pca_overview(adata, save = ".png")


mdata.update()
# save the scaled adata object

if use_muon:
    # for some reason the X_pca wasn't getting stored in mdata so update the object
    mdata.mod['rna'] = adata
    write_anndata(mdata, args.output_scaled_anndata, use_muon=use_muon, modality="all")
else:
    write_anndata(adata, args.output_scaled_anndata, use_muon=use_muon, modality="rna")
    

L.info("Done")

