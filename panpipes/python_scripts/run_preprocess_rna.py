## Processing mudata based on n_genes, percent mito and min cells/
## originally written by Tom Thomas (https://github.com/tomthomas3000/TAURUS)
## adapted for this pipeline by Charlotte Rich-Griffin 2020-09-30


# import numpy as np
import pandas as pd
import numpy as np
import scanpy as sc
# import scipy.io
# import matplotlib
# matplotlib.use("Agg")
# import matplotlib.pyplot as plt
# plt.ioff()
import os
import argparse
import re
import muon as mu


from panpipes.funcs.processing import check_for_bool
from panpipes.funcs.scmethods import X_is_raw

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
                    default='data/mudata-n5000_filt.h5ad',
                    help='')
parser.add_argument('--output_scaled_mudata',
                    default=None,
                    help='if not specified then the output will be written to output_logged_mudata path')
parser.add_argument('--use_muon',
                    default=False, type=check_for_bool,
                    help='')  
parser.add_argument('--fig_dir', 
                    default="./figures",
                    help="save plots here")
parser.add_argument('--exclude_file', default=None,
                    help='')
parser.add_argument('--exclude', default="exclude",
                    help='how are genes that need to be excl. from HVG are labeled')

# highly variable genes options
parser.add_argument('--flavor', default='seurat')
parser.add_argument('--n_top_genes', default=None)
parser.add_argument('--min_mean', default=0.0125)
parser.add_argument('--max_mean', default=3)
parser.add_argument('--min_disp', default=0.5)
parser.add_argument("--filter_by_hvg", default=False, type=check_for_bool)
parser.add_argument('--hvg_batch_key', default=None)
# regress out options
parser.add_argument('--regress_out', default=None)
# scale options
parser.add_argument('--scale', default=True, type=check_for_bool)
parser.add_argument('--scale_max_value', default=None)
# pca options
parser.add_argument("--n_pcs", default=50)
parser.add_argument("--pca_solver", default="arpack")
parser.add_argument("--color_by", default="batch") 

parser.set_defaults(verbose=True)
args, opt = parser.parse_known_args()
L.info(args)

figdir = args.fig_dir
if not os.path.exists(figdir):
    os.mkdir(figdir)

sc.settings.figdir = figdir
sc.set_figure_params(scanpy=True, fontsize=14, dpi=300, facecolor='white', figsize=(5,5))

mdata = mu.read(args.input_mudata)

# if we have actually loaded an anndata object, make into a mudata
if isinstance(mdata, sc.AnnData):
    mdata =  mu.MuData({'rna': mdata})

adata = mdata['rna']
# resolve multi-column batch for hvg batch key
if args.hvg_batch_key is not None:
    columns = [x.strip() for x in args.hvg_batch_key.split(",")]
    if len(columns) > 1: 
        L.info("combining batch columns into one column 'hvg_batch_key'")
        adata.obs["hvg_batch_key"] = adata.obs[columns].apply(lambda x: '|'.join(x), axis=1)
        # make sure that batch is a categorical
        adata.obs["hvg_batch_key"] = adata.obs["hvg_batch_key"].astype("category")
        hvg_batch_key="hvg_batch_key"
    else:
        hvg_batch_key=columns[0]
else:
    hvg_batch_key=None

# save raw counts as a layer
if X_is_raw(adata):
    adata.layers['raw_counts'] = adata.X.copy()
elif "raw_counts" in adata.layers :
    L.info("raw_counts layer already exists")
    adata.X = adata.layers['raw_counts'].copy()
else:
    L.error("X is not raw data and raw_counts layer not found")
    sys.exit("X is not raw data and raw_counts layer not found")


# sc.pp.highly variable genes Expects logarithmized data, 
# except when flavor='seurat_v3' in which count data is expected.
# change the order accordingly
if X_is_raw(adata):
    L.info("normalise, log and calculate highly variable genes")
    if args.flavor == "seurat_v3":
        if args.n_top_genes is None:
            raise ValueError("if seurat_v3 is used you must give a n_top_genes value")
        else:
            sc.pp.highly_variable_genes(adata, flavor="seurat_v3", 
                                        n_top_genes=int(args.n_top_genes),
                                        batch_key=hvg_batch_key)
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    else:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        L.debug(adata.uns['log1p'])
        sc.pp.highly_variable_genes(adata, flavor=args.flavor,
                                    min_mean=float(args.min_mean), 
                                    max_mean=float(args.max_mean),
                                    min_disp=float(args.min_disp),
                                    batch_key=hvg_batch_key)
        L.debug(adata.uns['log1p'])
if "highly_variable" in adata.var: 
    L.warning( "You have %s Highly Variable Features", np.sum(adata.var.highly_variable))
    sc.pl.highly_variable_genes(adata,show=False, save ="_genes_highlyvar.png")
else:
    sys.exit("I cannot find Highly Variable Features in your RNA")

if isinstance(mdata, mu.MuData):
    mdata.update()

# exclude hvgs if there is an exclude file
if args.exclude_file is not None:
    if os.path.exists(args.exclude_file):
        L.info("exclude genes from hvg")
        customgenes = pd.read_csv(args.exclude_file) 
        custom_cat = list(set(customgenes['group'].tolist()))
        cat_dic = {}
        for cc in custom_cat:
            cat_dic[cc] = customgenes.loc[customgenes["group"] == cc,"feature"].tolist()
        exclude_action = str(args.exclude)
        excl = cat_dic[exclude_action]
        L.info(len(excl))
        L.info("number of hvgs prior to filtering")
        L.info(adata.var.highly_variable.sum())
        adata.var['hvg_exclude'] = [True if gg in excl else False for gg in adata.var.index]
        L.debug("num to exclude %s" %  (adata.var['hvg_exclude'] & adata.var['highly_variable']).sum() )
        L.debug(adata)
        if any(adata.var['hvg_exclude'] & adata.var['highly_variable']):
            adata.var['highly_variable'].loc[adata.var.hvg_exclude] = False
            adata.var.loc[adata.var.hvg_exclude, 'highly_variable_rank'] = np.nan
            
            L.info("number of hvgs after filtering")
            L.info(adata.var.highly_variable.sum())
            sc.pl.highly_variable_genes(adata,show=False, save ="_exclude_genes_highlyvar.png")
    else:
        sys.exit("exclusion file %s not found, check the path and try again" % args.exclude_file)

if isinstance(mdata, mu.MuData):
    mdata.update()

L.debug(adata.uns['log1p'])
# filter by hvgs
filter_by_hvgs = args.filter_by_hvg

if filter_by_hvgs is True:
    L.info("filtering object to only include highly variable genes")
    mdata.mod["rna"] = adata[:, adata.var.highly_variable]
    if isinstance(mdata, mu.MuData):
        mdata.update()
    adata = mdata.mod["rna"]
    genes = adata.var
    genes['gene_name'] = adata.var.index
    genes.to_csv("filtered_genes.tsv", sep="\t", index=True)
 

adata.layers['logged_counts'] = adata.X.copy()
L.debug(adata.uns['log1p'])
# regress out
if args.regress_out is not None:
    regress_opts = args.regress_out.split(",")
    L.info("regressing out %s" % regress_opts)
    sc.pp.regress_out(adata, regress_opts)


if args.scale is True and args.scale_max_value is None:
    L.info("scaling data with default parameters")
    sc.pp.scale(adata)
elif args.scale is True:
    L.info("scaling data to max value %i" % int(args.scale_max_value))
    sc.pp.scale(adata, max_value=int(args.scale_max_value))
else:
    L.info("not scaling data")
    pass

L.debug(adata.uns['log1p'])
# run pca
L.info("running pca")

if adata.var.shape[0] < int(args.n_pcs):
    L.info("You have less features than number of PCs you intend to calculate")
    n_pcs = adata.var.shape[0] - 1
    L.info("Setting n PCS to %i" % int(n_pcs))    
else:
    n_pcs = int(args.n_pcs)
sc.tl.pca(adata, n_comps=n_pcs, 
                svd_solver=args.pca_solver, 
                random_state=0) 

# do some plots!
sc.pl.pca_variance_ratio(adata, log=True, n_pcs=n_pcs, save=".png")

col_variables = args.color_by.split(",")
col_variables = [a.strip() for a in col_variables]

col_use = [var for var in col_variables if var in adata.obs.columns]
sc.pl.pca(adata, color=col_use, save = "_vars.png")
sc.pl.pca_loadings(adata, components="1,2,3,4,5,6", save = ".png")
sc.pl.pca_overview(adata, save = ".png")


mdata.update()
# save the scaled adata object

mdata.write(args.output_scaled_mudata)
L.debug(adata.uns['log1p'])

L.info("done")
