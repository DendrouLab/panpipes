
import scanpy as sc
import os
import argparse
import muon as mu
from muon import atac as ac
from panpipes.funcs.processing import check_for_bool


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

parser.add_argument("--input_mudata",
                    default="adata_scaled.h5ad",
                    help="")
parser.add_argument('--modality',
                    default=False,
                    help='')
parser.add_argument("--fig_dir", default="figures/")
# pca options
parser.add_argument("--n_pcs", default=50)
parser.add_argument("--color_by", default="batch")

args, opt = parser.parse_known_args()

sc.settings.autoshow=False
sc.settings.figdir=args.fig_dir 
sc.settings.set_figure_params(dpi_save=300)


mdata = mu.read(args.input_mudata)


if args.color_by is not None:
    qc_vars=args.color_by.split(",")
    qc_vars= [a.strip() for a in qc_vars]
else:
    qc_vars = ["sample_id"]

my_dict = {"rna": 'gex', "prot": 'adt', "atac": 'atac', "rep": 'rep'}
for mod in mdata.mod.keys():
    sc.settings.figdir=os.path.join(args.fig_dir,my_dict[mod])
    pointsize = 120000 / adata.shape[0]
    if "X_pca" in mdata[mod].obsm.keys():
        L.info("Found a PCA in %s, Plotting it" % mod)
        sc.pl.pca(mdata[mod], color=qc_vars, size=pointsize)
        
    if "X_lsi" in mdata[mod].obsm.keys(): 
        L.info("Found a LSI in %s, Plotting it" % mod)
        ac.pl.lsi(mdata[mod],color=qc_vars)
        #https://satijalab.org/signac/reference/DepthCor.html
        #muon.pl.scatter(mdata[mod], x=)

L.info("Done")

