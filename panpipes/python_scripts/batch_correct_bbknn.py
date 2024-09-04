import pandas as pd
import scanpy as sc
import argparse
import os 
import muon as mu

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


# there is a problem with numvba pushong a warning about TBB which is causing the pipeline to crash
# we might need to need to suppress this warning

# parse arguments
parser = argparse.ArgumentParser()

parser.add_argument('--input_anndata',
                    default='adata_scaled.h5ad',
                    help='')
parser.add_argument('--dimred',
                    default='PCA',
                    help='which dimred to expect, relevant for ATAC')
parser.add_argument('--output_csv', default='batch_correction/umap_bc_bbknn.csv',
                    help='')
parser.add_argument('--integration_col', default='batch',
                    help='')
parser.add_argument('--neighbors_within_batch', default='3',
                    help='How many top neighbours to report for each batch; total number of neighbours will be this number times the number of batches.')
parser.add_argument('--neighbors_n_pcs', default='50',
                    help='How many top neighbours to report for each batch; total number of neighbours will be this number times the number of batches.')
parser.add_argument('--modality', default='',
                    help='which mod to run bbknn on')


args, opt = parser.parse_known_args()


L.info("Running with params: %s", args)


#adata = read_anndata(args.input_anndata, use_muon=use_muon, modality="rna")
L.info("Reading in MuData from '%s'" % args.input_anndata)
mdata = mu.read(args.input_anndata)
adata = mdata.mod[args.modality] 

nnb = int(args.neighbors_within_batch)
# bbknn can't integrate on 2+ variables, so create a fake column with combined information
columns = [x.strip() for x in args.integration_col.split(",")]

if args.dimred == "PCA":
    dimred = "X_pca"
elif args.dimred == "LSI":
    dimred = "X_lsi"

if dimred not in adata.obsm:
    L.warning("Dimred '%s' could not be found in adata.obsm. Computing PCA with default parameters." % dimred)
    dimred = "X_pca" 
    n_pcs = 50
    if adata.var.shape[0] < n_pcs:
        L.info("You have less features than number of PCs you intend to calculate")
        n_pcs = adata.var.shape[0] - 1
        L.info("Setting n PCS to %i" % int(n_pcs)) 
    L.info("Scaling data")   
    sc.pp.scale(adata)
    L.info("Computing PCA")
    sc.tl.pca(adata, n_comps=n_pcs, 
                    svd_solver='arpack', 
                    random_state=0) 

L.info("Preparing for integration")

if len(columns) > 1:
    L.info("Using 2 columns to integrate on more variables.")
    # comb_columns = "_".join(columns)
    adata.obs["comb_columns"] = adata.obs[columns].apply(lambda x: '|'.join(x), axis=1)

    # make sure that batch is a categorical
    adata.obs["comb_columns"] = adata.obs["comb_columns"].astype("category")
    # run bbknn
    L.info("Running BBKNN")
    adata = sc.external.pp.bbknn(adata, batch_key="comb_columns", copy=True, neighbors_within_batch=nnb)
else:
    adata.obs[args.integration_col] = adata.obs[args.integration_col].astype("category")
    # run bbknn
    L.info("Running BBKNN")
    adata = sc.external.pp.bbknn(adata,
                        use_rep=dimred,
                        batch_key=args.integration_col,
                        copy=True,
                        n_pcs = int(args.neighbors_n_pcs),
                        neighbors_within_batch=nnb)  # calculates the neighbours

L.info("Computing UMAP")
sc.tl.umap(adata)

# write out
L.info("Saving UMAP coordinates to csv file '%s" % args.output_csv)
umap = pd.DataFrame(adata.obsm['X_umap'], adata.obs.index)
umap.to_csv(args.output_csv)



# save full bbknn anndata in tmp, cause need more than just neighbors to work 
outfiletmp = ("tmp/bbknn_scaled_adata_" + args.modality + ".h5ad" )

L.info("Saving AnnData to '%s'" % outfiletmp)
write_anndata(adata, outfiletmp, use_muon=False, modality=args.modality)

L.info("Done")

