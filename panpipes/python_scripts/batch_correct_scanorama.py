##
# Scanorama is designed to be used in scRNA-seq pipelines downstream of noise-reduction methods,
# including those for imputation and highly-variable gene filtering.
##

import pandas as pd
import scanpy as sc
import os
# import scanorama
import scanpy.external as sce
import argparse
import muon as mu

from panpipes.funcs.processing import check_for_bool
from panpipes.funcs.io import read_anndata, write_anndata
from panpipes.funcs.scmethods import run_neighbors_method_choice

import sys
import logging
L = logging.getLogger()
L.setLevel(logging.INFO)
log_handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s: %(levelname)s - %(message)s')
log_handler.setFormatter(formatter)
L.addHandler(log_handler)

# check for threads number
import multiprocessing 
threads_available = multiprocessing.cpu_count()

# parse arguments
parser = argparse.ArgumentParser()

parser.add_argument('--input_anndata',
                    default='adata_scaled.h5ad',
                    help='')
parser.add_argument('--output_csv', default='batch_correction/umap_bc_scanorama.csv',
                    help='')
parser.add_argument('--integration_col', default='batch',
                    help='')
parser.add_argument('--n_threads', default=1,
                    help="num threads to use for neighbor computations")
parser.add_argument('--neighbors_n_pcs',
                    help="n_pcs")
parser.add_argument('--neighbors_method',
                    help="neighbours method, scanpy or hnsw")
parser.add_argument('--neighbors_k',
                    help="neighbors k")
parser.add_argument('--neighbors_metric',
                    help="neighbor metric, e.g. euclidean or cosine")
parser.add_argument('--batch_size',
                    help="batch size scanorama parameter, default 5000", default="5000")
parser.add_argument('--modality',
                    help="modality")



args, opt = parser.parse_known_args()



L.info("Running with options: %s", args)

# Scanorama is designed to be used in scRNA-seq pipelines downstream of noise-reduction methods,
# including those for imputation and highly-variable gene filtering.
L.info("reading data and starting integration pipeline with script: ")
L.info(os.path.basename(__file__))

#adata = read_anndata(args.input_anndata, use_muon=use_muon, modality="rna")
mdata = mu.read(args.input_anndata)
adata = mdata.mod[args.modality] 
bcs = adata.obs_names.tolist()

# scanorama can't integrate on 2+ variables, so create a fake column with combined information
columns = [x.strip() for x in args.integration_col.split(",")]
if len(columns) > 1:
    L.info("using 2 columns to integrate on more variables")
    # comb_columns = "_".join(columns)
    adata.obs["comb_columns"] = adata.obs[columns].apply(lambda x: '|'.join(x), axis=1)

    # make sure that batch is a categorical
    adata.obs["comb_columns"] = adata.obs["comb_columns"].astype("category")
    adata.obs['batch'] = adata.obs["comb_columns"]
else:
    if args.integration_col != "batch":
        adata.obs['batch'] = adata.obs[args.integration_col]

batches = adata.obs['batch'].unique()
L.info("define scanorama inputs based on batches and number of PCs")
# need contiguous batches
scanorama_input = adata[adata.obs.sort_values(by="batch").index.tolist(), :]

# filter by HVGs to make it equivalent to the old scripts,
# which inputted the scaled object after filtering by hvgs.
scanorama_input = scanorama_input[:, scanorama_input.var.highly_variable]
#also filter the X_PCA to be the number of PCs we actually want to use
scanorama_input.obsm['X_pca'] = scanorama_input.obsm['X_pca'][:,0:int(args.neighbors_n_pcs)]


# run scanoramam using the scanpy integrated approach
L.info("run_scanorama")

sce.pp.scanorama_integrate(scanorama_input, key='batch', batch_size=int(args.batch_size))
L.info("scanorama done")

# not integrated

# old method (simplified)
# alldata = {}
# for bb in batches:
#     alldata[bb] = scanorama_input[scanorama_input.obs['batch'] == bb]
#
# adatas = list(alldata.values())
# X_scanorama = scanorama.integrate_scanpy(adatas)
# scanorama_input.obsm['X_scanorama'] = np.concatenate(X_scanorama)

# put into the original order
scanorama_input = scanorama_input[bcs, :]
# check it worked
if all(scanorama_input.obs_names == bcs):
    L.info("barcode order is correct")
else:
    L.debug("barcode order in  batch corrected object is incorrect")
    sys.exit("barcode order in  batch corrected object is incorrect")

L.info(adata)
# add to the AnnData object
adata.obsm["X_scanorama"] = scanorama_input.obsm["X_scanorama"]

L.info("integration run now calculate umap")

if int(args.neighbors_n_pcs) > adata.obsm['X_scanorama'].shape[1]:
    L.warn(f"N PCs is larger than X_scanorama dimensions, reducing n PCs to  {adata.obsm['X_scanorama'].shape[1] -1}")
n_pcs= min(int(args.neighbors_n_pcs), adata.obsm['X_scanorama'].shape[1]-1)

# run neighbours and umap 
run_neighbors_method_choice(adata, 
    method=args.neighbors_method, 
    n_neighbors=int(args.neighbors_k), 
    n_pcs=n_pcs, 
    metric=args.neighbors_metric, 
    use_rep='X_scanorama',
    nthreads=max([threads_available, 6]))


sc.tl.umap(adata)
L.info("done umap, saving stuff")
# write out
umap = pd.DataFrame(adata.obsm['X_umap'], adata.obs.index)
umap.to_csv(args.output_csv)

# save the scanorama dim reduction in case scanorama is our favourite
adata.write("tmp/scanorama_scaled_adata.h5ad")

L.info("done")

