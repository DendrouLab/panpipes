import pandas as pd
import scanpy as sc
import muon as mu
import argparse
# import h5py
import os

from panpipes.funcs.processing import check_for_bool
from panpipes.funcs.io import read_anndata, write_anndata
from panpipes.funcs.scmethods import run_neighbors_method_choice

import multiprocessing 
threads_available = multiprocessing.cpu_count()

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

parser.add_argument('--input_anndata',
                    default='adata_log1p.h5ad',
                    help='')
parser.add_argument('--output_csv', default='batch_correction/umap_bc_combat.csv',
                    help='')
parser.add_argument('--integration_col', default='batch', help='')
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
parser.add_argument('--modality',
                    help="")


args, opt = parser.parse_known_args()



L.info("reading data and starting integration pipeline with script: ")
L.info(os.path.basename(__file__))
L.info("Running with options: %s", args)
L.info("threads available: %s", threads_available)
# this should be an object that contains the full log normalised data (adata_log1p.h5ad)
# prior to hvgs and filtering
#adata = read_anndata(args.input_anndata, use_muon=use_muon, modality="rna")
mdata = mu.read(args.input_anndata)
adata = mdata.mod[args.modality] 

# combat can't integrate on 2+ variables, so create a fake column with combined information
columns = [x.strip() for x in args.integration_col.split(",")]
if len(columns) > 1: 
    L.info("using 2 columns to integrate on more variables")
    #comb_columns = "_".join(columns)
    adata.obs["comb_columns"] = adata.obs[columns].apply(lambda x: '|'.join(x), axis=1)

    # make sure that batch is a categorical
    adata.obs["comb_columns"] = adata.obs["comb_columns"].astype("category")
    # run combat
    sc.pp.combat(adata, key="comb_columns")
    
else:

    # run combat
    sc.pp.combat(adata, key=args.integration_col)

L.info("integration run now calculate umap")
sc.pp.pca(adata, use_highly_variable=True, svd_solver='arpack') #use same args as plot_pca.py

# run neighbours and umap 
run_neighbors_method_choice(adata, 
    method=args.neighbors_method, 
    n_neighbors=int(args.neighbors_k), 
    n_pcs=int(args.neighbors_n_pcs), 
    metric=args.neighbors_metric, 
    use_rep='X_pca',
    nthreads=max([threads_available, 6]))


sc.tl.umap(adata)
L.info("done umap, saving stuff")
#write out
umap = pd.DataFrame(adata.obsm['X_umap'], adata.obs.index)
umap.to_csv(args.output_csv)

adata.write("tmp/combat_scaled_adata_" + args.modality + ".h5ad")


L.info("Done")

