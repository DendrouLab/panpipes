# check for threads number
import multiprocessing 
threads_available = multiprocessing.cpu_count()

# import numpy as np
import pandas as pd
import scanpy as sc
import os
import argparse
from muon import MuData
from panpipes.funcs.processing import check_for_bool
from panpipes.funcs.io import read_anndata, write_anndata
from panpipes.funcs.scmethods import run_neighbors_method_choice
import muon as mu
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
                    default='adata_scaled.h5ad',
                    help='')
parser.add_argument('--modality',
                    default='rna',
                    help='')
parser.add_argument('--dimred',
                    default='PCA',
                    help='which dimred to expect, relevant for ATAC')
parser.add_argument('--output_csv', default='batch_correction/umap_bc_none.csv',
                    help='')
parser.add_argument('--n_threads', default=1,
                    help="num threads to use for neighbor computations")
parser.add_argument('--integration_col', default='batch')
parser.add_argument('--neighbors_n_pcs',
                    help="n_pcs")
parser.add_argument('--neighbors_method',
                    help="neighbours method, scanpy or hnsw")
parser.add_argument('--neighbors_k',
                    help="neighbors k")
parser.add_argument('--neighbors_metric',
                    help="neighbor metric, e.g. euclidean or cosine")

args, opt = parser.parse_known_args()

L.info("reading data and starting integration pipeline with script: ")
L.info(os.path.basename(__file__))


# Scanorama is designed to be used in scRNA-seq pipelines downstream of noise-reduction methods,
# including those for imputation and highly-variable gene filtering.

# adata = sc.read(args.input_anndata)
#adata = read_anndata(args.input_anndata, use_muon=True, modality=args.modality)
adata = mu.read(args.input_anndata +"/" + args.modality)
# adata = mdata.mod[args.modality] 



columns = [x.strip() for x in args.integration_col.split(",")]
if len(columns)>1: 
    comb_columns = "|".join(columns)
    adata.obs[comb_columns] = adata.obs[columns].apply(lambda x: '|'.join(x), axis=1)
    columns += [comb_columns]

# write out batch
# adata.obs[columns].to_csv(os.path.join(os.path.dirname(args.output_csv), 'batch_'+ args.modality +'_mtd.csv'))
if args.dimred == "PCA":
    dimred = "X_pca"
elif args.dimred == "LSI":
    dimred = "X_lsi"


# run neighbours and umap without batch correction
pc_kwargs = {}
if int(args.neighbors_n_pcs) > 0:
    pc_kwargs['use_rep'] = dimred
    pc_kwargs['n_pcs'] = int(args.neighbors_n_pcs)
else:
    # If n_pcs==0 use .X if use_rep is None.
    pc_kwargs['use_rep'] = None
    pc_kwargs['n_pcs'] = int(args.neighbors_n_pcs)
    
run_neighbors_method_choice(adata, 
    method=args.neighbors_method, 
    n_neighbors=int(args.neighbors_k), 
    metric=args.neighbors_metric, 
    nthreads=max([threads_available, 6]), **pc_kwargs)

L.info("done n_neighbours, saving stuff")

sc.tl.umap(adata)

L.info("done umap, saving stuff")
#write out
umap = pd.DataFrame(adata.obsm['X_umap'], adata.obs.index)
umap.to_csv(args.output_csv)

L.info("done")


outfiletmp = ("tmp/no_correction_scaled_adata_" + args.modality + ".h5ad" )

L.info("saving harmony corrected adata")
write_anndata(adata, outfiletmp, use_muon=False, modality=args.modality)

L.info("done")

