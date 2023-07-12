import pandas as pd
import scanpy as sc
import argparse
import os
import harmonypy as hm
import muon as mu
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
                    default='adata_scaled.h5ad',
                    help='')
parser.add_argument('--modality',
                    default='rna',
                    help='')
parser.add_argument('--dimred',
                    default='PCA',
                    help='which dimred to expect, relevant for ATAC')
parser.add_argument('--output_csv', default='batch_correction/umap_bc_harmony.csv',
                    help='')
parser.add_argument('--integration_col',
                    help='')
parser.add_argument('--n_threads', default=1,
                    help="num threads to use for neighbor computations")
parser.add_argument('--harmony_npcs', default=30,
                    help="npcs for running harmony")        
parser.add_argument('--sigma_val', default=0.1,
                    help="sigma")
parser.add_argument('--theta_val', default=1.0,
                    help="theta")                   
parser.add_argument('--neighbors_n_pcs', default=30,
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
L.info("Running with options: %s", args)

# this should be an object that contains the full log normalised data (adata_log1p.h5ad)
# prior to hvgs and filtering
#adata = read_anndata(args.input_anndata, use_muon=use_muon, modality=args.modality)
mdata = mu.read(args.input_anndata)
adata = mdata.mod[args.modality] 

# Harmony can integrate on 2+ variables,
# but for consistency with other approaches create a fake column with combined information
columns = [x.strip() for x in args.integration_col.split(",")]
if args.dimred == "PCA":
    dimred = "X_pca"
elif args.dimred == "LSI":
    dimred = "X_lsi"

if len(columns)>1: 
    L.info("using 2 columns to integrate on more variables")
    #comb_columns = "_".join(columns)
    adata.obs["comb_columns"] = adata.obs[columns].apply(lambda x: '|'.join(x), axis=1)

    # make sure that batch is a categorical
    adata.obs["comb_columns"] = adata.obs["comb_columns"].astype("category")
    # run harmony

    ho = hm.run_harmony(adata.obsm[dimred][:,0:int(args.harmony_npcs)], adata.obs, ["comb_columns"], 
                                       sigma = float(args.sigma_val),theta = float(args.theta_val),verbose=True,max_iter_kmeans=30, 
                                       max_iter_harmony=40)

else:
    # make sure that batch is a categorical
    adata.obs[args.integration_col] = adata.obs[args.integration_col].astype("category")
    # run harmony
    ho = hm.run_harmony(adata.obsm[dimred][:,0:int(args.harmony_npcs)],
                        adata.obs,
                        [args.integration_col],
                        sigma=float(args.sigma_val),
                        theta = float(args.theta_val),
                        verbose=True,max_iter_kmeans=30,
                        max_iter_harmony=40)


L.info("integration run now calculate umap")

adjusted_pcs = pd.DataFrame(ho.Z_corr).T
adata.obsm['X_harmony']=adjusted_pcs.values
L.info("harmony co-ords derived")

if int(args.neighbors_n_pcs) >adata.obsm['X_harmony'].shape[1]:
    L.warn(f"N PCs is larger than X_harmony dimensions, reducing n PCs to  {adata.obsm['X_harmony'].shape[1] -1}")

n_pcs= min(int(args.neighbors_n_pcs), adata.obsm['X_harmony'].shape[1]-1)

# run neighbours and umap 
run_neighbors_method_choice(adata, 
    method=args.neighbors_method, 
    n_neighbors=int(args.neighbors_k), 
    n_pcs=n_pcs, 
    metric=args.neighbors_metric, 
    use_rep='X_harmony',
    nthreads=max([threads_available, 6]))


L.info("done n_neighbours, moving to umap")
sc.tl.umap(adata)
L.info("done umap, saving stuff")
# write out
umap = pd.DataFrame(adata.obsm['X_umap'], adata.obs.index)
umap.to_csv(args.output_csv)


#adata.write("tmp/harmony_scaled_adata_" + args.modality + ".h5ad")


outfiletmp = ("tmp/harmony_scaled_adata_" + args.modality + ".h5ad" )

L.info("saving harmony corrected adata")
write_anndata(adata, outfiletmp, use_muon=False, modality=args.modality)


L.info("Done")


