import argparse
import pandas as pd
import scanpy as sc
import numpy as np

from panpipes.funcs.processing import extract_parameter_from_fname, check_for_bool
from panpipes.funcs.io import read_anndata
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
                    default="/well/dendrou/projects/cart/data/by_sample/pipe_qc/pipe_int_gut_pctmt50/pipe_clust_bbknn/gut_nneigh15_pcs50_dir/gut_bcbbknn_meteuclidean_neighbors.h5ad",
                    help='')
parser.add_argument('--output_anndata', default='adata_cellxgene.h5ad',
                    help='')
parser.add_argument('--use_muon',
                    default=False,
                    help='')
parser.add_argument('--modality',
                    default="rna",
                    help='')
parser.add_argument('--sample_prefix', default='adata_cellxgene.h5ad',
                    help='')               
parser.add_argument('--clusters', default="/well/dendrou/projects/cart/data/by_sample/pipe_qc/pipe_int_gut_pctmt50/pipe_clust_bbknn/gut_nneigh15_pcs50_dir/all_res_clusters_list.txt.gz",
                    help='')
parser.add_argument('--umaps', default="/well/dendrou/projects/cart/data/by_sample/pipe_qc/pipe_int_gut_pctmt50/pipe_clust_bbknn/gut_nneigh15_pcs50_dir/gut_md0.1_umap.txt.gz,/well/dendrou/projects/cart/data/by_sample/pipe_qc/pipe_int_gut_pctmt50/pipe_clust_bbknn/gut_nneigh15_pcs50_dir/gut_md0.5_umap.txt.gz",
                    help="")
parser.add_argument('--best_md', default=None)
parser.add_argument('--best_cluster_col', default=None)



args, opt = parser.parse_known_args()
L.info(args)
use_muon = check_for_bool(args.use_muon)

adata = read_anndata(args.infile, use_muon=use_muon, modality=args.modality)


# add in umap coordinates
umap_files = args.umaps.split(",")
umap_names = ["X_umap_md" + str(extract_parameter_from_fname(x, "md", "")) for x in umap_files]

for i in range(len(umap_files)):
    umap_coords = pd.read_csv(umap_files[i], sep='\t', index_col=0)
    if all(umap_coords.index == adata.obs_names):
        adata.obsm[umap_names[i]] = np.array(umap_coords)
    else:
        L.error("umap barcodes do not match anndata")

# add in clusters
clusters = pd.read_csv(args.clusters, sep='\t', index_col=0)
clusters = clusters.apply(lambda x: x.astype('category'))
if all(clusters.index == adata.obs_names):
    obs = adata.obs.merge(clusters, left_index=True, right_index=True)
    if all(obs.index == adata.obs_names):
        adata.obs = obs
else:
    L.error("cluster barcodes do not match anndata")

# move log counts to a layer
adata.layers['log_normcounts'] = adata.raw.X.copy()
adata.raw = None

# set favourtie params
if args.best_md is not None:
    try:
        adata.obsm['X_umap_int'] = adata.obsm['X_umap']
    except KeyError:
        L.info("integration umap not found, X_umap_int not created")
    adata.obsm['X_umap'] = adata.obsm['X_umap_md' + str(args.best_md)]

if args.best_cluster_col is not None:
    adata.obs['clusters'] = adata.obs[args.best_cluster_col]


adata.write(args.output_anndata)


L.info("Done")

