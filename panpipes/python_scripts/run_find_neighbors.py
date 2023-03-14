import argparse
import os
import sys
import logging
import scanpy as sc
from muon import MuData, read
from panpipes.funcs.scmethods import run_neighbors_method_choice

L = logging.getLogger()
L.setLevel(logging.INFO)
log_handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s: %(levelname)s - %(message)s')
log_handler.setFormatter(formatter)
L.addHandler(log_handler)


parser = argparse.ArgumentParser()
parser.add_argument('--infile', default='./anndata.h5mu',
                    help="file name, format: anndata_scaled.h5ad")
parser.add_argument('--modality',
                    default="rna",
                    help='')
parser.add_argument('--outfile', default='./anndata_wneighbors.h5ad',
                    help="file name, format: .h5ad")
parser.add_argument('--integration_method', default=None,
                    help="helps find the correct dimension reduction for sc.pp.neighbors(), default='X_pca'")
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
args, opt = parser.parse_known_args()


sc.settings.n_jobs = int(args.n_threads)
L.info("Running with options: %s", args)

# read data
L.info("reading adata")
mdata = read(args.infile)
if type(mdata) is MuData:
    adata=mdata[args.modality]

if "X_pca" not in adata.obsm.keys():
    sc.tl.pca(adata)

# run command
opts = dict(method=args.neighbors_method,
            n_neighbors=int(args.neighbors_k),
            n_pcs=int(args.neighbors_n_pcs),
            metric=args.neighbors_metric,
            nthreads=args.n_threads)

dim_red_location = {
    "None": "X_pca",
    "combat": "X_pca",
    "scanorama": "X_scanorama",
    "harmony" : "X_harmony" ,
    "scvi": "X_scVI",
    "totalvi" : "X_totalVI",

}
opts['use_rep'] = dim_red_location[str(args.integration_method)]

if args.integration_method == 'bbknn':
    L.info("bbknn neighbours should already be included in the object, not rerunning.")
    try:
        adata.uns['neighbors']
    except KeyError:
        sys.exit("bbknn neighbours not found, make sure they are in the scaled object")
    pass
else:
    run_neighbors_method_choice(adata,
                                **opts)



L.info("saving data")
if os.path.exists(os.path.basename(args.outfile)) is False:
    os.makedirs(os.path.basename(args.outfile))
adata.write(args.outfile)


L.info("Done")

