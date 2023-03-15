import argparse
import os
import sys
import logging
import scanpy as sc
from muon import MuData, read
from panpipes.funcs.scmethods import run_neighbors_method_choice
from panpipes.funcs.io import read_yaml

L = logging.getLogger()
L.setLevel(logging.INFO)
log_handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s: %(levelname)s - %(message)s')
log_handler.setFormatter(formatter)
L.addHandler(log_handler)


parser = argparse.ArgumentParser()
parser.add_argument('--infile', default='./anndata.h5mu',
                    help="file name, format: anndata_scaled.h5ad")
parser.add_argument('--outfile', default='./anndata_wneighbors.h5ad',
                    help="file name, format: .h5ad")
parser.add_argument('--neighbor_dict', default=None,
                    help="helps find the correct dimension reduction for sc.pp.neighbors(), default='X_pca'")
parser.add_argument('--n_threads', default=1,help='number of threads available')

args, opt = parser.parse_known_args()
L.info(args)
# load the filtering dictionary 
neighbor_dict = args.neighbor_dict
if isinstance(args.neighbor_dict, dict):
    neighbor_dict = args.neighbor_dict
else:
    neighbor_dict = read_yaml(args.neighbor_dict) 

sc.settings.n_jobs = int(args.n_threads)
L.info("Running with options: %s", args)

# read data
L.info("reading mudata")
mdata = read(args.infile)



for mod in neighbor_dict.keys():
    if neighbor_dict[mod]['use_existing']:
        L.info('using existing neighbors graph for %s' % mod)
        pass
    else:
        L.info("computing new neighbors for %s" % mod)
        if type(mdata) is MuData:
            adata=mdata[mod]
        if (neighbor_dict[mod]['dim_red'] == "X_pca") and ("X_pca" not in adata.obsm.keys()):
            L.info("X_pca not found, computing it using default parameters")
            sc.tl.pca(adata)
        # run command
        opts = dict(method=neighbor_dict[mod]['method'],
                    n_neighbors=int(neighbor_dict[mod]['k']),
                    n_pcs=int(neighbor_dict[mod]['n_dim_red']),
                    metric=neighbor_dict[mod]['metric'],
                    nthreads=args.n_threads,
                    use_rep=neighbor_dict[mod]['dim_red'])


        run_neighbors_method_choice(adata,**opts)
        mdata.update()



L.info("saving data")
mdata.write(args.outfile)
L.info("done")