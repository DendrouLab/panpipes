import argparse
import os
import sys
import logging
import scanpy as sc
from muon import MuData, read

from panpipes.funcs.scmethods import run_neighbors_method_choice
from panpipes.funcs.io import read_yaml
from panpipes.funcs.scmethods import lsi

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
L.info("Running with params: %s", args)
# load the filtering dictionary 
neighbor_dict = args.neighbor_dict
if isinstance(args.neighbor_dict, dict):
    neighbor_dict = args.neighbor_dict
else:
    neighbor_dict = read_yaml(args.neighbor_dict) 

sc.settings.n_jobs = int(args.n_threads)

# read data
if ".zarr" in args.infile:
    import spatialdata as sd
    L.info("Reading in SpatialData from '%s'" % args.infile)
    sdata = sd.read_zarr(args.infile)
else: 
    L.info("Reading in MuData from '%s'" % args.infile)
    mdata = read(args.infile)



for mod in neighbor_dict.keys():
    if mod != "spatial": 
        if mod in mdata.mod.keys():
            if neighbor_dict[mod]['use_existing']:
                L.info('Using existing neighbors graph for %s' % mod)
                pass
            else:
                L.info("Computing new neighbors for modality %s on %s" % (mod, neighbor_dict[mod]['dim_red']))
                if type(mdata) is MuData:
                    adata=mdata[mod]
                if (neighbor_dict[mod]['dim_red'] == "X_pca") and ("X_pca" not in adata.obsm.keys()):
                    L.info("X_pca not found, computing it using default parameters")
                    sc.tl.pca(adata)
                    if (mod == "atac") and (neighbor_dict[mod]['dim_remove'] is not None):
                        dimrem = int(neighbor_dict[mod]['dim_remove'])
                        adata.obsm['X_pca'] = adata.obsm['X_pca'][:, dimrem:]
                        adata.varm["PCs"] = adata.varm["PCs"][:, dimrem:]
                if mod == "atac":
                    if (neighbor_dict[mod]['dim_red'] == "X_lsi") and ("X_lsi" not in adata.obsm.keys()):
                        L.info("X_lsi not found, computing it using default parameters")
                        lsi(adata=adata, num_components=50)
                        if neighbor_dict[mod]['dim_remove'] is not None:
                            L.info("Removing dimension %s from X_lsi" % neighbor_dict[mod]['dim_remove'])
                            dimrem = int(neighbor_dict[mod]['dim_remove'])
                            adata.obsm['X_lsi'] = adata.obsm['X_lsi'][:, dimrem:]
                            adata.varm["LSI"] = adata.varm["LSI"][:, dimrem:]
                            adata.uns["lsi"]["stdev"] = adata.uns["lsi"]["stdev"][dimrem:]

                # run command
                opts = dict(method=neighbor_dict[mod]['method'],
                            n_neighbors=int(neighbor_dict[mod]['k']),
                            n_pcs=int(neighbor_dict[mod]['n_dim_red']),
                            metric=neighbor_dict[mod]['metric'],
                            nthreads=args.n_threads,
                            use_rep=neighbor_dict[mod]['dim_red'])


                run_neighbors_method_choice(adata,**opts)
                mdata.mod[mod] = adata
                mdata.update()
    else:  
        if neighbor_dict[mod]['use_existing']:
            L.info('Using existing neighbors graph for %s' % mod)
            pass
        else:
            L.info("Computing new neighbors for modality %s on %s" % (mod, neighbor_dict[mod]['dim_red']))
            if (neighbor_dict[mod]['dim_red'] == "X_pca") and ("X_pca" not in sdata["table"].obsm.keys()):
                L.info("X_pca not found, computing it using default parameters")
                sc.tl.pca(sdata["table"])
            opts = dict(method=neighbor_dict[mod]['method'],
                        n_neighbors=int(neighbor_dict[mod]['k']),
                        n_pcs=int(neighbor_dict[mod]['n_dim_red']),
                        metric=neighbor_dict[mod]['metric'],
                        nthreads=args.n_threads,
                        use_rep=neighbor_dict[mod]['dim_red'])
            # run command
            run_neighbors_method_choice(sdata["table"],**opts)


if ".zarr" in args.infile:
    L.info("Saving updated SpatialData to '%s'" % args.outfile)
    sdata.write(args.outfile)
else: 
    L.info("Saving updated MuData to '%s'" % args.outfile)
    mdata.write(args.outfile)

L.info("Done")
