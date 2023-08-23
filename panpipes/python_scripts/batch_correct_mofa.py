import multiprocessing 

import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import os
import gc
import muon as mu
from cgatcore import pipeline as P
import panpipes.funcs as pp
from panpipes.funcs.processing import check_for_bool
from panpipes.funcs.io import read_anndata, write_anndata
from panpipes.funcs.scmethods import run_neighbors_method_choice, X_is_raw

import sys
import logging
L = logging.getLogger()
L.setLevel(logging.INFO)
log_handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s: %(levelname)s - %(message)s')
log_handler.setFormatter(formatter)
L.addHandler(log_handler)

# load arguments 
parser = argparse.ArgumentParser()
parser.add_argument('--scaled_anndata',
                    default='adata_scaled.h5ad',
                    help='a preprocessed mudata/anndata')
parser.add_argument('--output_csv', default='batch_correction/umap_bc_MultiVI.csv',
                    help='')
parser.add_argument('--use_gpu', default=False,
                    help='')
parser.add_argument('--integration_col_categorical', default='batch',
                    help='')
parser.add_argument('--n_factors', default="",
                    help='number of factors to train the model with - optional')
parser.add_argument('--n_iterations', default='',
                    help='upper limit on the number of iterations')
parser.add_argument('--convergence_mode', default='fast',
                    help='fast, medium or slow')
parser.add_argument('--save_parameters', default=bool,
                    help='whether to save parameters models')
parser.add_argument('--outfile_model', default="",
                    help="if args.save_parameters true, path to HDF5 file to store the model")
parser.add_argument('--figdir', default='./figures',
                    help='')
parser.add_argument('--neighbors_n_pcs',
                    help="n_pcs", default=50)
parser.add_argument('--neighbors_method', default="scanpy",
                    help="neighbours method, scanpy or hnsw")
parser.add_argument('--neighbors_k', default=30,
                    help="neighbors k")
parser.add_argument('--neighbors_metric',default="euclidean",
                    help="neighbor metric, e.g. euclidean or cosine")

args, opt = parser.parse_known_args()

L.info(args)
# scanpy settings
sc.set_figure_params(facecolor="white")
sc.settings.autoshow = False
sc.settings.figdir = args.figdir


# load parameters

threads_available = multiprocessing.cpu_count()
params = pp.io.read_yaml("pipeline.yml")

mdata = mu.read(args.scaled_anndata)

if params['multimodal']['mofa']['modalities'] is not None:
    modalities= params['multimodal']['mofa']['modalities']
    modalities = [x.strip() for x in modalities.split(",")]
    L.info(f"using modalities :{modalities}")
    removed_mods = None
    if all(x in modalities for x in mdata.mod.keys()):
        tmp = mdata.copy()
        L.info('using all modalities')
    else:
        tmp = mdata.copy()
        removed_mods = list(set(mdata.mod.keys()) - set(modalities))
        L.info(f"removing modalities {removed_mods}")
        for rmod in removed_mods:
            del tmp.mod[rmod]
else:
    L.warning("""you specified no modalities, so i will default to all available
                this may be a problem if you have repertoire in here""")
    removed_mods = None
    tmp = mdata.copy()  

print(tmp.mod.keys())

mu.pp.intersect_obs(tmp)

mofa_kwargs={}
#expected args:
# n_factors: 10
# n_iterations: 1000
# convergence_mode: fast, medium, slow
# save_parameters: False
# #if save_parameters True, set the following, otherwise leave blank
# outfile_model:

mofa_kwargs = params['multimodal']['mofa']
del mofa_kwargs['modalities']

if mofa_kwargs['filter_by_hvg']:
    mofa_kwargs['use_var'] = "highly_variable"
    del mofa_kwargs['filter_by_hvg']
    for mod in tmp.mod.keys():
        if "highly_variable" not in tmp[mod].var.columns:
            print(mod)
            tmp[mod].var["highly_variable"] = True
        
    tmp.update()
    


if args.integration_col_categorical is not None:
    if args.integration_col_categorical in tmp.obs.columns:
        mofa_kwargs['groups_label'] = args.integration_col_categorical

# default is to read yaml and parse directly to kwargs.
# if the defaults expected params are parsed by the script in some other way
# they will overwrite the initial reading of the yml

if mofa_kwargs['save_parameters'] is None:
    if args.save_parameters is not None:
        mofa_kwargs['save_parameters'] = check_for_bool(args.save_parameters)
        if args.outfile_model is not None:
            mofa_kwargs['outfile'] = args.outfile_model

if mofa_kwargs['n_factors'] is None:
    if args.n_factors is not None:
        mofa_kwargs['n_factors'] = int(args.n_factors)

if mofa_kwargs['n_iterations'] is None:
    if args.n_iterations is not None:
        mofa_kwargs['n_iterations'] = int(args.n_iterations)

if mofa_kwargs['convergence_mode'] is None:
    if args.convergence_mode is not None:
        mofa_kwargs['convergence_mode'] = args.convergence_mode
        
# are we using the gpu?
if args.use_gpu is True or args.use_gpu=='True':
     mofa_kwargs['gpu_mode'] = True

L.info(mofa_kwargs)
L.info(tmp.var.columns)


mu.tl.mofa(tmp, **mofa_kwargs)
#This adds X_mofa embeddings to the .obsm slot of the MuData object
#write the discovered latent rep to the original mudata object
mdata.obsm['X_mofa'] = tmp.obsm['X_mofa'].copy()   

if int(args.neighbors_n_pcs) > mdata.obsm['X_mofa'].shape[1]:
    L.warning(f"N PCs is larger than X_mofa dimensions, reducing n PCs to  {mdata.obsm['X_mofa'].shape[1]-1}")
n_pcs= min(int(args.neighbors_n_pcs), mdata.obsm['X_mofa'].shape[1]-1)


run_neighbors_method_choice(mdata, 
    method=args.neighbors_method, 
    n_neighbors=int(args.neighbors_k), 
    n_pcs=n_pcs, #this should be the # rows of var, not obs ???????
    metric=args.neighbors_metric, 
    use_rep='X_mofa',
    nthreads=max([threads_available, 6]))

sc.tl.umap(mdata, min_dist=0.4)
sc.tl.leiden(mdata,  key_added="leiden_mofa")

umap = pd.DataFrame(mdata.obsm['X_umap'], mdata.obs.index)
umap.to_csv(args.output_csv)

if removed_mods is not None:
    for rmd in removed_mods:
        tmp.mod[rmd] = mdata.mod[rmd].copy()

tmp.write("tmp/mofa_scaled_adata.h5mu")


L.info("Done")

