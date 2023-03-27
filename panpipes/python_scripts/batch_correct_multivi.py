import multiprocessing 

import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import scvi
import argparse
import os
import gc
import muon as mu
from cgatcore import pipeline as P
import anndata as ad
import panpipes.funcs as pp
import gc
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
                    help='')
parser.add_argument('--output_csv', default='batch_correction/umap_bc_MultiVI.csv',
                    help='')
parser.add_argument('--integration_col_categorical', default='batch',
                    help='')
parser.add_argument('--integration_col_continuous', default=None,
                    help='')
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

# scanpy settings
sc.set_figure_params(facecolor="white")
sc.settings.autoshow = False
sc.settings.figdir = args.figdir


# load parameters

threads_available = multiprocessing.cpu_count()
params = pp.io.read_yaml("pipeline.yml")

test_script=False
if test_script:
    L.info("this is a test run")
    params['MultiVI']['training_args']['max_epochs'] = 10

# ------------------------------------------------------------------
L.info("Running multivi")

# mdata = mu.read(args.scaled_anndata)
# rna = mdata['rna'].copy()
# atac = mdata['atac'].copy()
rna= mu.read(args.scaled_anndata +"/" + "rna")
atac= mu.read(args.scaled_anndata +"/" + "atac")


if check_for_bool(params["multimodal"]["MultiVI"]["lowmem"]):
    L.info("subsetting atac to top 25k HVF")
    sc.pp.highly_variable_genes(atac, n_top_genes=25000)
    atac = atac[:, atac.var.highly_variable]



gc.collect()

# if "modality" not in atac.var.keys():
#     atac.var["modality"] = "Peaks"

# if "modality" not in rna.var.keys():
#     rna.var["modality"] = "Gene Expression"

if rna.shape[0] == atac.shape[0]:
    n=int(rna.shape[0])
else:
    sys.exit("rna and atac have different number of cells, \
        can't deal with this in this version of MultiVI integration")

gc.collect()

n_genes = len(rna.var_names)
n_regions = len(atac.var_names)


if "raw_counts" in rna.layers.keys():
    L.info("raw counts found in layer for RNA")
elif X_is_raw(rna):
    # this means the X layer is already raw and we can make the layer we need
    L.info("raw counts found in X for RNA")
    rna.layers["raw_counts"] = rna.X.copy()
else:
    sys.exit("please provide a mu/anndata with raw counts for RNA")

if "raw_counts" in atac.layers.keys():
    L.info("raw counts found in layer for ATAC")
elif X_is_raw(atac):
    # this means the X layer is already raw and we can make the layer we need
    L.info("raw counts found in X for ATAC")
    atac.layers["raw_counts"] = atac.X.copy()
else:
    sys.exit("please provide a mu/anndata with raw counts for ATAC")

L.info("concatenating modalities to comply with multiVI")
# adata_paired = ad.concat([rna, atac], join="outer")
# adata_paired.var = pd.concat([rna.var,atac.var])

adata_paired = ad.concat([rna.T, atac.T]).T

rna_cols=rna.obs.columns
atac_cols=atac.obs.columns

rnaobs = rna.obs.copy()
rnaobs.columns= ["rna:"+ x for x in rna_cols]
atacobs = atac.obs.copy()
atacobs.columns= ["atac:"+ x for x in atac_cols]
adata_paired.obs = pd.merge(rnaobs, atacobs, left_index=True, right_index=True)

if "modality" not in adata_paired.obs.columns:
    adata_paired.obs["modality"] = "paired" 

print(adata_paired)

# adata = ad.concat([rna,atac],join="outer")
# adata.var = pd.concat([rna.var,atac.var])

del [rna , atac ]
gc.collect()

adata_mvi = scvi.data.organize_multiome_anndatas(adata_paired)


# MultiVI integrates by modality, to use batch correction you need a batch covariate to specify in
# categorical_covariate_keys
if args.integration_col_categorical is not None :
    cols = [x.strip() for x in args.integration_col_categorical.split(",")]
    columns = []
    for cc in cols:
        if cc in rna_cols:
            print("categorical batch covariate in rna")
            columns.append("rna:"+cc)
        elif cc in atac_cols:
            print("categorical batch covariate in atac")
            columns.append("atac:"+cc)
if args.integration_col_continuous is not None :
    if args.integration_col_continuous in rna_cols:
        args.integration_col_continuous = "rna:"+ args.integration_col_continuous
    elif args.integration_col_continuous in atac_cols:
        args.integration_col_continuous = "atac:"+ args.integration_col_continuous


kwargs = {}
# in case of more than 1 variable, create a fake column with combined information
if columns is not None:
    print(columns)
    if len(columns) > 1:
        L.info("using 2 columns to integrate on more variables")
        # bc_batch = "_".join(columns)
        adata_mvi.obs["bc_batch"] = adata_mvi.obs[columns].apply(lambda x: '|'.join(x), axis=1) 
        # make sure that batch is a categorical
        adata_mvi.obs["bc_batch"] = adata_mvi.obs["bc_batch"].astype("category")
    else:
        adata_mvi.obs['bc_batch'] = adata_mvi.obs[columns]
        adata_mvi.obs["bc_batch"] = adata_mvi.obs["bc_batch"].astype("category")

    batch_categories = list(adata_mvi.obs['bc_batch'].unique())
    kwargs['categorical_covariate_keys'] = ["bc_batch"]

if args.integration_col_continuous is not None :
    print(args.integration_col_continuous)
    adata_mvi.obs['bc_batch_continuous'] = adata_mvi.obs[args.integration_col_continuous]
    kwargs['continuous_covariate_keys'] = ["bc_batch_continuous"]

print(kwargs)

# 1 setup anndata
#scvi.model.MULTIVI.setup_anndata(adata_mvi, batch_key='modality', **kwargs)
scvi.model.MULTIVI.setup_anndata(
    adata_mvi, 
    batch_key='modality',
    layer =  "raw_counts",
    **kwargs) 
#2. setup model
if params["multimodal"]['MultiVI']['model_args'] is None:
    multivi_model_args =  {}
else:
    multivi_model_args =  {k: v for k, v in params["multimodal"]['MultiVI']['model_args'].items() if v is not None}

print(multivi_model_args)

mvi = scvi.model.MULTIVI(
    adata_mvi,
    n_genes=n_genes,
    n_regions=n_regions,
    **multivi_model_args
)

#3.train

if params["multimodal"]["MultiVI"]["training_args"] is None:
    multivi_training_args={}
else:
    multivi_training_args =  {k: v for k, v in params["multimodal"]['MultiVI']['training_args'].items() if v is not None}

print(multivi_training_args)

if params["multimodal"]['MultiVI']['training_plan'] is None:
    multivi_training_plan = {}
else:
    multivi_training_plan =  {k: v for k, v in params["multimodal"]['MultiVI']['training_plan'].items() if v is not None}

print(multivi_training_plan)

mvi.view_anndata_setup()

mvi.train( **multivi_training_args, **multivi_training_plan)

mvi.save(os.path.join("batch_correction", "multivi_model"), 
                  anndata=False)

plt.plot(mvi.history["elbo_train"], label="train")
plt.plot(mvi.history["elbo_validation"], label="validation")
plt.title("Negative ELBO over training epochs")

plt.legend()
plt.savefig(os.path.join(args.figdir, "multivi_elbo_plot.png"))

L.info("""We support the use of mudata as a general framework for multimodal data
        For this reason, the object we save is not the classical anndata
        you would find in scvitools tutorial.
        We chose to save the learned SC representation in the 
        Mudata obsm slot, and any single modality processing in
        its own modality slot
            """)

mdata = mu.read(args.scaled_anndata)
mdata.obsm["X_MultiVI"] = mvi.get_latent_representation()
#adata_mvi.obsm["X_MultiVI"] = mvi.get_latent_representation()

if int(args.neighbors_n_pcs) > mdata.obsm['X_MultiVI'].shape[1]:
    L.warn(f"N PCs is larger than X_MultiVI dimensions, reducing n PCs to  {mdata.obsm['X_MultiVI'].shape[1] -1}")
n_pcs= min(int(args.neighbors_n_pcs), mdata.obsm['X_MultiVI'].shape[1]-1)

run_neighbors_method_choice(mdata, 
    method=args.neighbors_method, 
    n_neighbors=int(args.neighbors_k), 
    n_pcs=n_pcs, 
    metric=args.neighbors_metric, 
    use_rep='X_MultiVI',
    nthreads=max([threads_available, 6]))
sc.tl.umap(mdata, min_dist=0.4)
sc.tl.leiden(mdata, key_added="leiden_multiVI")

umap = pd.DataFrame(mdata.obsm['X_umap'], mdata.obs.index)
umap.to_csv(args.output_csv)

write_anndata(mdata, "tmp/multivi_scaled_adata.h5mu",use_muon=False)

L.info("Done")

