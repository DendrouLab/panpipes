import multiprocessing 

import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import scvi
import argparse
import os
import gc
import muon as mu

from panpipes.funcs.io import read_anndata, read_yaml
from panpipes.funcs.scmethods import run_neighbors_method_choice, X_is_raw
from panpipes.funcs.processing import intersect_obs_by_mod

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
parser.add_argument('--raw_anndata',
                    default=None,
                    help='')
parser.add_argument('--output_csv', default='batch_correction/umap_bc_totalvi.csv',
                    help='')
parser.add_argument('--integration_col_categorical', default=None,
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
parser.add_argument('--scvi_seed',default=None,
                    help="set explicitly seed to make runs reproducible")
    
    


args, opt = parser.parse_known_args()
L.info("Running with params: %s", args)
# scanpy settings
sc.set_figure_params(facecolor="white")
sc.settings.autoshow = False
sc.settings.figdir = args.figdir


if args.scvi_seed is not None:
    scvi.settings.seed = int(args.scvi_seed)
else:
    scvi.settings.seed = 1492

# load parameters

threads_available = multiprocessing.cpu_count()

params = read_yaml("pipeline.yml")
params['sample_prefix']



test_script=False
if test_script:
    L.info("this is a test run")
    params['multimodal']['totalvi']['training_args']['max_epochs'] = 10

# ------------------------------------------------------------------
L.info("Reading in MuData from '%s'" % args.scaled_anndata)
mdata = mu.read(args.scaled_anndata)

# we need to intersect the prot and rna modalities.
L.info("Intersecting prot and rna modalities")
intersect_obs_by_mod(mdata, ['rna', 'prot'])

# this is a copy so we can subset vars and obs without changing the original object
rna = mdata['rna'].copy()
prot = mdata['prot'].copy()

kwargs={}
# in case of more than 1 variable, create a fake column with combined information
if args.integration_col_categorical is not None :
    columns = [x.strip() for x in args.integration_col_categorical.split(",")]
    columns =[ x.replace("rna:","") for x in columns]
    if len(columns) > 1:
        L.info("Using 2 columns to integrate on more variables")
        # bc_batch = "_".join(columns)
        rna.obs["bc_batch"] = rna.obs[columns].apply(lambda x: '|'.join(x), axis=1)
        # make sure that batch is a categorical
        rna.obs["bc_batch"] = rna.obs["bc_batch"].astype("category")
    else:
        rna.obs['bc_batch'] = rna.obs[columns] #since it's one
        rna.obs["bc_batch"] = rna.obs["bc_batch"].astype("category")
    batch_categories = list(rna.obs['bc_batch'].unique())
    kwargs["batch_key"] = "bc_batch"
else:
    batch_categories = None
if args.integration_col_continuous is not None :
        rna.obs['bc_batch_continuous'] = rna.obs[args.integration_col_continuous]
        kwargs["continuous_covariate_keys"] = ["bc_batch_continuous"]


# add in raw counts as a layer 
# try to find or load the raw counts
if "counts" in rna.layers.keys():
    L.info("Found raw RNA counts in .layers['counts']")
elif "raw_counts" in rna.layers.keys():
    L.info("Found raw RNA counts in .layers['raw_counts']")
    L.info("Saving raw RNA counts to .layers['counts']")
    rna.layers["counts"] = rna.layers["raw_counts"].copy() 
elif X_is_raw(rna):
    # this means the X layer is already raw and we can make the layer we need
    L.info("Found raw RNA counts in .X. Saving raw RNA counts to .layers['counts']")
    rna.layers["counts"] = rna.X.copy()
else:
    L.info("Merging in raw RNA counts from '%s" % args.raw_anndata)
    sc_raw = read_anndata(args.raw_anndata, use_muon=True, modality="rna")
    #filter by barcodes in the scaled object
    sc_raw = sc_raw[sc_raw.obs_names.isin(rna.obs_names),: ]
    L.info("Saving raw RNA counts to .layers['counts]")
    rna.layers["counts"] = sc_raw.X.copy()

# filter prot outliers 
if params['multimodal']['totalvi']['filter_prot_outliers']:
    # for this to work the user needs to (manually) make a column called prot_outliers
    # actually there is a thing in the qc pipe that calculates outliers, I don't like it very much though
    if "prot_outliers" in mdata['prot'].columns:
        L.info("Filtering out prot outliers")
        mu.pp.filter_obs(mdata, "prot_outliers")
    else:
        raise ValueError("'prot_outliers' column not found in mdata['prot'].obs")

# exluding isotypes
if 'isotype' in prot.var.columns:
    L.info("Excluding isotypes")
    prot = prot[:, ~prot.var.isotype]

# filter out mitochondria
if params['multimodal']['totalvi']['exclude_mt_genes']:
    L.info("Filtering out mitochondrial genes")
    rna = rna[:, ~rna.var[params['multimodal']['totalvi']['mt_column']]]
    
# filter by Hvgs

if params['multimodal']['totalvi']['filter_by_hvg']:
    L.info("Filtering by HVGs")
    rna = rna[:, rna.var.highly_variable]

mdata.update()

#make protein obsm pandas
# need to find the raw counts in prot
if X_is_raw(prot):
    L.info("Found raw prot counts in .X")
    X_array = prot.X.copy()
elif "raw_counts" in prot.layers.keys(): 
    L.info("Found raw prot counts in .layers['raw_counts']")
    X_array = prot.layers['raw_counts'].copy()
else:
    raise AttributeError("Raw counts not found in prot, \
                        Store in either .X or in 'raw_counts' layer")


X_df = pd.DataFrame(X_array.todense(), index=prot.obs_names, columns=prot.var_names)

if X_df.shape[0] == rna.X.shape[0]:
    L.info("Adding protein_expression to .obsm['protein_expression']")
    #check the obs are in the correct order
    X_df = X_df.loc[rna.obs_names,:]
    X_df = X_df.astype('int')
    rna.obsm['protein_expression'] = X_df 
else:
    L.error("Dimension 0 of protein expression matrix and RNA expression matrix do not match")
    sys.exit("Dimension 0 of protein expression matrix and RNA expression matrix do not match")

# clear up to preserve ram
del X_array, X_df
gc.collect()

mdata.update()

L.info("Setting up AnnData")
scvi.model.TOTALVI.setup_anndata(
    rna,
    layer="counts",
    protein_expression_obsm_key="protein_expression",
    # categorical_covariate_keys = "bc_batch"
    **kwargs
)

if params['multimodal']['totalvi']['model_args'] is None:
    totalvi_model_args =  {}
else:
    totalvi_model_args =  {k: v for k, v in params['multimodal']['totalvi']['model_args'].items() if v is not None}


if params['multimodal']['totalvi']['training_args'] is None:
    totalvi_training_args = {}
else:
    totalvi_training_args =  {k: v for k, v in params['multimodal']['totalvi']['training_args'].items() if v is not None}


if params['multimodal']['totalvi']['training_plan'] is None:
    totalvi_training_plan = {}
else:
    totalvi_training_plan =  {k: v for k, v in params['multimodal']['totalvi']['training_plan'].items() if v is not None}

print(totalvi_model_args)
print(totalvi_training_args)
print(totalvi_training_plan)

L.info("Defining model")
vae = scvi.model.TOTALVI(rna, **totalvi_model_args)
L.info("Running totalVI")
vae.train(**totalvi_training_args, plan_kwargs=totalvi_training_plan)

vae.save(os.path.join("batch_correction", "totalvi_model"), 
                  anndata=False, overwrite=True )


L.info("Plotting ELBO")
plt.plot(vae.history["elbo_train"], label="train")
plt.plot(vae.history["elbo_validation"], label="validation")
plt.title("Negative ELBO over training epochs")
# plt.ylim(1200, 1400)
plt.legend()
plt.savefig(os.path.join(args.figdir, "totalvi_elbo_plot.png"))

# scvi.model..view_anndata_setup(vae.adata)

# we want to put things back in the original mdata

L.info("""We support the use of mudata as a general framework for multimodal data
        For this reason, the object we save is not the classical anndata
        you would find in scvitools tutorial.
        We chose to save the learned SC representation in the 
        Mudata obsm slot, and any single modality processing in
        its own modality slot
            """)

L.info("Extracting latent space and saving latent to X_totalVI")
mdata.obsm["X_totalVI"] = vae.get_latent_representation()

if batch_categories is not None:
    L.debug(batch_categories)
    if type(batch_categories) is not list:
        batch_categories = [batch_categories]
    normX, protein = vae.get_normalized_expression(
        n_samples=25,
        return_mean=True,
        transform_batch=batch_categories
    )
    L.info("Saving denoised data to .obsm['totalvi_denoised_rna'] and .obsm['totalvi_denoised_protein']")
    # this is stored in obsm, because if we have filtered by hvg, it no longer fits specifications for a layer
    mdata['rna'].obsm["totalvi_denoised_rna"] = normX
    # reorder if rna and prot are in different orders
    mdata['prot'].obsm["totalvi_denoised_protein"] = protein.loc[mdata['prot'].obs_names,:]
    #
    df = vae.get_protein_foreground_probability(
        n_samples=25,
        return_mean=True,
        transform_batch=batch_categories
    )
    L.info("Saving protein foreground prob to .obsm['totalvi_protein_foreground_prob']")
    mdata['prot'].obsm["totalvi_protein_foreground_prob"] = df.loc[mdata['prot'].obs_names,:]

mdata.update()

if int(args.neighbors_n_pcs) > mdata.obsm['X_totalVI'].shape[1]:
    L.warning(f"N PCs is larger than X_totalVI dimensions, reducing n PCs to  {mdata.obsm['X_totalVI'].shape[1] -1}")

n_pcs= min(int(args.neighbors_n_pcs), mdata.obsm['X_totalVI'].shape[1]-1)


# neighbors function has been made to magically calculate mudata neighbors using sc.pp.neighbors
# because we want to use_rep = X_totalvi (as a pseudo single modality) not run multimodal mu.pp.neighbors
L.info("Computing neighbors")
run_neighbors_method_choice(mdata, 
    method=args.neighbors_method, 
    n_neighbors=int(args.neighbors_k), 
    n_pcs=n_pcs, 
    metric=args.neighbors_metric, 
    use_rep='X_totalVI',
    nthreads=max([threads_available, 6]))

L.info("Computing UMAP")
sc.tl.umap(mdata, min_dist=0.4)
L.info("Computing Leiden clustering")
sc.tl.leiden(mdata, key_added="leiden_totalVI")

L.info("Saving UMAP coordinates to csv file '%s" % args.output_csv)
umap = pd.DataFrame(mdata.obsm['X_umap'], mdata['rna'].obs.index)
umap.to_csv(args.output_csv)

L.info("Saving MuData to 'tmp/totalvi_scaled_adata.h5mu'")
mdata.write("tmp/totalvi_scaled_adata.h5mu")

L.info("Done")

