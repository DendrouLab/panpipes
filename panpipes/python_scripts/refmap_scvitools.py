import multiprocessing
from matplotlib import transforms 

import anndata as ad
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import scvi
import argparse
import os
import muon as mu
from mudata import MuData
import numpy as np
import panpipes.funcs as pp
from panpipes.funcs.processing import check_for_bool
from panpipes.funcs.io import read_yaml
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
parser.add_argument('--query_data',
                    default='adata_scaled.h5ad',
                    help='path to query data. can be a raw 10x dataset or a preprocessed anndata/mudata')
parser.add_argument('--query_celltype',default='celltype',
                    help='if query has already a column with cell labeling in obs. Default to "celltype" ')
parser.add_argument('--adata_reference',
                    default='adata_.h5ad',
                    help='path to reference ann/mudata. necessary only if you want to produce ref/ query map')
parser.add_argument('--reference_path',
                    default='reference model directory',
                    help='')
parser.add_argument('--reference_architecture', default='scanvi',
                    help='scvi, scanvi, totalvi supported')
# parser.add_argument('--generate_scanvi', default=False,
#                     help='create scanvi reference from scvi to allow for label transfer')
parser.add_argument('--predict_rf',default=False,
                    help='For totalvi and scvi, if reference model has random forest classifier, classify cells in query')
parser.add_argument('--impute_proteins', default=False,
                    help='if totalvi reference is used, allow to impute protein levels')
parser.add_argument('--transform_batch', default=None,
                    help='if totalvi reference is used, batch categories for the query dataset')
parser.add_argument('--neighbors_n_pcs',
                    help="n_pcs", default=50)
parser.add_argument('--neighbors_method', default="scanpy",
                    help="neighbours method, scanpy or hnsw")
parser.add_argument('--neighbors_k', default=30,
                    help="neighbors k")
parser.add_argument('--neighbors_metric',default="euclidean",
                    help="neighbor metric, e.g. euclidean or cosine")


args, opt = parser.parse_known_args()
sc.settings.figdir = "figures/"
sc.set_figure_params(figsize=(8, 6), dpi=300) 
L.info("running with args:")
args.predict_rf = check_for_bool(args.predict_rf)
args.impute_proteins = check_for_bool(args.impute_proteins)
L.info(args)

threads_available = multiprocessing.cpu_count()

# load parameters
params = read_yaml("pipeline.yml")
query_data = os.path.basename(args.query_data)

reference_architecture = str(args.reference_architecture)



mdata = mu.read(args.query_data)
if type(mdata) is mu.MuData:
    if "rna" not in mdata.mod.keys():
        sys.exit("we only support querying using RNA but your mdata doesn't contain rna")
    else:
        adata_query = mdata["rna"].copy()
        if "prot" in mdata.mod.keys() and reference_architecture=="totalvi":
            if X_is_raw(mdata['prot']):
                X_array = mdata['prot'].X.copy()
            elif "raw_counts" in mdata['prot'].layers.keys(): 
                X_array = mdata['prot'].layers['raw_counts'].copy()
        
        X_df = pd.DataFrame(X_array.todense(), index=mdata['prot'].obs_names, columns=mdata['prot'].var_names)
        if X_df.shape[0] == adata_query.X.shape[0]:
            L.info("adding protein_expression to obsm")
            #check the obs are in the correct order
            X_df = X_df.loc[adata_query.obs_names,:]
            adata_query.obsm['protein_expression'] = X_df 
        else:
            L.error("dimensions do not match, cannot create the raw obsm counts in query")
else:
    adata_query = mdata.copy()

L.info("this is your query anndata")
print(adata_query)

if "counts" not in adata_query.layers.keys():
    if "raw_counts" in adata_query.layers.keys():
        adata_query.layers["counts"] = adata_query.layers["raw_counts"].copy()
    elif X_is_raw(adata_query):
        adata_query.layers["counts"] = adata_query.X.copy()
    else:
        raise AttributeError("raw counts not found in query, \
                        store in either X or in 'counts/raw_counts' layer")


if reference_architecture in ['totalvi','scvi','scanvi']:
    L.info("running using %s" % reference_architecture)
else:
    sys.exit(" i don't recognise this architecture : %s") %reference_architecture

reference_path = args.reference_path
if not os.path.exists(reference_path):
    L.info("The reference path you provided does not exist")
    sys.exit("The reference path you provided does not exist")
else:
    reference_path = os.path.dirname(reference_path)
    L.debug("Reference path: %s" % reference_path)

train_kwargs = {}
if params["training_plan"][reference_architecture] is not None:
    train_kwargs = params["training_plan"][reference_architecture]
    if train_kwargs["max_epochs"] is not None:
        max_epochs = train_kwargs["max_epochs"]
        del train_kwargs["max_epochs"]
    else:
        max_epochs = 200
else:
    max_epochs = 200
    train_kwargs = {'weight_decay': 0.0}


if reference_architecture=="scvi":
    scvi.model.SCVI.prepare_query_anndata(adata_query, reference_path)
    vae_q = scvi.model.SCVI.load_query_data(
    adata_query,
    reference_path)
    latent_choice= "X_scVI"
    vae_q.train(max_epochs= max_epochs , plan_kwargs=train_kwargs)
    adata_query.obsm["X_scVI"] = vae_q.get_latent_representation()
    

if reference_architecture=="scanvi":
    # Notice that adata_query.obs["labels_scanvi"] does not exist. 
    # The load_query_data method detects this and fills it in adata_query with the unlabeled category (here "Unknown").
    scvi.model.SCANVI.prepare_query_anndata(adata_query, reference_path)
    vae_q = scvi.model.SCANVI.load_query_data( 
    adata_query,
    reference_path)
    latent_choice= "X_scANVI"
    #vae_q.train(**train_kwargs) this doesn't work anymore cause max_epochs is not recognised as part of the plan kwargs
    vae_q.train(max_epochs= max_epochs , plan_kwargs=train_kwargs)
    adata_query.obsm["X_scANVI"] = vae_q.get_latent_representation()
    adata_query.obs["predictions"] = vae_q.predict()
    if args.query_celltype is not None:
        L.info("Query has celltypes in column %i, i will plot what predictions look like from scanvi model" % args.query_celltype)
        df = adata_query.obs.groupby([str(args.query_celltype), "predictions"]).size().unstack(fill_value=0)
        norm_df = df / df.sum(axis=0)

        plt.figure(figsize=(8, 8))
        _ = plt.pcolor(norm_df)
        _ = plt.xticks(np.arange(0.5, len(df.columns), 1), df.columns, rotation=90)
        _ = plt.yticks(np.arange(0.5, len(df.index), 1), df.index)
        plt.xlabel("Predicted")
        plt.ylabel("Observed")
        file_name = "SCANVI_predicted_vs_observed_labels_query_data"
        plt.savefig(os.path.join("figures/", file_name + ".png"))
        
if reference_architecture=="totalvi":
    # temporary fix is disabled for now, need to modify to allow to fix query to match reference
    fix_query = False
    if fix_query:
        L.info(" will do some manipulation of query data to make it match to referece structure")
        if args.adata_reference is not None:
            reference_data = os.path.basename(args.adata_reference)
            mdata=mu.read(args.adata_reference)
            if type(mdata) is MuData:
                if "rna" not in mdata.mod.keys():
                    sys.exit("we only support querying using RNA but your mdata doesn't contain rna")
                else:
                    adata_ref = mdata["rna"].copy()
                    adata_ref.obsm = mdata.obsm.copy()
                    adata_ref.obsp = mdata.obsp.copy()
            else:
                adata_ref = mdata.copy()
            del mdata
            # match query and reference vars
            
            #adata_query.layers["counts"] = adata_query.X.copy()#already taken care of 
            # are the following necessary?
            sc.pp.normalize_total(adata_query, target_sum=1e4)
            sc.pp.log1p(adata_query)
            adata_query.raw = adata_query
            # subset to reference vars
            adata_query = adata_query[:, adata_ref.var_names].copy()
            
            adata_query.obs["celltype.l2"] = "Unknown"
            adata_query.obs["orig.ident"] = adata_query.obs["set"]
            #query.obsm["X_umap"] = query.obs[["UMAP1", "UMAP2"]]


            # obsm slot slot for protein counts can be called either prot_counts or protein_expression
            # in the reference anndata and model. allow a 

            if "protein_expression" in adata_ref.obsm.keys():
                pname = "protein_expression"
            
            elif "protein_counts" in adata_ref.obsm.keys():
                pname = "protein_counts"
            
            if pname not in adata_query.obsm.keys():
                X_df = pd.DataFrame(0, index=adata_ref.obs_names, columns=adata_ref.obsm[pname].columns)

                adata_query.obsm["protein_counts"] = query.obsm["pro_exp"].copy()
            

            for p in adata_ref.obsm[pname].columns:
                if p not in adata_query.obsm[pname].columns:
                    adata_query.obsm[pname][p] = 0.0
            # ensure columns are in same order
            adata_query.obsm[pname] = adata_query.obsm[pname].loc[:, adata_ref.obsm[pname].columns]

        else:
            
            # For now make this only work when reference anndata is available

            # if args.reference_var is not None:
            #     reference_var = pd.read(args.reference_var, sep="\t")
            # else: 
            #     sys.exit("I need var names from reference to make sure query has overlap")    
            # if args.reference_prot is not None:
            #     reference_prot = pd.read(args.reference_prot, sep="\t")
            # else: 
            #     sys.exit("I need var names from reference to make sure query has overlap")    

            # if args.reference_prot_assay is not None:
            #     L.info("name of obsm slot for reference is: %i" % args.reference_prot_assay)
            #     pname = args.reference_prot_assay
            sys.exit("need a reference dataset to check for matching query entries")        
        
    scvi.model.TOTALVI.prepare_query_anndata(adata_query, reference_path)
    vae_q = scvi.model.TOTALVI.load_query_data(
    adata_query,
    reference_path,
    freeze_expression=True)
    latent_choice= "X_totalvi"
    vae_q.train(max_epochs= max_epochs , plan_kwargs=train_kwargs)
    adata_query.obsm["X_totalvi"] = vae_q.get_latent_representation()
    #remove this after finishing
    adata_query.write(os.path.join( "query_temp_check_totvi.h5ad"))
        
    if args.predict_rf :
        L.info("predicting celltypes")
        predictions = (
           vae_q.latent_space_classifer_.predict(
            adata_query.obsm["X_totalvi"]
            )
        )
        adata_query.obs["predictions"] = predictions
        #remove this after finishing
        adata_query.write(os.path.join( "query_temp_check_predictions_totvi.h5ad"))
    
        if args.query_celltype is not None:
            L.info("Query has celltypes in column %s, i will plot what predictions look like from totalvi model" % args.query_celltype)
            df = adata_query.obs.groupby([str(args.query_celltype), "predictions"]).size().unstack(fill_value=0)
            norm_df = df / df.sum(axis=0)

            plt.figure(figsize=(8, 8))
            _ = plt.pcolor(norm_df)
            _ = plt.xticks(np.arange(0.5, len(df.columns), 1), df.columns, rotation=90)
            _ = plt.yticks(np.arange(0.5, len(df.index), 1), df.index)
            plt.xlabel("Predicted")
            plt.ylabel("Observed")
            file_name = "totalvi_predicted_vs_observed_labels_query_data"
            plt.savefig(os.path.join("figures/", file_name + ".png"))
        
    
L.info("Mapped Query to Reference, now some utilities")


if args.adata_reference is not None:
    reference_data = os.path.basename(args.adata_reference)
    mdata=mu.read(args.adata_reference)
    if type(mdata) is MuData:
        if "rna" not in mdata.mod.keys():
            sys.exit("we only support querying using RNA but your mdata doesn't contain rna")
        else:
            adata_ref = mdata["rna"].copy()
            adata_ref.obsm = mdata.obsm.copy()
            adata_ref.obsp = mdata.obsp.copy()
    else:
        adata_ref = mdata.copy()

    adata_ref.obs.loc[:, 'is_reference'] = 'Reference'
    adata_query.obs.loc[:, 'is_reference'] = 'Query'
    #expect the batch to be always encoded as `batch` in both Q and R
    adata_full = ad.concat( [adata_ref,adata_query])
    # add param to decide if to recompute the total embedding or to recalc the embedding
    #if "X_totalvi" not in adata_ref.obsm.keys():
    adata_full.obsm[latent_choice] = vae_q.get_latent_representation(adata_full)
    
    if args.impute_proteins: 
        if args.transform_batch is not None:
            transform_batch = str(args.transform_batch)
        else:
            transform_batch = None
        normX, protein = vae_q.get_normalized_expression(
            adata_full,
            n_samples=25,
            return_mean=True,
            transform_batch= transform_batch) 
        adata_full.obsm["totalvi_denoised_rna"], adata_full.obsm["totalvi_denoised_protein"] = normX, protein    

    if reference_architecture=="scanvi":
        full_predictions = vae_q.predict(adata_full)
        L.warn("Acc: {}".format(np.mean(full_predictions == adata_full.obs.celltype)))
        adata_full.obs["predictions"] = full_predictions
else:
    adata_query.obs.loc[:, 'is_reference'] = 'Query'
    adata_full = adata_query.copy()

if int(args.neighbors_n_pcs) > adata_full.obsm[latent_choice].shape[1]:
    L.warn("N PCs is larger than %i dimensions, reducing n PCs to " % adata_full.obsm[latent_choice].shape[1])


n_pcs= min(int(args.neighbors_n_pcs), adata_full.obsm[latent_choice].shape[1])


run_neighbors_method_choice(adata_full, 
            method=args.neighbors_method, 
            n_neighbors=int(args.neighbors_k), 
            n_pcs=n_pcs, 
            metric=args.neighbors_metric, 
            use_rep=latent_choice,
            nthreads=max([threads_available, 6]))
sc.tl.umap(adata_full, min_dist=0.4)
sc.tl.leiden(adata_full, key_added="leiden_" + latent_choice)

model_name = os.path.basename(os.path.dirname(args.reference_path))
model_name = ''.join(e for e in model_name if e.isalnum())


file_name= "umap_" + model_name + "_" + latent_choice

fig = sc.pl.embedding(adata_full, basis = "umap",color=["is_reference"],
         show=False, return_fig=True)

fig.savefig(os.path.join("figures/", file_name + ".png"))

umap = pd.DataFrame(adata_full.obsm['X_umap'], adata_full.obs.index)

umap.to_csv(os.path.join("refmap/", file_name + ".csv") )
file_name= "query_to_reference_" + model_name + "_" + latent_choice + ".h5mu" #change this to mudata

mdata_save = MuData({"rna":adata_full})

mdata_save.write(os.path.join( "refmap/" , file_name))

L.info('done')



