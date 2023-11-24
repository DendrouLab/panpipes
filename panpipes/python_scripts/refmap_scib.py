import muon as mu
import scib
import pandas as pd
import numpy as np
import seaborn as sns
import scanpy as sc
sc.set_figure_params(figsize=(8, 6), dpi=300)
import argparse
import gc
import os
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
parser.add_argument('--query',
                    default='',
                    help='the query anndata/mudata')
parser.add_argument('--repuse',
                    default=None,
                    help='which latent representation to use to calculate the metrics')
parser.add_argument('--covariate',
                    default=None,
                    help='the convariate in the query anndata/mudata you want to use to calculate the metrics')
parser.add_argument('--batch_key',
                    default=None,
                    help='the convariate in the query anndata/mudata you want to use to calculate the metrics')
parser.add_argument('--cluster_key',
                    default=None,
                    help='string of column in query anndata/mudata containing cluster assignments you want to use to calculate the metrics (ARI/NMI)')
parser.add_argument('--outdir', default=None,
                    help='directory to save the data to')

args, opt = parser.parse_known_args()

L.info(args)


if args.covariate is not None:
    covariates_use = args.covariate.split(",")
    covariates_use = [a.strip() for a in covariates_use]
else:
    sys.exit("i don't have covariates to calculate metrics on")


repuse = str(args.repuse)

mdata = mu.read(args.query)
if type(mdata) is mu.MuData:
    if "rna" not in mdata.mod.keys():
        sys.exit("we only support querying using RNA but your mdata doesn't contain rna")
    else:
        input_adata = mdata["rna"].copy()
del mdata


adata_query = input_adata[input_adata.obs['is_reference'] == 'Query'].copy()
L.info("repuse is %s" %(repuse))
if repuse not in adata_query.obsm.keys():
    sys.exit("the latent representation is not in the obsm.keys of this query")

L.info("query is:")
print(adata_query)

L.info("Calculating scib metrics using ground truth covariates:")
print(covariates_use)

if args.cluster_key is None:
    cluster_key = "leiden_" + repuse 
    L.info("you didn't specify cluster_key, so i'm using %s " % cluster_key)
else:
    cluster_key = str(args.cluster_key)
    L.info("cluster_key used is: %s" %cluster_key)

results = {}
for labelk in covariates_use:
    m={"ASW_scaled":scib.metrics.silhouette(adata_query, labelk, repuse, metric='euclidean', scale=True),
        "ARI":scib.metrics.ari(adata_query, label_key=labelk, cluster_key=cluster_key), #predictions or leiden clustering
        "NMI":scib.metrics.nmi(adata_query, label_key=labelk, cluster_key=cluster_key),
        "graph_connectivity":scib.metrics.graph_connectivity(adata_query, labelk)}    
    if args.batch_key is not None:
        m["clisi_graph_embed"]= scib.me.clisi_graph(adata_query, label_key=labelk, type_="knn", batch_key=str(args.batch_key))

    results[labelk] = m

file_name= os.path.splitext(os.path.basename(args.query).replace("query_to_reference_", "").replace(".h5mu", ""))[0]
#"query_to_reference_" + model_name + "_" + latent_choice + ".h5mu"
L.info("saving file output")
file_out = os.path.join(args.outdir,("scib.query_"+ file_name+".csv"))
print(file_out)
#save file to txt
pres = pd.DataFrame.from_dict(results,orient='index')
pres.to_csv(file_out, sep="\t")
L.info("Finished scib")