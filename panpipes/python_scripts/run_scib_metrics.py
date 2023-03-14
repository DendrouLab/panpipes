import scanpy as sc
import scIB as scib
import glob
import pandas as pd
import argparse
from panpipes.funcs.processing import check_for_bool
from panpipes.funcs.io import read_anndata
from panpipes.funcs.scmethods import merge_consensus_clust

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

parser.add_argument('--uncorrected_anndata',
                    default='adata_scaled.h5ad',
                    help='')
parser.add_argument('--batch_corrected_anndata', default='tmp/harmony_scaled.h5ad',
                    help='')
parser.add_argument('--use_muon',
                    default=False,
                    help='') 
parser.add_argument('--integration_method', default='batch',
                    help='')
parser.add_argument('--integration_col', default='batch',
                    help='')
parser.add_argument('--rough_consensus_clustering', default=None,
                    help='must match the barcodes in the anndata objects')
parser.add_argument('--reference_col', default=None,
                    help='')
parser.add_argument('--outfile', default='scib.tsv',
                    help='')
args = parser.parse_args()
args, opt = parser.parse_known_args()
use_muon = check_for_bool(args.use_muon)

L.info("reading data and starting integration pipeline with script: ")
L.info("Running with options: %s", args)

# which metrics to test?
batch_metrics = {
    "isolated_labels_asw_": False,
    "hvg_score_": True,
    "pcr_": True,
    "cell_cycle_": True,
    "ilisi_": True,

}
bio_metrics = {
    "clisi_": True,
    "nmi_": True,
    "ari_": True,
    "isolated_labels_f1_": True,
    "silhouette_": True,
    "graph_conn_": True,
    "kBET_": True,
}
adata = read_anndata(args.uncorrected_anndata, use_muon=use_muon, modality="rna")
# this is just an anndata object
adata_int = read_anndata(args.batch_corrected_anndata, use_muon=False, modality="rna")


# get rough consensus clustering
if args.reference_col is not None:
    if args.rough_consensus_clustering is not None:
        consensus_clust = pd.read_csv(args.rough_consensus_clustering, index_col=0)
        # merge in the consensus clustering
        adata = merge_consensus_clust(adata, consensus_clust)
        adata = adata[adata.obs[args.reference_col].notnull(),:]
        # merge in the consensus clustering
        adata_int = merge_consensus_clust(adata_int, consensus_clust)
        adata_int = adata_int[adata_int.obs[args.reference_col].notnull(),:]
    else:
        # check that reference col is found in the adata
        if args.reference_col not in adata.obs.columns:
            sys.exit("reference column not found in adata object")
    metrics_kwargs = {**batch_metrics, **bio_metrics}
    key_kwargs={"batch_key": args.integration_col, 
                "label_key": args.reference_col}  
else:
    metrics_kwargs = batch_metrics
    key_kwargs={"batch_key": args.integration_col,
            "label_key": "sample_id"}   

if args.integration_method == "harmony":
    x_embed = "X_harmony"
elif args.integration_method == "scanorama":
    x_embed = "X_scanorama"
elif args.integration_method == "scvi":
    x_embed = "X_scVI"
else:
    x_embed = "X_pca"

 
L.info("running scib on")
L.info(metrics_kwargs)
metrics = scib.metrics.metrics(adata=adata,
                        adata_int=adata_int, 
                        cluster_key="leiden",
                         embed=x_embed, 
                         organism="human", **key_kwargs, **metrics_kwargs)

metrics.to_csv(args.outfile)