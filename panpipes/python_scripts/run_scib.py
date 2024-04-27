import argparse
import os
import pandas as pd
from anndata import AnnData
from scib_metrics.benchmark import Benchmarker, BioConservation, BatchCorrection
from panpipes.funcs.io import read_yaml

import sys
import logging

L = logging.getLogger()
L.setLevel(logging.INFO)
log_handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter("%(asctime)s: %(levelname)s - %(message)s")
log_handler.setFormatter(formatter)
L.addHandler(log_handler)

# parse arguments
parser = argparse.ArgumentParser()

parser.add_argument('--combined_umaps_df')
parser.add_argument('--cell_meta_df')
parser.add_argument('--integration_dict')
parser.add_argument('--n_threads')
parser.add_argument("--fig_dir")
parser.add_argument("--rna_cell_type", default=None)
parser.add_argument("--prot_cell_type", default=None)
parser.add_argument("--atac_cell_type", default=None)

args = parser.parse_args()
L.info(args)

cell_meta_df = pd.read_csv(args.cell_meta_df, index_col=0)
batch_dict = read_yaml(args.integration_dict)
batch_dict = {modality: batch_col for modality, batch_col in batch_dict.items() if
              modality != "multimodal"}  # TODO: Check that multimodal should be included
umaps = pd.read_csv(args.combined_umaps_df, sep="\t", index_col=0)
cell_type_col = dict(rna=args.rna_cell_type, prot=args.prot_cell_type, atac=args.atac_cell_type)

for modality in batch_dict.keys():
    L.info("Computing scib metrics for modality: %s" % modality)
    # get one UMAP DataFrame per integration method for this modality
    modality_umaps = dict(list(umaps[umaps["mod"] == modality].groupby("method")))

    adata = AnnData(X=None)  # We only need the obs and obsm fields
    cell_ids = modality_umaps["none"].index
    adata.obs = cell_meta_df.loc[cell_ids, :]
    adata.obsm["Unintegrated"] = modality_umaps["none"].loc[:, ["umap_1", "umap_2"]].to_numpy()

    for method, umap_df in modality_umaps.items():
        adata.obsm[method] = umap_df.loc[cell_ids, ["umap_1", "umap_2"]].to_numpy()

    # Check if cell type information is available
    if cell_type_col[modality] is None:
        batch_correction_metrics = BatchCorrection(silhouette_batch=False,
                                                   ilisi_knn=True,
                                                   kbet_per_label=False,
                                                   graph_connectivity=False,
                                                   pcr_comparison=True)
        bio_conservation_metrics = BioConservation(isolated_labels=False,
                                                   nmi_ari_cluster_labels_leiden=False,
                                                   nmi_ari_cluster_labels_kmeans=False,
                                                   silhouette_label=False,
                                                   clisi_knn=False)

        # Add a dummy cell type column
        adata.obs["cell_type"] = "None"
        cell_type_col[modality] = "cell_type"
    else:
        batch_correction_metrics = None
        bio_conservation_metrics = None

    bm = Benchmarker(
        adata,
        batch_key=batch_dict[modality],
        label_key=cell_type_col[modality],
        embedding_obsm_keys=["Unintegrated"].extend([method for method in modality_umaps.keys() if method != "none"]),
        pre_integrated_embedding_obsm_key="Unintegrated",
        batch_correction_metrics=batch_correction_metrics,
        bio_conservation_metrics=bio_conservation_metrics,
        n_jobs=args.n_threads,
    )
    bm.benchmark()

    bm.plot_results_table(show=False, save_dir=os.path.join(args.fig_dir, modality))
    bm.plot_results_table(min_max_scale=False, show=False, save_dir=os.path.join(args.fig_dir, modality))
    df = bm.get_results(min_max_scale=False)
