import argparse
import pandas as pd
from anndata import AnnData
from scib_metrics.benchmark import Benchmarker
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

parser.add_argument('--integrated_anndata')
parser.add_argument("--fig_dir")

args = parser.parse_args()
L.info(args)


cell_meta_df = pd.read_csv(args.cell_meta_df, index_col=0)
batch_dict = read_yaml(args.integration_dict)
batch_dict = {modality: batch_col for modality, batch_col in batch_dict.items() if modality!="multimodal"} #TODO: Check that multimodal should be included
umaps = pd.read_csv(args.combined_umaps_df, sep="\t", index_col=0)


for modality in batch_dict.keys():
    L.info("Computing scib metrics for modality: %s" % modality)
    # get one UMAP DataFrame per integration method for this modality
    modality_umaps = dict(list(umaps[umaps["mod"] == modality].groupby("method")))

    adata = AnnData(X=None)
    cell_ids = modality_umaps["none"].index
    adata.obs = cell_meta_df.loc[cell_ids, :]
    adata.obsm["Unintegrated"] = modality_umaps["none"].loc[:, ["umap_1", "umap_2"]].to_numpy()

    for method, umap_df in modality_umaps.items():
        adata.obsm[method] = umap_df.loc[cell_ids, ["umap_1", "umap_2"]].to_numpy()


    bm = Benchmarker(
        adata,
        batch_key=batch_dict[modality],
        label_key="cell_type",
        embedding_obsm_keys=["Unintegrated", "Scanorama", "LIGER", "Harmony", "scVI", "scANVI"],
        n_jobs=6,
    )
    bm.benchmark()


bm.plot_results_table()
bm.plot_results_table(min_max_scale=False)
df = bm.get_results(min_max_scale=False)