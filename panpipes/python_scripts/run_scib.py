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
    # get one df per method for this modality
    splits = dict(list(umaps[umaps["mod"]==modality].groupby("method")))
    columns = batch_dict[modality]
    for k, xx in splits.items():
        L.info("computing lisi for correction %s" % k)
        # compute LISI for each batch
        umap_coords = xx.loc[:,["umap_1", "umap_2"]].to_numpy()
        # check it"s in the correct order
        batch_df = cell_meta_df.loc[xx.index,:]

bm = Benchmarker(
    adata,
    batch_key="batch",
    label_key="cell_type",
    embedding_obsm_keys=["Unintegrated", "Scanorama", "LIGER", "Harmony", "scVI", "scANVI"],
    n_jobs=6,
)
bm.benchmark()


bm.plot_results_table()
bm.plot_results_table(min_max_scale=False)
df = bm.get_results(min_max_scale=False)