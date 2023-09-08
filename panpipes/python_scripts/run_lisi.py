import argparse
import numpy as np
import pandas as pd
import harmonypy as hm
import seaborn as sns
import os
from panpipes.funcs.io import read_yaml
import matplotlib.pyplot as plt

import sys
import logging
L = logging.getLogger()
L.setLevel(logging.INFO)
log_handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter("%(asctime)s: %(levelname)s - %(message)s")
log_handler.setFormatter(formatter)
L.addHandler(log_handler)

parser = argparse.ArgumentParser()
parser.add_argument("--combined_umaps_df")
parser.add_argument("--cell_meta_df")
parser.add_argument("--integration_dict")
parser.add_argument("--fig_dir")
args = parser.parse_args()
L.info(args)

# load combined umaps and cellmtd file
# load metadata
cell_meta_df = pd.read_csv(args.cell_meta_df, index_col=0)
batch_dict = read_yaml(args.integration_dict)
umaps = pd.read_csv(args.combined_umaps_df, sep="\t", index_col=0)
L.info("fixing the batch_keys dictionary")
for k in batch_dict.keys():
    v=batch_dict[k]
    if k =="multimodal":
        mk = []
        exc = []
        for v2 in v:
            if v2 in cell_meta_df.columns:
                mk.append(v2)
            else:
                exc.append(v2)
                for k in batch_dict.keys():
                    if k!="multimodal":
                       tmp=[k + ":" +e for e in exc] #add the prepending modality
                       for t in tmp:
                           if t in cell_meta_df.columns:
                               mk.append(t)
        v = mk
        batch_dict[k] = v
   # else:
   #     v =[k + ":" +v1 if k!="multimodal" else v1 for v1 in v  ] #add the prepending modality
   #     batch_dict[k] = v
    if len(v) > 1:
        mod_meta_df= cell_meta_df[v].dropna()
        cell_meta_df[str(k)+ ":bc_batch"] = mod_meta_df.apply(lambda x: "|".join(str(x)), axis=1)
        batch_dict[k].append(k+ ":bc_batch")
    L.info("batch keys %s" %(batch_dict[k]))


for md in batch_dict.keys():
    L.info("Running lisi on modality: %s" % md)
    # get one df per method for this modality
    splits = dict(list(umaps[umaps["mod"]==md].groupby("method")))
    lisi_results = []
    columns = batch_dict[md]
    for k, xx in splits.items():
        L.info("computing lisi for correction %s" % k)
        # compute LISI for each batch
        umap_coords = xx.loc[:,["umap_1", "umap_2"]].to_numpy()
        # check it"s in the correct order
        batch_df = cell_meta_df.loc[xx.index,:]
        res = hm.compute_lisi(umap_coords, batch_df, columns)
        lisi_results.append(pd.DataFrame(res, index=batch_df.index, columns=columns))
    # put LISI scores into a pandas df
    lisi_df = pd.concat(lisi_results, keys=splits.keys()).reset_index().rename(columns={"level_0":"method", "level_1":"cellbarcode", "index":"cellbarcode"})
    lisi_df.to_csv(os.path.join(args.fig_dir, md, "LISI_scores.csv"), index=False)
    # # make a density  plot of LISI scores.
    L.info("plotting lisi density")
    plot_df = lisi_df.melt(id_vars=['cellbarcode', 'method'], var_name="integration_variable", value_name="LISI score")
    plot_df["integration_variable"] = plot_df["integration_variable"].astype("category")
    plot_df = plot_df.rename(columns={"method": "Correction"})
    plot_df["Correction"] = plot_df["Correction"].astype('category')
    # # reorder categories to match other plots
    methods = plot_df["Correction"].unique().tolist()
    if "none" in methods:
        methods.insert(0, methods.pop(methods.index("none")))
        plot_df["Correction"] = plot_df["Correction"].cat.reorder_categories(methods)
    # do plot
    fig, ax = plt.subplots(nrows=1, ncols=len(columns), figsize=(len(columns)*4, 4))
    if len(columns) > 1:
        for idx, col in enumerate(columns):
            sns.kdeplot(data=plot_df[plot_df.integration_variable==col],
                x="LISI score", hue="Correction", ax=ax[idx])
            ax[idx].set_title("integrated by :" + col) 
    else:
        sns.kdeplot(data=plot_df,
            x="LISI score", hue="Correction", ax=ax)
        ax.set_title("integrated by :" + columns[0])  
    fig.savefig(os.path.join(args.fig_dir, md, "LISI_scores.png"))
    plt.clf()



L.info("Done")