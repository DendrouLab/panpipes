import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import muon as mu
import os
import re
from matplotlib import cm
from matplotlib.colors import ListedColormap
import yaml 
from panpipes.funcs.io import read_yaml
from panpipes.funcs.plotting import batch_scatter_two_var
from pandas.api.types import is_string_dtype

import sys
import logging
L = logging.getLogger()
L.setLevel(logging.INFO)
log_handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s: %(levelname)s - %(message)s')
log_handler.setFormatter(formatter)
L.addHandler(log_handler)


parser = argparse.ArgumentParser()
parser.add_argument("--cell_meta_df")
parser.add_argument("--combined_umaps_tsv")
parser.add_argument("--fig_dir")
parser.add_argument("--batch_dict", default=None)
parser.add_argument("--qc_dict", default=None)

args, opt = parser.parse_known_args()


L.info(args)
palette_choice = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
# cmap_choice = "BuPu" # previous default = "viridis"
bupu = cm.get_cmap('BuPu', 512)
cmap_choice= ListedColormap(bupu(np.linspace(0.1, 1, 256)))


# load metadata
cell_meta_df = pd.read_csv(args.cell_meta_df, index_col=0)
umaps_df = pd.read_csv(args.combined_umaps_tsv, sep='\t', index_col=0)
umaps_df = pd.merge(umaps_df, cell_meta_df, left_index=True, right_index=True)

batch_dict = read_yaml(args.batch_dict)

# parse the qc dict
qc_dict = args.qc_dict
if isinstance(args.qc_dict, dict):
    qc_dict = args.qc_dict
else:
    qc_dict = read_yaml(args.qc_dict) 

# extract the grouping vars which are expected to be categorical
grpvar = [x.replace (" " ,"") for x in qc_dict['grouping_var'].split(",")]
grpvar = [value for value in grpvar if value in cell_meta_df.columns] 
del qc_dict['grouping_var']
L.info("additional grouping vars:")
# merge with the batch columns 
from itertools import chain
all_grpvar = list(batch_dict.values()) + [grpvar]
# unnest the list
grpvar = list(set(chain(*all_grpvar)))
L.info(grpvar)

# modality level metrics
# remove metrics from keys
qc_dict = {re.sub("_metrics", "", k) : v for k, v in qc_dict.items()}

# split the comma sep string of metrics (and remove Nones)
qc_dict = {k: v.split(',') for k, v in qc_dict.items() if v is not None}
# remove white spaces from list of the dictionary
qc_dict ={i: [a.strip(" ") for a in j ] for i,j in qc_dict.items()}
# make sure mod: prefaces each metric
qc_dict = {k: [k + ":" + vi if ":" not in vi else vi for vi in v ]  for k, v in qc_dict.items()}
L.info("additional qc vars:")
L.info(qc_dict)


# do plots
for mod in umaps_df['mod'].unique():
    # create directory if it doesn't exist
    if os.path.exists(os.path.join(args.fig_dir, mod)) is False:
        os.mkdir(os.path.join(args.fig_dir, mod))
    L.info("plotting modality: %s" % mod)
    plt_df = umaps_df[umaps_df['mod'] == mod].copy()
    pointsize = 150000 / plt_df.shape[0]
    plt_df["method"] = plt_df["method"].astype("category")
    # put none at the top of the list
    if mod != "multimodal":
        umaps_methods_list = plt_df['method'].cat.categories.tolist()
        umaps_methods_list.insert(0, umaps_methods_list.pop(umaps_methods_list.index("none")))
        plt_df['method'] = plt_df['method'].cat.reorder_categories(umaps_methods_list, ordered=True)
    #  get grouing vars
    if mod !="multimodal":
        if mod in batch_dict.keys() :
            columns = set(batch_dict[mod] + [mod +":"+ g for g in grpvar])
        else:
            # this means no batch correction was run. 
            columns = set([mod +":"+ g for g in grpvar])
    else:
        columns = grpvar
    # keep only the columns that are actually in plt_df
    columns = [cl for cl in columns if cl in plt_df.columns ]
    for col in columns:
        L.info("plotting grouping var %s" % col)
        # make sure everything is a category so that when they get plotted it's the same color
        plt_df[col] = plt_df[col].astype(str).astype('category')
        plt_df = plt_df.sort_values(by="umap_1")
        g = sns.FacetGrid(plt_df, col="method", col_wrap=3, sharex=False, sharey=False)
        g = (g.map(sns.scatterplot, "umap_1", "umap_2", col, s=pointsize, linewidth=0))
        g.add_legend() 
        g.savefig(os.path.join(args.fig_dir, mod, "umap_method_" + str(col) + ".png"))
        fig, ax = batch_scatter_two_var(plt_df, "method", col, palette_choice=palette_choice)
        if fig is not None:
            fig.savefig(os.path.join(args.fig_dir, mod,  "umap_method_facet_" + str(col) + ".png"), dpi=300)
        plt.clf()

    ncats  = len(plt_df['method'].unique())
    if ncats > 4:
        ncols=4
        nrows = int(np.ceil(ncats/4))
    else:
        ncols=int(ncats)
        nrows=1

    qcmetrics = []
    if mod in qc_dict.keys() and qc_dict[mod] is not None:
        qcmetrics =  qcmetrics + qc_dict[mod] 
    if "all" in qc_dict.keys() and qc_dict["all"] is not None:    
        qcmetrics = qcmetrics + qc_dict['all'] 
    if len(qcmetrics) == 0:
        next
    # but we'll still check
    # keep only the qc metrics that are actually in plt_df
    else:
        qcmetrics = [qc for qc in qcmetrics if qc in plt_df.columns ]
        if len(qcmetrics) > 0:
            for qc in qcmetrics:
                fig, axes = plt.subplots(ncols=ncols, nrows=nrows, figsize=(5*ncols,4*nrows))
                if nrows > 1 or ncols > 1:
                    axes = axes.ravel()
                else:
                    axes=[axes]
                if pd.api.types.is_numeric_dtype(plt_df[qc]):
                    L.info("plotting qc var (numeric) %s" % qc)
                    for idx, mm in enumerate(plt_df['method'].cat.categories):
                        im = axes[idx].scatter(data=plt_df[plt_df['method']== mm], 
                                        x="umap_1",
                                        y="umap_2",
                                        c=qc, 
                                        s=pointsize,
                                        cmap=cmap_choice)
                        axes[idx].set_title(mm)
                        axes[idx].axis('off')
                    # delte the last few empty plots.
                    for idx2 in range(idx,ncols*nrows):
                        axes[idx2].axis('off')
                    if pd.api.types.is_numeric_dtype(plt_df[qc]):
                        fig.subplots_adjust(top=0.925,right=0.8)
                        cbar_ax = fig.add_axes([0.85, 0.35, 0.025, 0.35])
                        fig.colorbar(im, cbar_ax)
                    fig.suptitle(qc)
                    plt.savefig(os.path.join(args.fig_dir, mod , "umap_method_" + qc + ".png"), dpi = 300)
                    plt.clf()
                else:
                    # this is a categorical colored plot
                    plt_df[qc] = plt_df[qc].fillna('nan')
                    if len(plt_df[qc].unique()) < 40:
                        L.info("plotting qc var (categorical) %s"%qc)
                        g = sns.FacetGrid(plt_df, col="method", col_wrap=3, sharex=False, sharey=False)
                        g = (g.map(sns.scatterplot, "umap_1", "umap_2", qc, s=pointsize, linewidth=0))
                        g.add_legend() 
                        g.savefig(os.path.join(args.fig_dir, mod, "umap_method_" + qc + ".png"), dpi = 300)
                        plt.clf()
                    else:
                        L.info('skipping plot as too many categorys %s' % qc )

L.info('done')