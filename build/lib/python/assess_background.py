import argparse
import pandas as pd
import os
import numpy as np
import scanpy as sc
import muon as mu
from matplotlib import use
use('Agg')
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Agg")
import seaborn as sns
import panpipes.funcs as pnp


import sys
import logging
L = logging.getLogger()
L.setLevel(logging.INFO)
log_handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s: %(levelname)s - %(message)s')
log_handler.setFormatter(formatter)
L.addHandler(log_handler)
L.debug("test logging message")
# parse arguments
parser = argparse.ArgumentParser()


parser.add_argument('--channel_col',
                    default="sample_id",
                    help='')
parser.add_argument('--filtered_mudata',
                    default=None,
                    help='')  
parser.add_argument('--raw_mudata',
                    default=None,
                    help='')    
parser.add_argument('--figpath',
                    default="./figures",
                    help='')  


args, opt = parser.parse_known_args()
L.info(args)
# load filtered data - if use_umon is True this will return a mudata object, else returns an mudata object

L.info("reading filtered mudata object")
try:
    mdata = mu.read(args.filtered_mudata)
    # need to delete rep else we lose too many cells in intersect_obs later
    if "rep" in mdata.mod.keys():
        del mdata.mod["rep"]
except FileNotFoundError:
    sys.exit("filtered mudata object not found")


# mu.pp.intersect_obs(mdata)
# load raw data

L.info("reading raw mudata object")
try:
    mdata_raw = mu.read(args.raw_mudata)
except FileNotFoundError:
    sys.exit("raw_mudata object not found")


if 'rna' in mdata.mod.keys():
    # remove cells that passed cellranger filtering from raw
    mdata_raw['rna'].obs['is_empty'] = ~mdata_raw['rna'].obs_names.isin(mdata['rna'].obs_names)
    mu.pp.filter_obs(mdata_raw['rna'], 'is_empty')

# mdata_raw.update()
# mu.pp.intersect_obs(mdata_raw)

# calculate qc metrics
if 'rna' in mdata.mod.keys():
    # run stadard qc metrics
    sc.pp.calculate_qc_metrics(mdata['rna'], inplace=True, percent_top=None, log1p=True)
    sc.pp.calculate_qc_metrics(mdata_raw['rna'], inplace=True, percent_top=None, log1p=True)



if 'prot' in mdata.mod.keys():
    sc.pp.calculate_qc_metrics(mdata['prot'], inplace=True, percent_top=None, log1p=True)
    sc.pp.calculate_qc_metrics(mdata_raw['prot'], inplace=True, percent_top=None, log1p=True)

mdata.update()

mdata_raw.update()


## QC for gex and protein fcomparing foreground and background
if os.path.exists(args.figpath) is False:
    os.mkdirs(args.figpath)

if 'rna' in mdata.mod.keys():
    pnp.pl.scatter_fg_vs_bg(mdata, mdata_raw,x="rna:log1p_n_genes_by_counts", y="rna:log1p_total_counts", facet_row=args.channel_col)
    plt.subplots_adjust(right=0.7)
    plt.savefig(os.path.join(args.figpath,"scatter_bg_fg_rna_nGene_rna_nUMI.png"),transparent=False)
    if 'prot' in mdata.mod.keys() :
        pnp.pl.scatter_fg_vs_bg(mdata, mdata_raw,x="prot:log1p_total_counts", y="rna:log1p_n_genes_by_counts", facet_row=args.channel_col)
        plt.subplots_adjust(right=0.7)
        plt.savefig(os.path.join(args.figpath, "scatter_bg_fg_prot_nUMI_rna_nGene.png"),transparent=False)
        pnp.pl.scatter_fg_vs_bg(mdata, mdata_raw,x="rna:log1p_total_counts", y="prot:log1p_total_counts", facet_row=args.channel_col)
        plt.subplots_adjust(right=0.7)
        plt.savefig(os.path.join(args.figpath,"scatter_bg_fg_rna_nUMI_prot_nUMI.png"),transparent=False)

# quantifying the top background features
if 'rna' in mdata.mod.keys():
    sc.pl.highest_expr_genes(mdata_raw['rna'],n_top=30, save="_rna_background.png")
    top_genes = pnp.scmethods.get_top_expressed_features(mdata_raw['rna'], n_top=30, group_by=args.channel_col)
    bg_df = pnp.scmethods.get_mean_background_fraction(mdata_raw['rna'], top_background=top_genes, group_by=args.channel_col)
    fig, ax= plt.subplots(figsize=(12,8))
    sns.heatmap(bg_df, ax=ax)
    fig.tight_layout()
    plt.savefig(os.path.join(args.figpath,"heatmap_background_" + args.channel_col + "_rna_top_expressed.png"))


## Repeat for protein (if it exists)
if 'prot' in mdata.mod.keys():
    # this time we'll just use all the adts.
    top_genes = list(mdata_raw['prot'].var_names)
    bg_df = pnp.scmethods.get_mean_background_fraction(mdata_raw['prot'], top_background=top_genes, group_by=args.channel_col)
    if len(top_genes) > 50:
        # split the plot into two
        split_int = int(np.ceil(137/2))
        fig, ax = plt.subplots(nrows=2,ncols=1,figsize=(24,10), facecolor="white")
        sns.heatmap(bg_df.iloc[:,1:split_int], ax=ax[0])
        sns.heatmap(bg_df.iloc[:,split_int:len(top_genes)], ax=ax[1])
        fig.suptitle("mean exprs (raw counts) of ADT in empty drops")  
        fig.tight_layout()
    else:
        fig, ax = plt.subplots(figsize=(12,10), facecolor="white")
        sns.heatmap(bg_df, ax=ax)
        fig.tight_layout()
    plt.savefig(os.path.join(args.figpath,"heatmap_background_" + args.channel_col + "_prot_top_expressed.png"))

