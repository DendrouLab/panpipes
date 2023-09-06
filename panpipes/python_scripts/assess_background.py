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
from scipy.stats import rankdata

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
parser.add_argument('--bg_mudata',
                    default=None,
                    help='')    
parser.add_argument('--figpath',
                    default="./figures/background",
                    help='')  


args, opt = parser.parse_known_args()
L.info(args)

sc.settings.figdir = args.figpath
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
# load bg data

L.info("reading bg mudata object")
try:
    mdata_bg = mu.read(args.bg_mudata)
except FileNotFoundError:
    sys.exit("bg_mudata object not found")

# barcode rank plots
df = mdata_bg['rna'].obs.copy()
df['total_counts'] = np.ravel(mdata_bg['rna'].X.sum(axis=1))
df['log_total_counts'] = np.log(df['total_counts'])
df['ranks'] = df.groupby('sample_id', as_index=False)['total_counts'].transform(lambda x: rankdata(x*-1, method="ordinal"))
df['log_ranks'] = np.log(df['ranks'])
df['is_cell'] = df.index.isin(mdata['rna'].obs_names.tolist())
df['is_cell'] = df['is_cell'].astype('category')
df['is_cell'] = df['is_cell'].cat.set_categories([True, False])

fig, ax= plt.subplots(figsize=(6,6))
sns.scatterplot(data=df, x="log_ranks", y="log_total_counts", 
                hue=args.channel_col, size="is_cell",  sizes=(2, 20),
                ax=ax, linewidth=0)
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.savefig(os.path.join(args.figpath,'rna_barcode_ranks.png'), bbox_inches="tight")
plt.clf()

n_samples = len(df[args.channel_col].unique())
if n_samples > 1:
    g = sns.FacetGrid(df, col=args.channel_col, 
                    height=2.5,  
                    col_wrap=np.min([n_samples,8]))
    g.map(sns.scatterplot, "log_ranks", "log_total_counts", linewidth=0, s=1.5)
    g.set_titles(col_template="{col_name}")
    plt.savefig(os.path.join(args.figpath,'rna_facet_barcode_ranks.png'))
    plt.clf()



if 'rna' in mdata.mod.keys():
    # remove cells that passed cellranger filtering from bg
    mdata_bg['rna'].obs['is_empty'] = ~mdata_bg['rna'].obs_names.isin(mdata['rna'].obs_names)
    mu.pp.filter_obs(mdata_bg['rna'], 'is_empty')

# mdata_bg.update()
# mu.pp.intersect_obs(mdata_bg)

# calculate qc metrics
if 'rna' in mdata.mod.keys():
    # run stadard qc metrics
    sc.pp.calculate_qc_metrics(mdata['rna'], inplace=True, percent_top=None, log1p=True)
    sc.pp.calculate_qc_metrics(mdata_bg['rna'], inplace=True, percent_top=None, log1p=True)



if 'prot' in mdata.mod.keys():
    sc.pp.calculate_qc_metrics(mdata['prot'], inplace=True, percent_top=None, log1p=True)
    sc.pp.calculate_qc_metrics(mdata_bg['prot'], inplace=True, percent_top=None, log1p=True)


mdata.update()
mdata_bg.update()

n_samples_rna=len(mdata['rna'].obs['sample_id'].unique())
n_samples_prot=len(mdata['prot'].obs['sample_id'].unique())

# quantifying the top background features
if os.path.exists(args.figpath) is False:
    os.mkdirs(args.figpath)

if 'rna' in mdata.mod.keys():
    sc.pl.highest_expr_genes(mdata_bg['rna'],n_top=30, save="_rna_background.png")
    top_genes = pnp.scmethods.get_top_expressed_features(mdata_bg['rna'], n_top=30, group_by=args.channel_col)
    bg_df = pnp.scmethods.get_mean_background_fraction(mdata_bg['rna'], top_background=top_genes, group_by=args.channel_col)
    if bg_df.shape[0] >1:
        fig, ax= plt.subplots(figsize=(12,8))
        sns.heatmap(bg_df, ax=ax)
        fig.tight_layout()
        plt.savefig(os.path.join(args.figpath,"heatmap_background_" + args.channel_col + "_rna_top_expressed.png"))
    else:
        fig, ax= plt.subplots(figsize=(12,8))
        plt_df = bg_df.stack().reset_index().rename(columns={'level_1':'feature', 0:'exprs'})
        sns.barplot(data=plt_df, x='feature', y='exprs', ax=ax)
        pnp.plotting.adjust_x_axis(ax)
        fig.tight_layout()
        plt.savefig(os.path.join(args.figpath,"barplot_background_" + args.channel_col + "_rna_top_expressed.png"))


# quantifying the top background features
## Repeat for protein (if it exists)
if 'prot' in mdata.mod.keys():
    # this time we'll just use all the adts.
    top_genes = list(mdata_bg['prot'].var_names)
    bg_df = pnp.scmethods.get_mean_background_fraction(mdata_bg['prot'], top_background=top_genes, group_by=args.channel_col)
    if bg_df.shape[0] >1:
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
    else:
        plt_df = bg_df.stack().reset_index().rename(columns={'level_1':'feature', 0:'exprs'})
        if len(top_genes) > 50:
            split_int = int(np.ceil(137/2))
            fig, ax = plt.subplots(nrows=2,ncols=1,figsize=(24,10), facecolor="white")
            sns.heatmap(plt_df.iloc[1:split_int,:], ax=ax[0])
            pnp.plotting.adjust_x_axis(ax[0])
            sns.heatmap(plt_df.iloc[split_int:len(top_genes),:], ax=ax[1])
            pnp.plotting.adjust_x_axis(ax[1])
            fig.suptitle("mean exprs (raw counts) of ADT in empty drops")  
            fig.tight_layout()
        else:
            fig, ax= plt.subplots(figsize=(12,8))
            sns.barplot(data=plt_df, x='feature', y='exprs', ax=ax)
            pnp.plotting.adjust_x_axis(ax)
            fig.tight_layout()
        plt.savefig(os.path.join(args.figpath,"barplot_background_" + args.channel_col + "_prot_top_expressed.png"))


## QC for gex and protein fcomparing foreground and background

if 'rna' in mdata.mod.keys():
    pnp.pl.scatter_fg_vs_bg(mdata, mdata_bg,x="rna:log1p_n_genes_by_counts", y="rna:log1p_total_counts", facet_row=args.channel_col)
    plt.subplots_adjust(right=0.7)
    plt.savefig(os.path.join(args.figpath,"scatter_bg_fg_rna_nGene_rna_nUMI.png"),transparent=False)
    if 'prot' in mdata.mod.keys():
        if n_samples_rna != n_samples_prot:
            L.info("n_samples are not equal in rna and prot,taking intrsect for mdata and mdata_bg")
            mu.pp.intersect_obs(mdata_bg)
            mu.pp.intersect_obs(mdata)
            pnp.pl.scatter_fg_vs_bg(mdata, mdata_bg,x="prot:log1p_total_counts", y="rna:log1p_n_genes_by_counts", facet_row=args.channel_col)
            plt.subplots_adjust(right=0.7)
            plt.savefig(os.path.join(args.figpath, "scatter_bg_fg_prot_nUMI_rna_nGene.png"),transparent=False)
            pnp.pl.scatter_fg_vs_bg(mdata, mdata_bg,x="rna:log1p_total_counts", y="prot:log1p_total_counts", facet_row=args.channel_col)
            plt.subplots_adjust(right=0.7)
            plt.savefig(os.path.join(args.figpath,"scatter_bg_fg_rna_nUMI_prot_nUMI.png"),transparent=False)
        else:
            L.info("n_samples are equal in rna and prot")
            pnp.pl.scatter_fg_vs_bg(mdata, mdata_bg,x="prot:log1p_total_counts", y="rna:log1p_n_genes_by_counts", facet_row=args.channel_col)
            plt.subplots_adjust(right=0.7)
            plt.savefig(os.path.join(args.figpath, "scatter_bg_fg_prot_nUMI_rna_nGene.png"),transparent=False)
            pnp.pl.scatter_fg_vs_bg(mdata, mdata_bg,x="rna:log1p_total_counts", y="prot:log1p_total_counts", facet_row=args.channel_col)
            plt.subplots_adjust(right=0.7)
            plt.savefig(os.path.join(args.figpath,"scatter_bg_fg_rna_nUMI_prot_nUMI.png"),transparent=False)



plt.clf()


L.info("Done")

