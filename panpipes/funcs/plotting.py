import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scanpy.get import obs_df
from matplotlib.cm import get_cmap

import logging
from matplotlib import use
use('Agg')
plt.ioff()

def scatter_one(group_choice, col_choice, plot_df, axs=None, colour="#1f77b4", title=""):
    plot_df['plot_col'] = plot_df[col_choice]
    # highlight one choice
    other_choices = plot_df['plot_col'].unique().tolist()
    other_choices.remove(group_choice)
    # " nan" is an essential  cheat to get the nans to plot at the back
    plot_df['plot_col'] = plot_df['plot_col'].replace(other_choices, " nan")
    plot_df['plot_col'] = plot_df['plot_col'].astype("category")
    plot_df['plot_col'] = plot_df['plot_col'].cat.set_categories([group_choice, " nan"], ordered=True)
    # reorder plot so nans are plotted first
    plot_df = plot_df.sort_values(by='plot_col', ascending=False)
    # do plot
    if axs is None:
        fig, axs = plt.subplots(figsize=(3, 3), facecolor='w', edgecolor='k')
    sns.scatterplot(data=plot_df, x=plot_df.columns[0], y=plot_df.columns[1], hue="plot_col", palette=[colour, "#dddddd"],
                    s=2, alpha=.3, linewidth=0,
                    legend=False, ax=axs)
    axs.axis("off")
    axs.set_title(title)
    return axs


def batch_scatter_two_var(plot_df, method, batch, palette_choice=None):
    plot_df = plot_df[['umap_1', 'umap_2', method, batch]]
    # n_vars = len(df[var_choice].unique())
    nrows = len(plot_df[batch].unique())
    ncols = len(plot_df[method].unique())
    fig, axs = plt.subplots(nrows, ncols, figsize=(4*ncols, 3*nrows), facecolor='w', edgecolor='k')
    fig.subplots_adjust(hspace=.2, wspace=.2)
    axs = axs.ravel(order="F")
    if palette_choice is None:
        palette_choice = ['#1f77b4']*ncols
    # for j in range(ncols):
    idx=0
    for j, method_choice in enumerate(plot_df[method].cat.categories.tolist()):
        plot_df2 = plot_df[plot_df[method] == method_choice].copy()
        for i in range(nrows):
            group_choice = plot_df2[batch].unique()[i]
            scatter_one(group_choice, batch, plot_df2, axs[idx], colour=palette_choice[j], title=method_choice + "|" + group_choice)
            idx+=1

def facet_scatter(x, y, c, **kwargs):
    """Draw scatterplot with point colors from a faceted DataFrame columns."""
    kwargs.pop("color")
    plt.scatter(x, y, c=c, **kwargs)



def adjust_x_axis(ax):
    ax.set_xticklabels(ax.get_xticklabels(),
                       rotation= 90, 
                       horizontalalignment='center')
    ax.set(xlabel=None)


def _check_col_from_any_assay(mdata, var_col):
    if var_col in mdata.obs.columns:
        pass
    elif var_col in mdata['rna'].obs.columns: 
    # facet_row not found in obs, checking rna obs
        var_col = "rna:" + var_col
    elif var_col in mdata['prot'].obs.columns: 
    # facet_row not found in obs, checking rna obs
        var_col = "prot:" + var_col
    # finally check that it is not all NAs
    if mdata.obs[var_col].isna().sum() == mdata.shape[0]:
        logging.error("%s is all NaN, cannot plot" % var_col)
    

        
        

def scatter_fg_vs_bg(mdata, mdata_raw, x, y, facet_row=None):
    facet_row = _check_col_from_any_assay(mdata_raw, facet_row)
    logging.debug(facet_row)
    
    data_cols = [x for x in [x,y,facet_row] if x is not None]
    logging.debug(data_cols)
    df = pd.concat(
        [obs_df(mdata_raw, data_cols),
         obs_df(mdata, data_cols)],
        keys=["background", "cells"]).reset_index(col_level=0)
    logging.debug(df.head())
    df.rename(columns ={"level_0": "droplet_type", "level_1" : "barcode"}, inplace=True)
    g=sns.FacetGrid(data=df,  row=facet_row, col="droplet_type", height=3, aspect=1.6, margin_titles=True)
    # g.map(sns.scatterplot, x, y, s=0.2, color=".15")
    g.map(sns.histplot, x, y,cmap=sns.cubehelix_palette(start=.5, rot=-.5, as_cmap=True), bins=100 )
    ### Tidy up the labels
    # Iterate thorugh each axis
    for ax in g.axes.flat:
        # Make x and y-axis labels slightly larger
        ax.set_xlabel(ax.get_xlabel(), fontsize='large')
        ax.set_ylabel(ax.get_ylabel(), fontsize='large')

        # Make title more human-readable and larger
        if ax.get_title():
            ax.set_title(ax.get_title().split('=')[1],
                         fontsize='large')

        # Make right ylabel more human-readable and larger
        # Only the 2nd and 4th axes have something in ax.texts
        if ax.texts:
            # This contains the right ylabel text
            txt = ax.texts[0]
            ax.text(txt.get_unitless_position()[0], 
                    txt.get_unitless_position()[1],
                    txt.get_text().split('=')[1],
                    transform=ax.transAxes,
                    va='center',
                    fontsize='large'
                   )
            # Remove the original text
            ax.texts[0].remove()
    return g



from sklearn.neighbors import KernelDensity
def _kde_curve(xx, df, bandwidth=0.1):
    X=df.loc[:,xx].values[:, np.newaxis]
    kde = KernelDensity(kernel='gaussian', bandwidth=bandwidth).fit(X) 
    x=np.linspace(0,X.max().max(),100)[:, np.newaxis]
    log_density_values=kde.score_samples(x)
    density=np.exp(log_density_values).ravel()
    return x.ravel(), density

def ridgeplot(adata, features, layer=None, splitplot=3, bandwidth=0.1):
    """
    # code based on https://github.com/rougier/scientific-visualization-book/blob/master/code/anatomy/zorder-plots.py
    """
    # get data 
    if layer is not None:
        df = obs_df(adata, keys=features, layer=layer)
    else:
        df = obs_df(adata, keys=features)

    ncols = splitplot
    nrows = int(np.ceil(len(features)/ncols))

    np.random.seed(123)
    cmap = get_cmap("Spectral")
    
    fig_width = max([5, splitplot*5])
    fig_height = max([6, int(nrows/3)])
    fig = plt.figure(figsize=(fig_width, fig_height ))

    ax = None
    for n in range(ncols):
        ax = plt.subplot(1, ncols, n + 1, frameon=False, sharex=ax)
        n_feat = len(features)
        features_sub = features[n*nrows: n*nrows+nrows]
        # print(features_sub)
        for i, feat in enumerate(features_sub):
            X, Y = _kde_curve(feat, df, bandwidth=bandwidth)
            ax.plot(X, 1.5 * Y + i, color="k", linewidth=0.75, zorder=100 - i)
            color = cmap(i / len(features_sub))
            ax.fill_between(X, 1.5 * Y + i, i, color=color, zorder=100 - i)
        # if layer == "clr":
        #     ax.set_xlim(0, 3)

        ax.yaxis.set_tick_params(tick1On=False)
        ax.set_ylim(-1, len(features_sub) + 3)

        ax.yaxis.set_tick_params(labelleft=True)
        ax.set_yticks(np.arange(len(features_sub)))
        ax.set_yticklabels(features_sub)
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(10)
            tick.label.set_verticalalignment("bottom")

    fig.tight_layout()
    return fig, ax