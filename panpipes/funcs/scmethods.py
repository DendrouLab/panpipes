import pandas as pd
import numpy as np
# from scanpy.pp import normalize_total
from scipy.sparse import issparse
from scanpy.get import obs_df as get_obs_df
from scanpy.pp import normalize_total
import scanpy as sc
import warnings
import logging
from typing import Optional, Literal
import sys
import muon as mu
from muon import MuData
from anndata import AnnData
import os
from matplotlib.pyplot import savefig
from statsmodels.distributions.empirical_distribution import ECDF

from .plotting import ridgeplot
from .io import write_10x_counts
from .processing import check_for_bool, which_ind
import matplotlib
import matplotlib.pyplot as plt



def cell2loc_filter_genes(adata, fig_path, cell_count_cutoff=15, cell_percentage_cutoff2=0.05, nonz_mean_cutoff=1.12):
# Adjusted from: https://github.com/BayraktarLab/cell2location/blob/master/cell2location/utils/filtering.py
    adata.var["n_cells"] = np.array((adata.X > 0).sum(0)).flatten()
    adata.var["nonz_mean"] = np.array(adata.X.sum(0)).flatten() / adata.var["n_cells"]

    cell_count_cutoff = np.log10(cell_count_cutoff)
    cell_count_cutoff2 = np.log10(adata.shape[0] * cell_percentage_cutoff2)
    nonz_mean_cutoff = np.log10(nonz_mean_cutoff)

    gene_selection = (np.array(np.log10(adata.var["n_cells"]) > cell_count_cutoff2)) | (
        np.array(np.log10(adata.var["n_cells"]) > cell_count_cutoff)
        & np.array(np.log10(adata.var["nonz_mean"]) > nonz_mean_cutoff)
    )
    gene_selection = adata.var_names[gene_selection]
    adata_shape = adata[:, gene_selection].shape

    fig, ax = plt.subplots()
    ax.hist2d(
        np.log10(adata.var["nonz_mean"]),
        np.log10(adata.var["n_cells"]),
        bins=100,
        norm=matplotlib.colors.LogNorm(),
        range=[[0, 0.5], [1, 4.5]],
    )
    ax.axvspan(0, nonz_mean_cutoff, ymin=0.0, ymax=(cell_count_cutoff2 - 1) / 3.5, color="darkorange", alpha=0.3)
    ax.axvspan(
        nonz_mean_cutoff,
        np.max(np.log10(adata.var["nonz_mean"])),
        ymin=0.0,
        ymax=(cell_count_cutoff - 1) / 3.5,
        color="darkorange",
        alpha=0.3,
    )
    plt.vlines(nonz_mean_cutoff, cell_count_cutoff, cell_count_cutoff2, color="darkorange")
    plt.hlines(cell_count_cutoff, nonz_mean_cutoff, 1, color="darkorange")
    plt.hlines(cell_count_cutoff2, 0, nonz_mean_cutoff, color="darkorange")
    plt.xlabel("Mean non-zero expression level of gene (log)")
    plt.ylabel("Number of cells expressing gene (log)")
    plt.title(f"Gene filter: {adata_shape[0]} cells x {adata_shape[1]} genes")
    plt.savefig(fig_path)

    return gene_selection




def findTopFeatures_pseudo_signac(adata, min_cutoff):
    # Adapted from:
    # https://stuartlab.org/signac/reference/findtopfeatures
    """
    :param adata:
    :param min_cutoff:
        "q[x]": "q" followed by the minimum percentile, e.g. q5 will set the top 95% most common features as higly variable
        "c[x]": "c" followed by a minimum cell count, e.g. c100 will set features present in > 100 cells as highly variable
        "tc[x]": "tc" followed by a minimum total count, e.g. tc100 will set features with total counts > 100 as highly variable
        "NULL": All features are assigned as highly variable
        "NA": Highly variable features won't be changed
    :return: Percentile rank of each feature stored in "adata.var["percentile"], highly variable features stored in adata.var["highly_variable"].
    "total_counts" and "n_cells_by_counts" also saved, if not present.
    Filtering is done according to the "min_cutoff" parameter. If min_cutoff == NA, adata.var["highly_variable"] won't be changed
    """
    if "total_counts" not in adata.var:
        adata.var["total_counts"] = np.ravel(adata.layers["raw_counts"].sum(axis=0))
    if "n_cells_by_counts" not in adata.var:
        if issparse(adata.layers["raw_counts"]):
            adata.var["n_cells_by_counts"] = adata.layers["raw_counts"].getnnz(axis=0)
        else:
            adata.var["n_cells_by_counts"] = np.count_nonzero(adata.layers["raw_counts"], axis=0)


    if "percentile" not in adata.var:
        # fit ecdf
        ecdf = ECDF(adata.var["total_counts"])
        extracted = adata.var[["total_counts"]]
        extracted = extracted.sort_values("total_counts", ascending=False)
        # add decreasingly ordered percentiles to dataframe:
        extracted["percentile"] = -np.sort(-ecdf.y)[:-1]  # last percentile always corresponds to the percentile of -Inf
        # merge with adata.var
        adata.var = pd.merge(adata.var, extracted.drop("total_counts", axis=1), left_index=True, right_index=True)

    if min_cutoff.startswith("q"):  # quantile filtering
        min_cutoff = int(min_cutoff[1:]) / 100
        adata.var['highly_variable'] = [True if percentile > min_cutoff else False for percentile in adata.var["percentile"]]
        return
    if min_cutoff.startswith("c"):  # filtering by minimum number of cells
        min_cutoff = int(min_cutoff[1:])
        adata.var['highly_variable'] = [True if ncells > min_cutoff else False for ncells in adata.var["n_cells_by_counts"]]
        return
    if min_cutoff.startswith("tc"):  # filtering by total counts
        min_cutoff = int(min_cutoff[2:])
        adata.var['highly_variable'] = [True if totalcounts > min_cutoff else False for totalcounts in adata.var["total_counts"]]
        return
    if min_cutoff == "NA":  # don't change variable features
        return
    if min_cutoff == "NULL":  # assign all features as highly variable
        adata.var["highly_variable"] = True
        return



def exp_mean_sparse(x):
    """
    This returns the log of the mean of the not-logged data
    this version of the function works directly on the sparse matrix but is slow.
    """
    # convert out of compressed sparse matrix
    return np.log(x.expm1().mean(1)+1).T


def exp_mean_dense(x):
    """
    This returns the log of the mean of the not-logged data
    this version of the function requires a dense matrix, so might be memory hungry?
    But it is super fast
    """
    # convert out of compressed sparse matrix
    return np.log((np.sum(np.exp(x)-1)/x.shape[1]) + 1)

def find_all_markers_pseudo_seurat(
        adata, 
        groups,
        groupby,
        layer=None,
        method=None,
        n_genes=float("inf"), 
        corr_method="bonferroni",
        arg_minpct=0.1,
        arg_mindiffpct=-float("inf"), 
        arg_logfcdiff=0.25):
    # add replace X with layer
    if layer is not None:
        adata.X = adata.layers[layer]
    # need to check is the assay layer is dense or not
    assay_is_sparse = issparse(adata.X)
    use_dense = assay_is_sparse==False 
    if groups == 'all':
        groups = adata.obs[groupby].unique().tolist()
    markers_dict = {}
    filter_dict = {}
    for cv in groups:
        # \ set up idenst as cv ==1 and everything else = 0
        adata.obs['idents'] = ['1' if x == cv else '0' for x in adata.obs[groupby]]
        filter_dict[cv] = pseudo_seurat(adata, use_dense=use_dense,arg_minpct=arg_minpct,
                  arg_mindiffpct=arg_mindiffpct, 
                  arg_logfcdiff=arg_logfcdiff )
        logging.info("number of genes remaining after filtering:  %i\n" % filter_dict[cv]['background'].sum())
        adata_rg = adata[:, filter_dict[cv]['background'].tolist()].copy()
        sc.tl.rank_genes_groups(adata_rg, layer=layer,
                                groupby="idents", groups=["1"],  
                                reference="0",
                                method=method, 
                                n_genes=float("inf"), 
                                corr_method="bonferroni")
        markers_dict[cv] = sc.get.rank_genes_groups_df(adata_rg, group="1")
        # remove adata from mem
        adata_rg = None
    markers = pd.concat(markers_dict.values(), keys=markers_dict.keys())
    filter_stats = pd.concat(filter_dict.values(), keys=filter_dict.keys())
    return markers, filter_stats

def pseudo_seurat(adata, 
                  arg_minpct=0.1,
                  arg_mindiffpct=-float("inf"), 
                  arg_logfcdiff=0.25, 
                  use_dense=False):
    """
    alternative method that"s more like seurat (pseudo seurat if you will)
    In that you filter genes before running rank genes
    ---
    1.  define idents
    2. define cluster cells and other cells
    3. check min cells
    4. compute percentages and difference from equiv data slot to s.obj@assays$RNA@data
    5. computer mean expression levels and differences
    6. save stats
    7. define background based on min_pct (pct of cells that have to express marker)
    8. filter stats for testing and save universe
    genes|cluster_mean|other_mean|diff_mean|cluster_pct|other_pct|max_pct|min_pct|diff_pct|background
    9. Filter adata.raw based on these genes?

    9. Find markers based on cluster cells and other cells,
    values from adata.raw ??
    gene: only expressed
    cells: we need to define manually based on above stats.
    10. save results.
    """
    # define cells
    cluster_cells_ind = which_ind(adata.obs["idents"] == "1")
    other_cells_ind = which_ind(adata.obs["idents"] == "0")

    # compute perecentage expressed
    # from normnalised but not scaled data
    # remember cells are rows and genes are columns

    # note: I don't know why norm_counts[cluster_cell_ind:, col_ind] deosn"t work, but it doesn't
    cluster_pct = (adata.X[cluster_cells_ind, :] > 0).sum(axis=0) / len(cluster_cells_ind)
    other_pct = (adata.X[other_cells_ind, :] > 0).sum(axis=0) / len(other_cells_ind)

    pcts = pd.DataFrame(np.vstack((cluster_pct,  other_pct)).transpose())
    max_pct = pcts.max(axis=1)
    min_pct = pcts.min(axis=1)
    diff_pct = max_pct - min_pct
    take_diff_pct = diff_pct > arg_mindiffpct
    # remove genes that are not expressed higher than 0.1 in one of the groups
    take_min_pct = max_pct > arg_minpct

    # KEEP IN CASE NP.ARRAY METHOD USES TOO MUCH MEMORY
    # import time
    # this has the potential to be very slow. Transposeing it speeds it up a bit.
    # I need to undertand sparse matrices better to make it work
    if use_dense:
        logging.info("using dense matrix")
        # extract the counts for cluster cells and calculate exp means on each row
        nct = adata.X.T[:, cluster_cells_ind]
        cluster_mean = np.apply_along_axis(exp_mean_dense, 1, nct.todense())

        # likewise for non-cluster cells
        nct = adata.X.T[:, other_cells_ind]
        other_mean = np.apply_along_axis(exp_mean_dense, 1, nct.todense())
        diff_mean = abs(cluster_mean - other_mean)
    else:
        logging.info("using sparse matrix")
        cluster_mean = exp_mean_sparse(adata.X.T[:, cluster_cells_ind])
        other_mean = exp_mean_sparse(adata.X.T[:, other_cells_ind])
        diff_mean = abs(cluster_mean - other_mean).A1

    # remove genes with less than threshold difference
    take_thresh = diff_mean > arg_logfcdiff
    # take = if a cell passes all the tests then it is to be kept.
    take = [a and b and c for a, b, c in zip(take_thresh, take_min_pct, take_diff_pct)]
    print("saving universe for fisher test")
    stats_df = pd.DataFrame(np.vstack((adata.var_names, cluster_mean, other_mean, diff_mean,
                                       cluster_pct, other_pct, max_pct, min_pct, diff_pct, take)).transpose(),
                            columns=["gene", "cluster_mean", "other_mean", "diff_mean",
                                     "cluster_pct", "other_pct",
                                     "max_pct", "min_pct", "diff_pct", "background"])
    return stats_df



def run_neighbors_method_choice(adata, method, n_neighbors, n_pcs, metric, use_rep, nthreads=1):
    # This works with both Anndata and MuData inputs 
    # useful if we are dealing with a MuData object but we want to use single rep, e.g.
    # calculating neighbors on a totalVI latent rep
    if method == "scanpy":
        logging.info("Computing neighbors using scanpy")
        from scanpy.pp import neighbors
        neighbors(adata,
                        n_pcs=n_pcs,
                        n_neighbors=n_neighbors,
                        metric=metric,
                        use_rep=use_rep)
    elif method == "hnsw":
        from scvelo.pp import neighbors
        logging.info("Computing neighbors using hnswlib (with scvelo a la pegasus!)")
        # we use the neighbors function from scvelo (thanks!)
        # with parameters from pegasus (for a more exact result).
        # code snippet from Steve Sansom, via COMBAT project
        neighbors(adata,
                n_neighbors=int(n_neighbors),
                n_pcs=n_pcs,
                use_rep=use_rep,
                knn=True,
                random_state=0,
                method='hnsw',
                metric=metric,
                metric_kwds={"M": 20,
                            "ef": 200,
                            "ef_construction": 200},
                num_threads=int(nthreads))


def merge_consensus_clust(adata, consensus_clust, ref_col="rough_ref"):
    if ref_col in adata.obs.columns:
        out = adata
    else:
        xx = consensus_clust[consensus_clust.index.isin(adata.obs_names)]
        adata.obs = adata.obs.merge(xx,how="left", left_index=True, right_index=True)
        adata.obs[ref_col] = adata.obs[ref_col].astype('category') 
        out = adata
    return(out)


# assessing background

def _calc_top_n_genes(adata, n_top=50):
    norm_dict = normalize_total(adata, target_sum=100, inplace=False)
    
    if issparse(norm_dict['X']):
            mean_percent = norm_dict['X'].mean(axis=0).A1
            top_idx = np.argsort(mean_percent)[::-1][:n_top]
            counts_top_genes = norm_dict['X'][:, top_idx].A
    else:
        mean_percent = norm_dict['X'].mean(axis=0)
        top_idx = np.argsort(mean_percent)[::-1][:n_top]
        counts_top_genes = norm_dict['X'][:, top_idx]
    columns = (
        adata.var_names[top_idx]
    )
    counts_top_genes = pd.DataFrame(
        counts_top_genes, index=adata.obs_names, columns=columns
    )
    
    return list(counts_top_genes.columns)

def get_top_expressed_features(adata, n_top=50, group_by=None):
    """
    Get highest expressed dfeatures from an andata object.
    If group by is specfied the union of the top features 
    from each grouping var will be returned
    
    Parameters
    ----------
    adata : AnnData object
    n_top : int
              No. of features to return (per grouping)
    group_by : str
              Designates a column from adata.obs on which to group
              vars.
    """
    if group_by is not None:
        groups = adata.obs[group_by].unique()
        top_genes_list = [_calc_top_n_genes(adata[adata.obs[group_by] == x, :], n_top=n_top) for x in groups]
        from itertools import chain
        top_genes = list(set(chain.from_iterable(top_genes_list)))
    else:
        top_genes = _calc_top_n_genes(adata, n_top=n_top)
    return(top_genes)

    
def get_mean_background_fraction(adata, top_background, group_by=None):
    norm_dict = normalize_total(adata, target_sum=100, inplace=False)
    gene_idx = [idx for idx, gene in  enumerate(adata.var_names ) if gene in top_background]
    if issparse(norm_dict['X']):
        mean_percent = norm_dict['X'].mean(axis=0).A1
        counts_top_genes = norm_dict['X'][:, gene_idx].A
    else:
        mean_percent = norm_dict['X'].mean(axis=0)
        counts_top_genes = norm_dict['X'][:, gene_idx]

    columns = (
        adata.var_names[gene_idx]
    )
    counts_top_genes = pd.DataFrame(
        counts_top_genes, index=adata.obs_names, columns=columns
    )
    if group_by is None:
        means_df = counts_top_genes.mean()
    else:
        obs_df = get_obs_df(adata, [group_by])
        counts_top_genes = pd.merge(counts_top_genes, obs_df, left_index=True, right_index=True)
        means_df = counts_top_genes.groupby(group_by).mean()
    return means_df


def identify_isotype_outliers(prot, isotypes, quantile_val=0.9, n_isotypes_pass=2, groupby=None,layer=None, inplace=True):
    # make sure we don't mess
    if inplace is False:
        prot=prot.copy()
    if groupby is None:
            isotypes_counts = get_obs_df(prot, isotypes, layer=layer)
            # get quantile thresholds
            thres = isotypes_counts.quantile(quantile_val)
            # for each isotype test whether the threshold is exceeded
            for iso in isotypes:
                isotypes_counts['outliers_' + iso]= isotypes_counts[iso] > thres[iso] 
    else:
        # we are grouping by an extra variable
        isotypes_counts = get_obs_df(prot, [groupby] + isotypes, layer=layer)
        thres = isotypes_counts.groupby(groupby).quantile(0.9)
        for iso in isotypes:
            tmp = pd.merge(isotypes_counts[[groupby, iso]], thres[[iso]], left_on=groupby, right_index=True)
            isotypes_counts['outliers_' + iso] = tmp.apply(lambda x: x[1] > x[2] , axis=1)
    # now work out which ones to exclude
    isotypes_counts['outliers'] = isotypes_counts.iloc[:,isotypes_counts.columns.str.startswith("outliers")].sum(axis=1)
    isotypes_counts['exclude'] = isotypes_counts['outliers'] > n_isotypes_pass
    # add into obs
    prot.obs['isotype_exclude_outliers'] = isotypes_counts['exclude']
    prot.obs['n_isotype_in_90_percentile'] = isotypes_counts['outliers'].astype('category')
    if inplace is False:
        return prot

def X_is_raw(adata):
    # we assume if all the rows sums are integers then the matrix is populated by integers
    return np.array_equal(adata.X.sum(axis=0).astype(int), adata.X.sum(axis=0))

# run clr 
    

def run_adt_normalise(mdata, 
                      mdata_bg, 
                      method: Literal['dsb', 'clr'] = 'clr',
                      clr_margin: Literal[0, 1] = '0',
                      isotypes=None):
    # rplace bg counts back in X
    # check if integers in mdata['prot'].X
    if X_is_raw(mdata['prot']) is False:
        raise ValueError("mdata['prot'].X does not contain raw counts, cannot normalise")
    if method == "clr":
        mu.prot.pp.clr(mdata["prot"], inplace=True, axis=int(clr_margin))
        mdata["prot"].layers["clr"] = mdata["prot"].X.copy()
    elif method == "dsb":
        # check mdata_bg actually contains raw counts
        mu.prot.pp.dsb(mdata,
                    mdata_bg, 
                    empty_counts_range=(1.0, 2.7), 
                    isotype_controls=isotypes, 
                    random_state=4)

        mdata["prot"].layers["dsb"] = mdata["prot"].X.copy()
    else:
        raise ValueError("normalisation method %s not recognised, please choose dsb or clr" % method)
    

from scipy.sparse import issparse
def quantile_clipping(mdata, modality=None, layer=None,  qmax=0.995, qmin=0.005, inplace=False):
    if isinstance(mdata, MuData):
        # need to pull out the modality and layer
        if modality is None:
            raise ValueError("if passing a MuData object, please specify a modality name e.g. prot")
        if layer is None:
            arr = mdata[modality].X
        else:
            arr = mdata[modality].layers[layer]
        if issparse(arr):
            arr = arr.todense()
        df = pd.DataFrame(data=arr, index=mdata[modality].obs_names.tolist(), columns=mdata[modality].var_names.tolist())
    elif isinstance(mdata, AnnData):
        if layer is None:
            arr = mdata.X
        else:
            arr = mdata.layers[layer] 
        if issparse(arr):
            arr = arr.todense()
        df = pd.DataFrame(data=arr, index=mdata.obs_names.tolist(), columns=mdata.var_names.tolist())
    # do quantile clipping
    qlow = np.quantile(arr, qmin, axis=0)
    qhigh = np.quantile(arr, qmax, axis=0)  
    df_clipped = df.clip(lower = qlow, upper = qhigh, axis=1).to_numpy()
    # return results
    if inplace:
        if isinstance(mdata, MuData):
        # need to pull out the modality and layer
            if layer is None:
                mdata[modality].X = df_clipped
            else:
                mdata[modality].layers[layer] = df_clipped
        elif isinstance(mdata, AnnData):
            if layer is None:
                mdata.X = df_clipped
            else:
                mdata.layers[layer] = df_clipped
    else:
        return df_clipped
    