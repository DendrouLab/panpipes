import pandas as pd
import numpy as np
# from scanpy.pp import normalize_total
from scipy.sparse import issparse
from scanpy.get import obs_df as get_obs_df
from scanpy.pp import normalize_total
from typing import Optional, Literal
import muon as mu
from muon import MuData
from anndata import AnnData
import os
from matplotlib.pyplot import savefig

from .plotting import ridgeplot
from .io import write_10x_counts
from .processing import check_for_bool, which_ind



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



def pseudo_seurat(adata, arg_minpct=0.1, arg_mindiffpct=-float("inf"), arg_logfcdiff=0.25, use_dense=False):
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
        print("using dense matrix")
        # extract the counts for cluster cells and calculate exp means on each row
        nct = adata.X.T[:, cluster_cells_ind]
        cluster_mean = np.apply_along_axis(exp_mean_dense, 1, nct.todense())

        # likewise for non-cluster cells
        nct = adata.X.T[:, other_cells_ind]
        other_mean = np.apply_along_axis(exp_mean_dense, 1, nct.todense())
        diff_mean = abs(cluster_mean - other_mean)
    else:
        print("using sparse matrix")
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
    if method == "scanpy":
        print("Computing neighbors using scanpy")
        from scanpy.pp import neighbors
        neighbors(adata,
                        n_pcs=n_pcs,
                        n_neighbors=n_neighbors,
                        metric=metric,
                        use_rep=use_rep)
    elif method == "hnsw":
        from scvelo.pp import neighbors
        print("Computing neighbors using hnswlib (with scvelo a la pegasus!)")
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
        return adata

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
                      mdata_raw, 
                      method: Literal['dsb', 'clr'] = 'clr',
                      clr_margin: Literal[0, 1] = '0',
                      isotypes=None):
    # rplace raw counts back in X
    # check if integers in mdata['prot'].X
    if X_is_raw(mdata['prot']) is False:
        raise ValueError("mdata['prot'].X does not contain raw counts, cannot normalise")
    if method == "clr":
        mu.prot.pp.clr(mdata["prot"], inplace=True, axis=int(clr_margin))
        mdata["prot"].layers["clr"] = mdata["prot"].X.copy()
    elif method == "dsb":
        # check mdata_raw actually contains raw counts
        mu.prot.pp.dsb(mdata,
                    mdata_raw, 
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
    