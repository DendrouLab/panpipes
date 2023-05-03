from collections import UserString
import sys
import warnings
import logging
import os 
import re
import numpy as np
import pandas as pd
import muon as mu

from pandas import concat
from pandas import merge as pdmerge
from itertools import compress, chain
from scanpy.get import obs_df
from scanpy.pp import subsample
from muon import MuData
from anndata import AnnData
# from muon.pp import intersect_obs
from random import sample
from functools import reduce


def is_float_try(string):
    """
    quick test to see if a string can be turned intoa float
    """
    try:
        float(string)
        return True
    except ValueError:
        return False


def splitall(path):
    """
    Function to completely split up a path into it's compoenents
    e.g. /path/to/cheese => ['path', 'to', 'cheese'] 
    inspired by https://stackoverflow.com/questions/3167154/how-to-split-a-dos-path-into-its-components-in-python
    """
    allparts = []
    # I don't know how this works but it does, magic?
    while True:
        parts = os.path.split(path)
        if parts[0] == path:  # sentinel for absolute paths
            break
        elif parts[1] == path: # sentinel for relative paths
            allparts.insert(0, parts[1])
            break
        else:
            path = parts[0]
            allparts.insert(0, parts[1])
    return allparts


def extract_parameter_from_fname(fname, parameter, prefix):
    """
    Extract parameters from filename under the following circumstances:
    INPUT fname: res0.6_cluster.txt.gz or <parametervalue>_<analysis_type><file suffix>
    INPUT parameters: "res"
    OUTPUT: "0.6"
    Currently it doesn't work if the paramtervalue combo is next to the file suffix. (e.g. res0.6.txt.gz)
    """
    # first remove the sample prefix from the fname
    # so that any character strings that match the parameters are not picked up instead
    fname = re.sub("^" + prefix, "", fname)
    # split by underscores
    fname_split0 = [re.split('_', x) for x in splitall(fname)]  # split by path and _
    fname_split1 = list(chain.from_iterable(fname_split0))
    # extract the string that matches the requested parameter
    fname_grep = list(filter(lambda x: parameter in x, fname_split1))[0]
    # extract the specific value
    value = fname_grep.replace(parameter, "")
    # if the value is a number, convert to float or integer as appropriate
    if is_float_try(value):
        value = float(value)
        # convert integer values to integers e.g. 1.0 to 1 to match seurat columns
        if int(value) == float(value):
            value = int(value)
    return value

def intersect_obs_by_mod(mdata: MuData, mods:list =None):
    """
    Subset observations (samples or cells) in-place by intersect
    taking observations present only in modalities listed in mods, or all modalities if mods is None.
    Parameters
    ----------
    mdata: MuData
            MuData object
    """
    if mdata.isbacked:
        warnings.warn(
            "MuData object is backed. It might be required to re-read the object with `backed=False` to make the intersection work."
        )
    if mods is None:
        # we assume all modalities if mods is None
        mods = mdata.mod.keys()

    if len(mods) < 2:
        raise ValueError("cannot run intersect_obs on one modality")

    common_obs = reduce(np.intersect1d, [mdata[x].obs_names for x in mods])
    for mod in mdata.mod:
        mu.pp.filter_obs(mdata.mod[mod], common_obs)
    mdata.update_obs()
    return


def setdiff_obs_by_mod(mdata: MuData, x: str, y:str):
    """
    Subset observations (samples or cells) in-place by set difference
    taking observations present only in the mdata[x] and filtering mdata[y] to only contain the obs in x
    Use case example. You want to filter the rep modality to only include obsnames from rna, but you do not want
    to filter rna to only contain obs in the rep
    Parameters
    ----------
    mdata: MuData
            MuData object
    x:
    """
    if mdata.isbacked:
        warnings.warn(
            "MuData object is backed. It might be required to re-read the object with `backed=False` to make the intersection work."
        )
    
    # get the intersect between the two mods
    common_obs = reduce(np.intersect1d, [mdata[x].obs_names for x in [x,y]])
    # only filter the y mod
    mu.pp.filter_obs(mdata.mod[y], common_obs)
    mdata.update_obs()
    return


def downsample_adata(adata, nn, cat_col=None, seed=0):
    # if no category col, but adata in a list
    if cat_col is None:
        adatas = [adata]
    else:
        #just in case, remove any empty categories
        adata.obs[cat_col] = adata.obs[cat_col].cat.remove_unused_categories()
        # copy the obs for later
        # orig_obs = adata.obs.copy()
        adatas = [adata[adata.obs[cat_col].isin([clust])] for clust in adata.obs[cat_col].cat.categories]
    # do the subsample
    for dat in adatas:
        if dat.n_obs > nn:
            subsample(dat, n_obs=nn, random_state=seed)
    # the bit below is to be compatible with muon, otherwise we could just return adata_downsampled
    # get barcodes
    keep_bc = chain.from_iterable([list(x.obs_names) for x in adatas])
    adata = adata[list(keep_bc), ].copy()
    return(adata)


def downsample_mudata(mdata, nn, cat_col=None, mods: list = ['rna'], inplace=True, seed=0):
    # if inplace is False then we need to copy the object
    if inplace is False:
        mdata = mdata.copy()
    # do the downsample for 1 of the modalities.
    mdata.mod[mods[0]] = downsample_adata(mdata[mods[0]], cat_col=cat_col, nn=nn, seed=seed)
    if len(mods) > 1:
        intersect_obs_by_mod(mdata, mods=mods )
    return mdata




# def remove_double_TRB(filt_contig):
#     # remove clonotypes with more than one TRB
#     df = filt_contig[filt_contig['chain']=="TRB"].copy()
#     df['n_TRB'] = [x + 1 for x in list(df.groupby(['barcode']).cumcount())]
#     keep_bcs = list(df[df['n_TRB'] !=1].barcode.unique())
#     logging.info("removing %i annotations due to having more than 1 TRB chain" % len(keep_bcs))
#     filt_contig = filt_contig[~filt_contig["barcode"].isin(keep_bcs)]
#     return filt_contig, keep_bcs


# def define_clone_id(filt_contigs):
#     #filter out low confidence
#     df = filt_contigs[filt_contigs.is_cell & filt_contigs.high_confidence]
#     df = df[["barcode",'length', 'chain', 'v_gene', 'j_gene', 'c_gene', 'cdr3', 'cdr3_nt', "umis", "reads"]]
#     df['idx'] = [x + 1 for x in list(df.groupby(['barcode', 'chain']).cumcount())]
#     df['idx'] = df['idx'].astype("int64")
#     df['chain2'] = df['chain'] + df['idx'].astype("string")
#     df = df.pivot(index='barcode',columns='chain2')
#     df.columns = ['_'.join(col).strip() for col in df.columns.values]
#     df['clone_id'] = df.loc[:, df.columns.str.startswith('cdr3_TR')].fillna('').apply(lambda x: '-'.join(x),axis=1)
#     return df.loc[:,~df.columns.str.startswith("idx")]




def which_ind(bool_list):
    return list((compress(range(len(bool_list)), bool_list)))


def which_val(bool_list, val_list):
    return list((compress(val_list, bool_list)))


def check_for_bool(input_string):
    out = None
    if isinstance(input_string, str):
        if input_string == 'True':
            out = True
        elif input_string == 'False':
            out = False
        else:
            logging.info(input_string)
            raise TypeError("type error: please specify --%s as True or False" % input_string)
    elif isinstance(input_string, bool):
        return input_string
    else:
        # in this case the input is not a bool or a string
        raise TypeError("type error : please specify --%s as True or False (as a bool or a string)" %  input_string)
    return out

def test_file_or_value(test):
    try:
        float(test)
        out = "value"
    except:
        if os.path.isfile(test):
            out = "file"
        else:
            raise ValueError("no file or value found")
    return out

# this and the next function are very similar in that they essentialy are merge pandas functions
# either by index or by column.
def merge_with_adata_obs(adata, df, on_col="sample_id", inplace=False):
    """
    Merge a dataframe with adata.obs where the indexes do not match 
    between adata.obs and the df
    """
    if not isinstance(adata, AnnData):
        raise TypeError("adata is not an AnnData object")
    if not isinstance(df, pd.DataFrame):
        raise TypeError("df is not an Pandas dataframe")
    if on_col not in adata.obs.columns:
        raise KeyError("%s is not in adata.obs.columns" % on_col)
    if on_col not in df.columns:
        raise KeyError("%s is not in df.columns" % on_col)
    # preserve the index
    adata.obs['barcode_id'] = adata.obs.index
    new_df = adata.obs.merge(df, how="left", on=on_col)
    # reenter index
    new_df.index = new_df['barcode_id']
    new_df = new_df.drop(columns=['barcode_id'])
    new_df.index.name = None
    if inplace:
        adata.obs = new_df
    else:
        return new_df


def add_var_mtd(prot, df, left_on="adt_id", right_on="adt_id"):
    if not isinstance(prot, AnnData):
        raise TypeError("prot is not an AnnData object")
    if not isinstance(df, pd.DataFrame):
        raise TypeError("df is not an Pandas dataframe")
    # if the column we need is currently the index, make it a column
    if df.index.name == right_on:
        df[right_on] = df.index
        df.index.name=None
    # if the indexes of the adata and the df are equivalent go ahead with merge
    if len(set(prot.var[left_on]).intersection(set(df[right_on]))) == len(df[right_on]):
        # merge using indexes keeping only the vars from the orig prot
        prot_var = prot.var.merge(df, 
        left_on=left_on, right_on=right_on, 
        how="left")
        # check that the order is still ok (although the new df has no index, 
        # so we need to drop the index on the prev one.)
        if prot_var[left_on].equals(prot.var[left_on].reset_index(drop=True)):
            orig_index= prot.var.index
            prot.var = prot_var
            prot.var.index = orig_index
            logging.debug(prot.var.head())
        else:
            warnings.warn(UserWarning("prot.var is at risk of getting disrodered, not merging the df"))
    



def update_var_index(adata, new_index_col):
    adata.var = adata.var.reset_index().set_index(new_index_col)


def concat_adatas(adatas, batch_key, join_type="inner"):
    if isinstance(adatas, list) and len(adatas)!=1:
        # L.info("concatenating ...")
        # pull out the sample names in the correct order for the adaatas list
        batches = [x.obs[batch_key][0] for x in adatas]
        adata = adatas[0].concatenate(adatas[1:], 
                                        batch_key=batch_key, 
                                        batch_categories=batches, 
                                        join=join_type)
    elif isinstance(adatas, list) and len(adatas)==1:
        adata = adatas[0]
    else:
        raise TypeError("adatas is not a list")
    return(adata)


def concat_mdatas(mdata_list, batch_key, join_type="inner"):
    logging.debug(batch_key)
    if len(mdata_list)!=1:
        slots_filled = set(chain(*[list(m.mod.keys()) for m in mdata_list]))
        logging.debug(slots_filled)
        # need to concatenate each slot separately then recreate mdata
        concat_adatas={}
        for sf in slots_filled:
            logging.debug(sf)
            # extract if modality exists (it's possuible that not all samples have all modalities (thinking specifically of rep))
            adata_list = [x[sf]for x in mdata_list if sf in x.mod.keys()]
            [x.var_names_make_unique() for x in adata_list]
            batches = [x.obs[batch_key][0] for x in adata_list]
            concat_adatas[sf] = adata_list[0].concatenate(adata_list[1:],
                                    batch_key=batch_key,
                                    batch_categories=batches, 
                                    join=join_type)
            
        mdata = MuData(concat_adatas)
        # for sf in slots_filled:
        #     mdata[sf].obs[batch_key] ==concat_adatas[sf].obs[batch_key]
        
        # make sure sample id is in the top obs (updated to not put extra rep columns in top obs)
        # changed to deal with circumstances where a category is missing from one modality
        all_sample_id_df = mdata.obs.iloc[:,mdata.obs.columns.str.endswith(":sample_id")]
        # first it needs to not be a category for the merge.
        all_sample_id_df = all_sample_id_df.apply(lambda x: x.astype('object'))
        mdata.obs['sample_id']  = all_sample_id_df.fillna('').astype(str).apply(lambda x: ' '.join(set(' '.join(x).split())), axis=1)
        # then turn it back into a category
        mdata.obs['sample_id'] = mdata.obs['sample_id'].astype('category')
        return mdata
    else:
        return mdata_list[0]

def concat_adata_list(adata_list, use_muon, **kwargs):
    if use_muon:
        out = concat_mdatas(adata_list, **kwargs)
    else:
        out = concat_adatas(adata_list, **kwargs)
    return out

def dedup(value):
    words = set(value.split(' '))
    return ' '.join(words)


def intersection(lst1, lst2): 
    lst3 = [value for value in lst1 if value in lst2] 
    return lst3 


def remove_unused_categories(df: pd.DataFrame):
    if not isinstance(df, pd.DataFrame):
        raise TypeError('df is not a pandas dataframe')
    for c in df.columns:
        if pd.api.types.is_categorical_dtype(df[c]):
            df[c] = df[c].cat.remove_unused_categories()
    

def mu_get_obs(mdata, features=[],modalities=[], layers=None):
    """
    returns pandas dataframe of features, having searched all layers for said features
    """
    if modalities is None:
        # assume we want all modalities
        mods = mdata.mod.keys()
    else:
        mods = modalities
    if layers is None:
        # we assume we want the X layer
            mod_dict = dict(zip(mods, [None]*len(mods)))
    else:
        if len(mods) == len(layers):
            mod_dict = dict(zip(mods, layers))
        else:
            logging.info("make sure mods is the same length as layers")

    out = {}
    for mod, layer in mod_dict.items():
        out[mod] = obs_df(mdata[mod], keys=features, layer=layer)
    df = concat(out, axis=1)
    df.columns = ['-'.join(col).strip() for col in df.columns.values]
    return df
