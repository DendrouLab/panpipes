import muon as mu
import pandas as pd 
import numpy as np
import argparse
import sys
import logging
import re
from muon import MuData
from anndata import AnnData


L = logging.getLogger()
L.setLevel(logging.INFO)
log_handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s: %(levelname)s - %(message)s')
log_handler.setFormatter(formatter)
L.addHandler(log_handler)

# parse argumetns
parser = argparse.ArgumentParser()
parser.add_argument("--input_mudata",
                    default="mdata.h5mu",
                    help="file name, format: .h5mu")
parser.add_argument("--clusters_files_csv",
                    default=None,
                    help="comma separated table of mod,compiled clusters dataframes")
parser.add_argument("--umap_files_csv",
                    default=None,
                    help = "comma separated  table of mod,txt file containing umap coordinates")
parser.add_argument("--output_mudata",
                    default="mdata.h5mu",
                    help="file name, format: .h5mu")
args, opt = parser.parse_known_args()
L.info(args)
mdata = mu.read(args.input_mudata)

L.info("loading clusters")
cf = pd.read_csv(args.clusters_files_csv)

if isinstance(mdata, MuData):
    pass
elif isinstance(mdata, AnnData):
    mds = cf['mod'].unique().tolist()
    if len(mds)>1:
        sys.exit("You have clustered multiple modalities but are providing only a unimodal anndata")
    else:
        L.warn("found one modality, converting to mudata: %s " % mds[0] )    
        tmp = MuData({mds[0]:mdata})
        del mdata
        mdata = tmp
        del tmp

# add in the clusters



for i in range(cf.shape[0]):
    cf_df = pd.read_csv(cf['fpath'][i], sep='\t', index_col=0) 
    cf_df['clusters'] = cf_df['clusters'].astype('str').astype('category')
    cf_df = cf_df.rename(columns={"clusters":cf['new_key'][i]})
    
    if cf['mod'][i] != "multimodal":
        mdata[cf['mod'][i]].obs = mdata[cf['mod'][i]].obs.merge(cf_df, left_index=True, right_index=True)
    else:
        mdata.obs = mdata.obs.merge(cf_df, left_index=True, right_index=True)

uf = pd.read_csv(args.umap_files_csv)


for i in range(uf.shape[0]):
    uf_df = pd.read_csv(uf['fpath'][i], sep='\t', index_col=0) 
    mod = uf['mod'][i]
    new_key = uf['new_key'][i]
    if uf['mod'][i] != "multimodal":
        if all(mdata[mod].obs_names == uf_df.index):
            mdata[mod].obsm[new_key] =  uf_df.to_numpy()
        else:
            L.warn("cannot integrate %s into mdata as obs_names mismatch" % uf.iloc[i,:] )
    else:
        # check the observations are the same
        if set(mdata.obs_names).difference(uf_df.index) == set():
            # put the observations in the same order
            uf_df = uf_df.loc[mdata.obs_names,:]
            mdata.obsm[new_key] =  uf_df.to_numpy()
        else:
            L.warning("cannot integrate %s into mdata as obs_names mismatch" % uf.iloc[i,:] )

mdata.write(args.output_mudata)
mdata.obs.to_csv(re.sub(".h5mu", "_cell_metdata.tsv", args.output_mudata), sep='\t')

L.info("done")
