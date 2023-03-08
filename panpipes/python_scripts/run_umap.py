"""
# reads in an anndata object with neighbours already computed, runs scanpy umap
# save umap coordinates (for plotting elsewhere
# CRG 2020-06-12
"""
import scanpy as sc
import numpy as np
import pandas as pd
import argparse
import muon as mu
from anndata import AnnData

import sys
import logging
L = logging.getLogger()
L.setLevel(logging.INFO)
log_handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s: %(levelname)s - %(message)s')
log_handler.setFormatter(formatter)
L.addHandler(log_handler)

parser = argparse.ArgumentParser()
parser.add_argument("--infile",
                    default="mdata.h5mu",
                    help="file name, format: .h5mu")
parser.add_argument('--modality',
                    default=None,
                    help='')
parser.add_argument("--outfile", 
                    default="umap.txt.gz", 
                    help="file name, format: .h5ad")
parser.add_argument("--min_dist", 
                    default=0.1, 
                    help="no. neighbours parameters for sc.pp.neighbors()")
parser.add_argument("--neighbors_key", 
                    default="neighbors", help="algortihm choice from louvain and leiden")

args, opt = parser.parse_known_args()
L.info(args)

# read data
mdata = mu.read(args.infile)
if type(mdata) is AnnData:
    adata = mdata
elif args.modality is not None:
    adata = mdata[args.modality]
else:
    adata = mdata
    

# set seed
# seed = int(200612)

uns_key=args.neighbors_key
# check sc.pp.neihgbours has been run
if uns_key not in adata.uns.keys():
    # sys.exit("Error: sc.pp.neighbours has not been run on this object")
    L.warning("running neighbors with default parameters since no neighbors graph found in this data object")
    sc.pp.neighbors(adata)
    uns_key="neighbors"

# code to run multiple umaps in one script, currently not used
# manipulate min dist to list of floats
# if type(args.min_dist) is str:
#     min_dist = [float(md) for md in args.min_dist.split(',')]
# elif type(args.min_dist) is float:
#     min_dist = [args.min_dist]


# what parameters?
if uns_key =="wnn":
    mu.tl.umap(adata, min_dist=float(args.min_dist), neighbors_key=uns_key)
else:
    sc.tl.umap(adata, min_dist=float(args.min_dist), neighbors_key=uns_key)

# extract umap coordinates for plotting (in R??)
umap_coords = pd.DataFrame(adata.obsm['X_umap'])

# add in the rownames 
umap_coords.index = adata.obs_names

# save coordinates to file
# (note this saves values values up to 6 significant figures, because why save 20 for a plot
umap_coords.to_csv(args.outfile, sep = '\t')
