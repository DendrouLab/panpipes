import muon as mu
import scanpy as sc
import argparse
import pandas as pd
import os
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
                    default="sc_preprocess.h5ad", help="file name, format: .h5ad")
parser.add_argument('--modality',
                    default=None,
                    help='')
parser.add_argument("--outfile", 
                    default="clusters.txt", help="file name, format: .txt")
parser.add_argument("--resolution",
                    default=0.5, help="no. neighbours parameters for sc.pp.neighbors()")
parser.add_argument("--algorithm", 
                    default="leiden", help="algortihm choice from louvain and leiden")
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

uns_key=args.neighbors_key
# check sc.pp.neihgbours has been run
if uns_key not in adata.uns.keys():
    # sys.exit("Error: sc.pp.neighbours has not been run on this object")
    sc.pp.neighbors(adata)
    uns_key="neighbors"


# run command
if args.algorithm == "louvain":
    sc.tl.louvain(adata, resolution=float(args.resolution), key_added='clusters', neighbors_key=uns_key)
elif args.algorithm == "leiden":
    sc.tl.leiden(adata, resolution=float(args.resolution), key_added='clusters', neighbors_key=uns_key)
else:
    sys.exit("algorithm not found: please specify 'louvain' or 'leiden'")

#mdata.update()
## write out clusters as text file
clusters = pd.DataFrame(adata.obs['clusters'])
clusters.to_csv(args.outfile, sep='\t')

tmp = clusters['clusters'].value_counts().to_frame("cell_num").reset_index().rename(columns={"index":"cluster"})
tmp.to_csv(os.path.dirname(args.outfile) + "/cellnum_per_cluster.csv")

if "sample_id" in adata.obs.columns:
    ccounts = pd.DataFrame(adata.obs[['sample_id', 'clusters']])
    tmp = ccounts.value_counts().to_frame("cell_num").reset_index()
    tmp.to_csv(os.path.dirname(args.outfile) + "/cellnum_per_sample_id_per_cluster.csv")

L.info("Done")

