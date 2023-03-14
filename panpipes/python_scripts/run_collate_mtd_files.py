import argparse
import numpy as np
import pandas as pd
import muon as mu
import os
import re
import sys
import logging
import yaml

L = logging.getLogger()
L.setLevel(logging.INFO)
log_handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s: %(levelname)s - %(message)s')
log_handler.setFormatter(formatter)
L.addHandler(log_handler)


parser = argparse.ArgumentParser()
parser.add_argument("--input_mudata")
parser.add_argument("--input_umap_files")
parser.add_argument('--rna_integration_col', default=None, help='the column in the adata.obs on which you have integrated')
parser.add_argument('--prot_integration_col', default=None, help='the column in the adata.obs on which you have integrated')
parser.add_argument('--atac_integration_col', default=None, help='the column in the adata.obs on which you have integrated')
parser.add_argument('--multimodal_integration_col', default=None, help='the column in the adata.obs on which you have integrated')
parser.add_argument("--output_cell_metadata_csv")
parser.add_argument("--output_combined_umaps_tsv")
parser.add_argument("--output_batch_yml")
args, opt = parser.parse_known_args()
L.info(args)
L.info("reading in files")

cell_meta_df = mu.read(args.input_mudata).obs
mtd_columns = cell_meta_df.columns.to_list()

# get all the batch columns 
# load up all the batch columns
batch_dict = dict(rna=args.rna_integration_col, 
                  prot=args.prot_integration_col, 
                  atac=args.atac_integration_col, 
                  multimodal=args.multimodal_integration_col)
batch_dict = {k:v.split(',') for k,v in batch_dict.items() if v is not None}

# format the column names corrrectly, and add any merge columns to the dataset
for k in batch_dict.keys():
    v=batch_dict[k]
    v =[k + ":" +v1 if k!="multimodal" else v1 for v1 in v  ]
    batch_dict[k] = v
    if len(v) > 1:
        mod_meta_df= cell_meta_df[v].dropna()
        cell_meta_df[str(k)+ ":bc_batch"] = mod_meta_df.apply(lambda x: '|'.join(x), axis=1)
        batch_dict[k].append(k+ ":bc_batch")

cell_meta_df.to_csv(args.output_cell_metadata_csv)

with open(args.output_batch_yml, 'w') as outfile:
    yaml.dump(batch_dict, outfile, default_flow_style=False)

L.info("reading in all umaps")
# load files
input_files = args.input_umap_files.split(",")
umaps_list = [pd.read_csv(x, index_col=0) for x in input_files]
umaps_index = [re.sub("umap_|.csv", "", os.path.basename(x)) for x in input_files]
umaps_df = pd.concat(umaps_list, axis=0, keys=umaps_index).reset_index(level=0)
umaps_df.columns = ["method", "umap_1", "umap_2"]
umaps_df = umaps_df.merge(cell_meta_df, left_index=True, right_index=True)

# reorder methods so None is at the most left on plot
umaps_df['mod'] = [x.split('_')[0] for x in umaps_df.method]
umaps_df['method'] = [x.split('_')[1] for x in umaps_df.method]

# tidy up umaps df before writing out
def move_column_inplace(df, col, pos):
    col = df.pop(col)
    df.insert(pos, col.name, col)
move_column_inplace(umaps_df, "mod", 0)
move_column_inplace(umaps_df, "method", 1)
umaps_df = umaps_df.sort_values(by=['mod','method'])

# save umap df to file
L.info("save umap df to file")
umaps_df = umaps_df.iloc[:,0:4]
umaps_df['cellbarcode'] = umaps_df.index
umaps_df.to_csv(args.output_combined_umaps_tsv, sep="\t")

L.info("done")

