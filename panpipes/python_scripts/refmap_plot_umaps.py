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
parser.add_argument("--output_combined_umaps_tsv")
args, opt = parser.parse_known_args()
L.info(args)
L.info("reading in files")

# cell_meta_df = mu.read(args.input_mudata+'/rna').obs
# mtd_columns = cell_meta_df.columns.to_list()



L.info("reading in all umaps")
# load files
input_files = args.input_umap_files.split(",")
umaps_list = [pd.read_csv(x, index_col=0) for x in input_files]
umaps_index = [re.sub("umap_|.csv", "", os.path.basename(x)) for x in input_files]
umaps_df = pd.concat(umaps_list, axis=0, keys=umaps_index).reset_index(level=0)
umaps_df.columns = ["method", "umap_1", "umap_2"]
# umaps_df = umaps_df.merge(cell_meta_df, left_index=True, right_index=True)

# reorder methods so None is at the most left on plot
umaps_df['method'] = [x.split('_')[0] for x in umaps_df.method]

umaps_df = umaps_df.sort_values(by='method')

# save umap df to file
L.info("save umap df to file")
umaps_df.iloc[:,0:4].to_csv(args.output_combined_umaps_tsv, sep="\t")

L.info("done")

