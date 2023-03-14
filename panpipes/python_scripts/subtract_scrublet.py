
import scanpy as sc
import argparse
import pandas as pd
import os


from panpipes.funcs.processing import test_file_or_value, merge_with_adata_obs, check_for_bool
from panpipes.funcs.io import read_anndata, write_anndata 

import sys
import logging
L = logging.getLogger()
L.setLevel(logging.INFO)
log_handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s: %(levelname)s - %(message)s')
log_handler.setFormatter(formatter)
L.addHandler(log_handler)

sc.settings.verbosity = 0  # verbosity: errors (0), warnings (1), info (2), hints (3)
# comment this because of numba issues
# sc.logging.print_versions()

# parse arguments
parser = argparse.ArgumentParser()

parser.add_argument('--input_anndata',
                    default='data/anndata-n5000-filt.h5ad',
                    help='')
parser.add_argument('--output_anndata',
                    default='data/anndata-n5000-filt.h5ad',
                    help='')
parser.add_argument('--use_muon',
                    default=False,type=check_for_bool,
                    help='')
parser.add_argument('--output_prefix',
                    default='',
                    help='prefix to prepend to saved files. here just used for metadata file')
parser.add_argument('--scrublet_cutoff', default=None,
                    help='this is either a value or a file')
parser.add_argument('--batch_col', default="sample_id",
                    help='')

L.warning("Running subtract scrublet")

parser.set_defaults(verbose=True)
args, opt = parser.parse_known_args()

use_muon = args.use_muon

# load data
# if muon is True then rna only is loaded
# if muon is False then the modality argument is ignored and the normal anndata h5ad format is loaded
adata = read_anndata(args.input_anndata, use_muon=use_muon, modality="rna")

# test what kind of input we are dealing with, value or file?
run = test_file_or_value(args.scrublet_cutoff)


# adata=sc.read_h5ad("./data/run_integration/taurus_test_filt.h5ad")

# run = test_file_or_value("./data/run_integration/scrublet_cutoffs.csv")
# subtract scrublet
if run == 'value':
    # if args.scrublet_cutoff is not None:
    # if it is a value
    adata = adata[adata.obs['doublet_scores'] < float(args.scrublet_cutoff), :]
elif run == 'file':
    # if it is a file
    scr_df = pd.read_csv(args.scrublet_cutoff, names=["sample_id", "scrublet_cutoff"])
    adata.obs = merge_with_adata_obs(adata.obs, scr_df, on_col="sample_id")
    adata = adata[adata.obs['doublet_scores'] < adata.obs['scrublet_cutoff'],]
else:
    pass

L.info("saving anndata and obs in a metadata tsv file")
metafile = adata.obs
metafile["cellbarcode"] = adata.obs.index

savename= args.output_prefix + "_filtered_cell_metadata.tsv"
metafile.to_csv(savename, sep='\t', index=True)  

write_anndata(adata, args.input_anndata, use_muon=use_muon, modality="rna")

L.info("done")