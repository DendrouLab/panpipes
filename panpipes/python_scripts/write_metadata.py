import argparse
import pandas as pd
import scanpy as sc

from muon import read as muread
from panpipes.funcs.processing import check_for_bool
from panpipes.funcs.io import read_anndata

import sys
import logging
L = logging.getLogger()
L.setLevel(logging.INFO)
log_handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s: %(levelname)s - %(message)s')
log_handler.setFormatter(formatter)
L.addHandler(log_handler)

parser = argparse.ArgumentParser()
parser.add_argument('--infile',
                    default="./mdata.h5mu",
                    help="mudata object ")
parser.add_argument('--outfile', default='sc_preprocess.txt',
                    help="path for text file output")
args = parser.parse_args()

args, opt = parser.parse_known_args()

# read file
mdata = muread(args.infile)
mdata.obs.to_csv(args.outfile, sep='\t')