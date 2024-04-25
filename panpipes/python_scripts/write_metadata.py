import argparse
import pandas as pd
import scanpy as sc

from muon import read 

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
L.info("Running with params: %s", args)

# read file
L.info("Reading in MuData from '%s'" % args.infile)
mdata = read(args.infile)

L.info("Writing obs to '%s'" % args.outfile)
mdata.obs.to_csv(args.outfile, sep='\t')