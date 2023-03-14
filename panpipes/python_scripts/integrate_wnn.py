import multiprocessing 

import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import os
import muon as mu
from cgatcore import pipeline as P


import sys
import logging
L = logging.getLogger()
L.setLevel(logging.INFO)
log_handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s: %(levelname)s - %(message)s')
log_handler.setFormatter(formatter)
L.addHandler(log_handler)

# load arguments 
parser = argparse.ArgumentParser()
parser.add_argument('--input_mudata',
                    default='mudata_scaled.h5ad',
                    help='')
parser.add_argument('--use_muon',
                    default=False,
                    help='')
parser.add_argument('--output_csv', default='batch_correction/umap_bc_totalvi.csv',
                    help='')
parser.add_argument('--integration_col', default='batch',
                    help='')
parser.add_argument('--figdir', default='./figures',
                    help='')
parser.add_argument('--neighbors_n_pcs',
                    help="n_pcs", default=50)
parser.add_argument('--neighbors_k', default=30,
                    help="neighbors k")
parser.add_argument('--neighbors_metric',default="euclidean",
                    help="neighbor metric, e.g. euclidean or cosine")

args, opt = parser.parse_known_args()

# scanpy settings
sc.set_figure_params(facecolor="white")
sc.settings.autoshow = False
sc.settings.figdir = args.figdir


# load parameters
threads_available = multiprocessing.cpu_count()
params = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "../pipeline.yml",
     "pipeline.yml"])


mdata = mu.read(args.input_mudata)
# we assume here there is a mudata object that already has a neighbours for each modality


mu.pp.neighbors(mdata)
mu.tl.umap(mdata)


L.info("Done")

