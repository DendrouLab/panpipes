import muon as mu
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import os
import itertools
import re

from panpipes.funcs.io import read_yaml
from panpipes.funcs.plotting import get_layer, plot_scatters


import sys
import logging
L = logging.getLogger()
L.setLevel(logging.INFO)
log_handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s: %(levelname)s - %(message)s')
log_handler.setFormatter(formatter)
L.addHandler(log_handler)
L.debug("test logging message")


# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('--mdata_object',
                    default='',
                    help='')
parser.add_argument('--layers_dict',
                    default='',
                    help='')
parser.add_argument('--scatters_csv',
                    default='',
                    help='')

parser.set_defaults(verbose=True)
args, opt = parser.parse_known_args()

L.info(args)

sc.settings.figdir  = "./scatters/"
sc.set_figure_params(fontsize=12)


mdata = mu.read(args.mdata_object)

layers = read_yaml(args.layers_dict)


df = pd.read_csv(args.scatters_csv)
L.debug(df)

layers = df.applymap(lambda x: get_layer(x,  mdata=mdata, layers=layers))
L.info(layers)

for idx in range(df.shape[0]):
    features_list = df.iloc[idx,:].tolist()
    # convert NAs and X to None
    features_list = [x  if pd.notna(x) else None for x in features_list]
    
    layers_list = layers.iloc[idx,:].tolist()
    # exclude any Nones to get the layers list to be the same length as features
    layers_list = [x for idx, x in enumerate(layers_list) if features_list[idx] is not None]
    
    fname = os.path.join("scatters", "_".join([x for x in features_list if x is not None]) + ".png")
    plot_scatters(mdata, features_list, layers_list)
    # plt.legend(bbox_to_anchor=(1.04,0.5), loc="center left", borderaxespad=0)
    plt.subplots_adjust(right=0.85)
    plt.savefig(fname, bbox_inches="tight")

L.info("Done")

