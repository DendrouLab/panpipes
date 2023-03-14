'''
scanpy QC script GEX
order of QC:
- GEX
- ADT
- Repertoire
- ATAC
'''
import sys
import logging
L = logging.getLogger()
L.setLevel(logging.INFO)
log_handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s: %(levelname)s - %(message)s')
log_handler.setFormatter(formatter)
L.addHandler(log_handler)

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import argparse
import pandas as pd
import muon as mu
import scirpy as ir
import matplotlib.pyplot as plt
plt.ioff()
from panpipes.funcs.io import write_obs
import yaml

parser = argparse.ArgumentParser()
# required option
parser.add_argument("--sampleprefix",
                    default="",
                    help="prefix to prepend when saving the metadata file")
parser.add_argument("--input_mudata",
                    default="mdata_unfilt.h5mu",
                    help="")
parser.add_argument("--output_mudata",
                    default="mdata_unfilt.h5mu",
                    help="")
parser.add_argument("--figdir",
                    default="./figures/",
                    help="path to save the figures to")

parser.add_argument("--distance_metrics",
                    default=None, type=yaml.safe_load,
                    help="what metrics to calculate sequence distance metric? \
                    any arguments from scirpy.pp.ir_dist as a dict in format'{arg: value}' ")
parser.add_argument("--clonotype_metrics",
                    default=None,type=yaml.safe_load,
                    help="what metrics to calculate sequence distance metric? \
                    any arguments from scirpy.pp.define_clonotypes as a dict in format'{arg: value}' ")
args, opt = parser.parse_known_args()
L.info(args)
mdata = mu.read(args.input_mudata)
rep = mdata['rep']
# chain qc
ir.tl.chain_qc(rep)

# remove nones, so defaults are used
if args.distance_metrics is not None:
    # remove nones, so defaults are used
    dist_args = {k: v for k, v in args.distance_metrics.items() if v != "None"}
    dist_args = {k: v for k, v in dist_args.items() if v}
else:
    dist_args = {}

L.info("distance args:")
L.info(dist_args)

ir.pp.ir_dist(rep, **dist_args)

# pull function arguments from args. 
if args.clonotype_metrics is not None:
    # remove nones, so defaults are used
    clonotype_args = {k: v for k, v in args.clonotype_metrics.items() if v != "None"}
    clonotype_args = {k: v for k, v in clonotype_args.items() if v}
else:
    clonotype_args={}
# define clonotypes
L.info("clonotypes args:")
L.info(clonotype_args)

ir.tl.define_clonotypes(rep, **clonotype_args)
ir.tl.clonal_expansion(rep)

    
category_cols = ['has_ir']
for cc in category_cols:
    if pd.api.types.infer_dtype(rep.obs[cc]) != "categorical":
        rep.obs[cc] = rep.obs[cc].astype('string').astype('category')


if "BCR" in rep.obs['receptor_type'].values:
    bcr = rep[rep.obs.receptor_type =="BCR", :].copy()
else:
    bcr=None
    
if "TCR" in rep.obs.receptor_type.values:
    tcr = rep[rep.obs.receptor_type =="TCR", :].copy()
else:
    tcr=None

# plot group abundances
fig, ax = plt.subplots()
ax = ir.pl.group_abundance(rep, groupby="receptor_type", ax=ax)
fig.savefig(args.figdir +"/barplot_group_abundance_receptor_type.png", bbox_inches="tight")

if bcr is not None:
    fig, ax = plt.subplots()
    ax = ir.pl.group_abundance(bcr, groupby="chain_pairing", ax=ax)
    fig.tight_layout()
    fig.savefig(args.figdir +"/barplot_group_abundance_bcr_receptor_subtype.png", bbox_inches="tight")

if tcr is not None:
    fig, ax = plt.subplots()
    ax = ir.pl.group_abundance(tcr, groupby="chain_pairing",  ax=ax)
    fig.tight_layout()
    fig.savefig(args.figdir +"/barplot_group_abundance_tcr_receptor_subtype.png", bbox_inches="tight")

# clonal expansion

clone_counts = rep.obs.groupby("receptor_type").clone_id.value_counts().to_frame("cell_counts").reset_index()

if bcr is not None:
    ax = ir.pl.clonal_expansion(bcr, groupby="sample_id")
    ax.set_title("bcr clonal expansion")
    plt.savefig(args.figdir +"/barplot_clonal_expansion_bcr.png", bbox_inches="tight" )

if tcr is not None:
    ax = ir.pl.clonal_expansion(tcr, groupby="sample_id")
    ax.set_title("tcr clonal expansion")
    plt.savefig(args.figdir +"/barplot_clonal_expansion_tcr.png", bbox_inches="tight")


L.info("saving anndata and obs in a metadata tsv file")
mdata.update()
write_obs(mdata, output_prefix=args.sampleprefix, 
    output_suffix="_cell_metadata.tsv")
mdata.write(args.output_mudata)

L.info("Done")

