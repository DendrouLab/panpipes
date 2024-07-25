## Filtering mudata based on n_genes, percent mito and min cells/
## originally written by Tom Thomas (https://github.com/tomthomas3000/TAURUS)
## adapted and augmented for this pipeline by Charlotte Rich-Griffin 2020-09-30

import scanpy as sc
import argparse
import pandas as pd
import re
import muon as mu
from anndata import AnnData
from muon import MuData
import yaml
import os

# import scpipelines.funcs as scp
from panpipes.funcs.processing import intersect_obs_by_mod, remove_unused_categories
from panpipes.funcs.io import write_obs, read_yaml, dictionary_stripper
import collections.abc
import sys
import logging
L = logging.getLogger()
L.setLevel(logging.INFO)
log_handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s: %(levelname)s - %(message)s')
log_handler.setFormatter(formatter)
L.addHandler(log_handler)


sc.settings.verbosity = 2  # verbosity: errors (0), warnings (1), info (2), hints (3)
# comment this because of numba issues
# sc.logging.L.info_versions()



def map_nested_dicts_remove_none(ob):
    if isinstance(ob, collections.abc.Mapping):
        return {k: map_nested_dicts_remove_none(v) for k, v in ob.items() if v is not None}
    else:
        return ob


def test_matching_df_ignore_cat(new_df, old_df):
    remove_unused_categories(old_df)
    old_df = old_df.loc[new_df.index, :]
    remove_unused_categories(old_df)
    return new_df.equals(new_df)

# parse arguments
parser = argparse.ArgumentParser()

parser.add_argument('--input_mudata',
                    default='gut_minus1_amp.h5ad',
                    help='')
parser.add_argument('--output_mudata',
                    default='',
                    help='')
parser.add_argument('--filter_dict',
                    default='',
                    help='this is pull')
# cross modalities args
parser.add_argument('--keep_barcodes', default=None,
                    help='1 column list of barcodes to keep, note that they should match the mudata input, this filtering happens first')
parser.add_argument('--intersect_mods', default=None,
                    help='comma separated string of modalities we want to intersect_obs on')                 


# load options
parser.set_defaults(verbose=True)
args, opt = parser.parse_known_args()
L.info("Running with params: %s", args)
 
# load the filtering dictionary 
filter_dict = args.filter_dict
if isinstance(args.filter_dict, dict):
    filter_dict = args.filter_dict
else:
    filter_dict = read_yaml(args.filter_dict) 

# remove Nones for each level of filter dict
filter_dict = map_nested_dicts_remove_none(filter_dict)
filter_dict = dictionary_stripper(filter_dict)
L.info("Filter dictionary:\n %s" %filter_dict)

# load mudata
L.info("Reading in MuData from '%s'" % args.input_mudata)
mdata = mu.read(args.input_mudata)

if isinstance(mdata, AnnData):
    raise TypeError("Input '%s' should be of MuData format, not Anndata"  % args.input_mudata)

orig_obs = mdata.obs.copy()

L.info("Before filtering: "+ str(mdata.n_obs) + " cells and " + str(mdata.n_vars) + " features")

# filter based on provided barcodes -----
if args.keep_barcodes is not None:
    if os.path.exists(args.keep_barcodes):
        L.info("Reading in keep_barcodes file from '%s'" % args.keep_barcodes)
        keep_bc = pd.read_csv(args.keep_barcodes,header=None)
    else:
        L.error("The path of the keep_barcodes file '%s' could not be found" % args.keep_barcodes)
        sys.exit("The path of the keep_barcodes file '%s' could not be found" %  args.keep_barcodes)

    L.info("Filtering all modalities by keep_barcodes file")
    mdata = mdata[mdata.obs_names.isin(keep_bc[0]),:].copy()
    remove_unused_categories(mdata.obs)
    mdata.update()
    L.info("Remaining cells %d" % mdata.n_obs)

# L.debug(mdata.obs['sample_id'].value_counts())

# filter more than
if filter_dict['run']:
    # this will go through the modalities one at a time,
    # then the categories max, min and bool
    for mod in mdata.mod.keys():
        L.info("Filtering modality '%s'" % mod)
        if mod in filter_dict.keys():
            for marg in filter_dict[mod].keys():
                if marg == "obs":
                    if "max" in filter_dict[mod][marg].keys():
                        for col, n in filter_dict[mod][marg]['max'].items():
                            L.info("Filtering cells of modality '%s' by '%s' in obs to less than %s" % (mod, col, n))
                            mu.pp.filter_obs(mdata.mod[mod], col, lambda x: x <= n)
                            L.info("Remaining cells %d" % mdata[mod].n_obs)
                    if "min" in filter_dict[mod][marg].keys():
                        for col, n in filter_dict[mod][marg]['min'].items():
                            L.info("Filtering cells of modality '%s' by '%s' in obs to more than %s" % (mod, col, n))
                            mu.pp.filter_obs(mdata.mod[mod], col, lambda x: x >= n)
                            L.info("Remaining cells %d" % mdata[mod].n_obs)
                    if "bool" in filter_dict[mod][marg].keys():
                        for col, n in filter_dict[mod][marg]['bool'].items():
                            L.info("Filtering cells of modality '%s' by '%s' in obs marked %s" % (mod, col, n))
                            mu.pp.filter_obs(mdata.mod[mod], col, lambda x: x == n)
                            L.info("Remaining cells %d" % mdata[mod].n_obs)
                if marg == "var":
                    if "max" in filter_dict[mod][marg].keys():
                        for col, n in filter_dict[mod][marg]['max'].items():
                            L.info("Filtering features of modality '%s' by '%s' in .var to less than %s" % (mod, col, n))
                            mu.pp.filter_var(mdata.mod[mod], col, lambda x: x <= n)
                            L.info("Remaining features: %d" % mdata[mod].n_vars)

                    if "min" in filter_dict[mod][marg].keys():
                        for col, n in filter_dict[mod][marg]['min'].items():
                            L.info("Filtering features of modality '%s' by '%s' in .var to more than %s" % (mod, col, n))
                            mu.pp.filter_var(mdata.mod[mod], col, lambda x: x >= n)
                            L.info("Remaining features: %d" % mdata[mod].n_vars)

                    if "bool" in filter_dict[mod][marg].keys():
                        for col, n in filter_dict[mod][marg]['bool'].items():
                            L.info("Filtering features of modality '%s' by '%s' in .var marked %s" % (mod, col, n))
                            mu.pp.filter_var(mdata.mod[mod], col, lambda x: x == n)
                            L.info("Remaining features: %d" % mdata[mod].n_vars)
                            

mdata.update()

# intersect specific modalities
# we don't want to just use mu.pp.intersect_obs, 
# because we don't want to necessarily filter by repertoire
if args.intersect_mods is not None:
    intersect_mods = args.intersect_mods.split(',')
    intersect_mods= [a.strip() for a in intersect_mods]
    L.info("Intersecting barcodes of modalities %s" % args.intersect_mods)
    intersect_obs_by_mod(mdata, mods = intersect_mods)

L.info("After filtering: "+ str(mdata.n_obs) + " cells and " + str(mdata.n_vars) + " features across all modalities")

remove_unused_categories(mdata.obs)

# run quick test before saving out.
assert test_matching_df_ignore_cat(mdata.obs, orig_obs)  

# write out obs
output_prefix = re.sub(".h5mu", "", args.output_mudata)

L.info("Saving updated obs in a metadata tsv file to '" + output_prefix + "_filtered_cell_metadata.tsv'")
write_obs(mdata, output_prefix=output_prefix, output_suffix="_filtered_cell_metadata.tsv")

# write out the per sample_id cell numbers 
cell_counts_dict={}
for mm in mdata.mod.keys():
    cell_counts_dict[mm] = mdata[mm].obs.sample_id.value_counts().to_frame('n_cells')

cell_counts = pd.concat(cell_counts_dict).reset_index().rename(
    columns={"level_0": "modality", "level_1": "sample_id"})

L.info("cell_counts\n%s" %cell_counts)
L.info("Saving cell counts in a metadata csv file to '" + output_prefix + "_cell_counts.csv'")
cell_counts.to_csv(output_prefix + "_cell_counts.csv", index=None)


mdata.update()
L.info("Saving updated MuData to '%s'" % args.output_mudata)
mdata.write(args.output_mudata)

L.info("Done")

