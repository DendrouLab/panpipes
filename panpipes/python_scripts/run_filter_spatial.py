import scanpy as sc
import argparse
import pandas as pd
import re
import muon as mu
from anndata import AnnData
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


# load options
parser.set_defaults(verbose=True)
args, opt = parser.parse_known_args()
L.info('running with args:')
L.info(args)
 
# load the filtering dictionary 
filter_dict = args.filter_dict
if isinstance(args.filter_dict, dict):
    filter_dict = args.filter_dict
else:
    filter_dict = read_yaml(args.filter_dict) 

# remove Nones for each level of filter dict
filter_dict = map_nested_dicts_remove_none(filter_dict)
filter_dict = dictionary_stripper(filter_dict)
L.info("filter dictionary:\n %s" %filter_dict)

# load mudata
mdata = mu.read(args.input_mudata)

if isinstance(mdata, AnnData):
    raise TypeError("input_mudata should be of MuData format, not Anndata")

orig_obs = mdata.obs.copy()

L.info("Before filtering %d" % mdata.n_obs)
# filter based on provided barcodes -----
if args.keep_barcodes is not None:
    L.info("filtering all modalities by keep_barcodes file")
    keep_bc = pd.read_csv(args.keep_barcodes,header=None)
    mdata = mdata[mdata.obs_names.isin(keep_bc[0]),:].copy()
    remove_unused_categories(mdata.obs)
    mdata.update()
    L.info("Remaining cells %d" % mdata.n_obs)



# filter more than
if filter_dict['run']:
    # this will go through the modalities one at a time,
    # then the categories max, min and bool
    for mod in mdata.mod.keys():
        L.info(mod)
        if mod in filter_dict.keys():
            for marg in filter_dict[mod].keys():
                if marg == "obs":
                    if "max" in filter_dict[mod][marg].keys():
                        for col, n in filter_dict[mod][marg]['max'].items():
                            L.info("filter %s %s to less than %s" % (mod, col, n))
                            mu.pp.filter_obs(mdata.mod[mod], col, lambda x: x <= n)
                            L.info("Remaining cells %d" % mdata[mod].n_obs)
                            L.info("Remaining cells %d" % mdata[mod].shape[0])
                    if "min" in filter_dict[mod][marg].keys():
                        for col, n in filter_dict[mod][marg]['min'].items():
                            L.info("filter %s %s to more than %s" % (mod, col, n))
                            mu.pp.filter_obs(mdata.mod[mod], col, lambda x: x >= n)
                            L.info("Remaining cells %d" % mdata[mod].n_obs)
                    if "bool" in filter_dict[mod][marg].keys():
                        for col, n in filter_dict[mod][marg]['bool'].items():
                            L.info("filter %s %s marked %s" % (mod, col, n))
                            mu.pp.filter_obs(mdata.mod[mod], col, lambda x: x == n)
                            L.info("Remaining cells %d" % mdata[mod].n_obs)
                if marg == "var":
                    if "max" in filter_dict[mod][marg].keys():
                        for col, n in filter_dict[mod][marg]['max'].items():
                            L.info("filter %s %s to less than %s" % (mod, col, n))
                            mu.pp.filter_var(mdata.mod[mod], col, lambda x: x <= n)
                            L.info("Remaining cells %d" % mdata[mod].n_obs)

                    if "min" in filter_dict[mod][marg].keys():
                        for col, n in filter_dict[mod][marg]['min'].items():
                            L.info("filter %s %s to more than %s" % (mod, col, n))
                            mu.pp.filter_var(mdata.mod[mod], col, lambda x: x >= n)
                            L.info("Remaining cells %d" % mdata[mod].n_obs)

                    if "bool" in filter_dict[mod][marg].keys():
                        for col, n in filter_dict[mod][marg]['bool'].items():
                            L.info("filter %s %s marked %s" % (mod, col, n))
                            mu.pp.filter_var(mdata.mod[mod], col, lambda x: x == n)
                            L.info("Remaining cells %d" % mdata[mod].n_obs)
                            

mdata.update()

# intersect specific modalities
# we don't want to just use mu.pp.intersect_obs, 
# because we don't want to necessarily filter by repertoire
#if args.intersect_mods is not None:
#    intersect_mods = args.intersect_mods.split(',')
#    intersect_mods= [a.strip() for a in intersect_mods]
#    L.info("intersecting barcodes in %s" % args.intersect_mods)
#    intersect_obs_by_mod(mdata, mods = intersect_mods)
#    L.info("Remaining cells %d" % mdata.n_obs)
#    L.info(mdata.shape)



remove_unused_categories(mdata.obs)

# run quick test before saving out.
assert test_matching_df_ignore_cat(mdata.obs, orig_obs)  

# write out obs
output_prefix = re.sub(".h5mu", "", os.path.basename(args.output_mudata))
write_obs(mdata, output_prefix=os.path.join("tables/",output_prefix), output_suffix="_filtered_cell_metadata.tsv")

# write out the per sample_id cell numbers 
cell_counts_dict={}
for mm in mdata.mod.keys():
    cell_counts_dict[mm] = mdata[mod].obs.sample_id.value_counts().to_frame('n_cells')

cell_counts = pd.concat(cell_counts_dict).reset_index().rename(
    columns={"level_0": "modality", "level_1": "sample_id"})

L.info("cell_counts\n%s" %cell_counts)
cell_counts.to_csv("tables/" + output_prefix + "_cell_counts.csv", index=None)


mdata.update()
L.info("writing mdata h5mu")
L.info(args.output_mudata)
mdata.write(args.output_mudata)

L.info("Done")

