import scanpy as sc
import argparse
import pandas as pd
import re
import muon as mu
from anndata import AnnData
import spatialdata as sd
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

parser.add_argument('--input_spatialdata',
                    default='gut_minus1_amp.h5ad',
                    help='')
parser.add_argument('--output_spatialdata',
                    default='',
                    help='')
parser.add_argument('--filter_dict',
                    default='',
                    help='this is pull')
# cross modalities args
parser.add_argument('--keep_barcodes', default=None,
                    help='1 column list of barcodes to keep, note that they should match the spatialdata input, this filtering happens first')


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

# load spatialdata

L.info("Reading in SpatialData from '%s'" % args.input_spatialdata)
sdata = sd.read_zarr(args.input_spatialdata)
#mdata = mu.read(args.input_spatialdata)

#if isinstance(mdata, AnnData):
#    raise TypeError("Input '%s' should be of spatialdata format, not Anndata"  % args.input_spatialdata)

orig_obs = sdata["table"].obs.copy()

L.info("Before filtering: "+ str(sdata["table"].n_obs) + " cells and " + str(sdata["table"].n_vars) + " features")

# filter based on provided barcodes -----
if args.keep_barcodes is not None:
    L.info("Filtering SpatialData by keep_barcodes file")
    keep_bc = pd.read_csv(args.keep_barcodes,header=None)
    sdata["table"] = sdata["table"][sdata["table"].obs_names.isin(keep_bc[0]),:].copy()
    remove_unused_categories(sdata["table"].obs)
    #mdata.update()
    L.info("Remaining cells: %d" % sdata["table"].n_obs)



# filter more than
if filter_dict['run']:
    for marg in filter_dict["spatial"].keys():
        if marg == "obs":
            if "max" in filter_dict["spatial"][marg].keys():
                for col, n in filter_dict["spatial"][marg]['max'].items():
                    L.info("Filtering cells of modality '%s' by '%s' in .obs to less than %s" % ("spatial", col, n))
                    mu.pp.filter_obs(sdata["table"], col, lambda x: x <= n)
                    L.info("Remaining cells: %d" % sdata["table"].n_obs)
            if "min" in filter_dict["spatial"][marg].keys():
                for col, n in filter_dict["spatial"][marg]['min'].items():
                    L.info("Filtering cells of modality '%s' by '%s' in .obs to more than %s" % ("spatial", col, n))
                    mu.pp.filter_obs(sdata["table"], col, lambda x: x >= n)
                    L.info("Remaining cells: %d" % sdata["table"].n_obs)
            if "bool" in filter_dict["spatial"][marg].keys():
                for col, n in filter_dict["spatial"][marg]['bool'].items():
                    L.info("Filtering cells of modality '%s' by '%s' in .obs marked %s" % ("spatial", col, n))
                    mu.pp.filter_obs(sdata["table"], col, lambda x: x == n)
                    L.info("Remaining cells: %d" % sdata["table"].n_obs)
        if marg == "var":
            if "max" in filter_dict["spatial"][marg].keys():
                for col, n in filter_dict["spatial"][marg]['max'].items():
                    L.info("Filtering features of modality '%s' by '%s' in .var to less than %s" % ("spatial", col, n))
                    mu.pp.filter_var(sdata["table"], col, lambda x: x <= n)
                    L.info("Remaining features: %d" % sdata["table"].n_vars)

            if "min" in filter_dict["spatial"][marg].keys():
                for col, n in filter_dict["spatial"][marg]['min'].items():
                    L.info("Filtering features of modality '%s' by '%s' in .var to more than %s" % ("spatial", col, n))
                    mu.pp.filter_var(sdata["table"], col, lambda x: x >= n)
                    L.info("Remaining features: %d" % sdata["table"].n_vars)

            if "bool" in filter_dict["spatial"][marg].keys():
                for col, n in filter_dict["spatial"][marg]['bool'].items():
                    L.info("Filtering features of modality '%s' by '%s' in .var marked %s" % ("spatial", col, n))
                    mu.pp.filter_var(sdata["table"], col, lambda x: x == n)
                    L.info("Remaining features: %d" % sdata["table"].n_vars)



#mdata.update()

L.info("After filtering: "+ str(sdata["table"].n_obs) + " cells and " + str(sdata["table"].n_vars) + " features")

remove_unused_categories(sdata["table"].obs)

# run quick test before saving out.
assert test_matching_df_ignore_cat(sdata["table"].obs, orig_obs)  

# write out obs
output_prefix = re.sub(".zarr", "", os.path.basename(args.output_spatialdata))

L.info("Saving updated obs in a metadata tsv file to './tables/" + output_prefix + "_filtered_cell_metadata.tsv'")
write_obs(sdata["table"], output_prefix=os.path.join("tables/",output_prefix), output_suffix="_filtered_cell_metadata.tsv")

# write out the per sample_id cell numbers 
cell_counts_dict={}
#for mm in mdata.mod.keys():
cell_counts_dict["spatial"] = sdata["table"].obs.sample_id.value_counts().to_frame('n_cells')

cell_counts = pd.concat(cell_counts_dict).reset_index().rename(
    columns={"level_0": "modality", "level_1": "sample_id"})

L.info("cell_counts: \n%s" %cell_counts)
L.info("Saving cell counts in a metadata csv file to './tables/" + output_prefix + "_cell_counts.csv'")
cell_counts.to_csv("tables/" + output_prefix + "_cell_counts.csv", index=None)

#mdata.update()

L.info("Saving updated SpatialData to '%s'" % args.output_spatialdata)
sdata.write(args.output_spatialdata)

L.info("Done")

