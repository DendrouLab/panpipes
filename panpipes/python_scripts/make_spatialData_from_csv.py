import argparse
import yaml
# import scanpy as sc
import pandas as pd
# import numpy as np
# from scipy.sparse import csr_matrix
import muon as mu
import warnings
from muon._atac.tools import add_peak_annotation, locate_fragments
import squidpy as sq
import spatialdata_io as sd_io
from mudata import MuData
import os
"""
this script copies the make_adata_from_csv.py that creates
ONE MUDATA PER SAMPLE, with in each ONE LAYER per modality
for cell-suspension, saves them to temp. 
concatenation of the mudatas saved in tmp happens 
in the concat_anndata.py script
"""

import sys
import logging

L = logging.getLogger()
L.setLevel(logging.INFO)
log_handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s: %(levelname)s - %(message)s')
log_handler.setFormatter(formatter)
L.addHandler(log_handler)

from panpipes.funcs.io import check_filetype, read_yaml


parser = argparse.ArgumentParser()
#parser.add_argument('--mode_dictionary',default=mode_dictionary,type=json.loads)
parser.add_argument('--mode_dictionary',default="",
                    help="accepts a yaml speficication of modalities dictionary", 
                    type=yaml.safe_load)
parser.add_argument('--sample_id',
                    default=None,
                    help='')
parser.add_argument('--output_file',
                    default=None,
                    help='')
# start of spatial args
parser.add_argument('--spatial_infile', 
                    default=None,
                    help='')
parser.add_argument('--spatial_filetype', 
                    default=None,
                    help='')
parser.add_argument('--spatial_counts', 
                    default=None,
                    help='')
parser.add_argument('--scalefactors_file', 
                    default=None,
                    help='')
parser.add_argument('--fullres_image_file', 
                    default=None,
                    help='')
parser.add_argument('--tissue_positions_file', 
                    default=None,
                    help='')
parser.add_argument('--spatial_metadata', 
                    default=None,
                    help='')
parser.add_argument('--spatial_transformation', 
                    default=None,
                    help='')

parser.set_defaults(verbose=True)
args, opt = parser.parse_known_args()
L.info("Running with params: %s", args)

# unimodal mu (check if all the modalities)
#if isinstance(args.mode_dictionary, dict):
#    mode_dictionary = args.mode_dictionary
#else:
#    mode_dictionary = read_yaml(args.mode_dictionary) 
#{'spatialT': True}

#permf = [key for key, value in mode_dictionary.items() if value == True]
all_files = {
            "spatial":[args.spatial_infile, #path
                        args.spatial_filetype, #needed for the load_adata_in function to call one of vizgen,visium
                        args.spatial_counts, #name of the counts file, mandatory for squidpy
                        args.fullres_image_file, # visium
                        args.tissue_positions_file, #visium
                        args.scalefactors_file]} # visium 
#                        args.spatial_metadata, #name of the metadata file, mandatory for squidpy
#                        args.spatial_transformation]}
#subset to the modalities we want from permf (in this case only spatial)
#all_files = {nm: x  for (nm, x) in all_files.items() if nm in permf}

#[check_filetype(x[0], x[1]) for x in all_files.values()]
# read the spatial data with one of the functions inside
# load_mdata_from_multiple_files
#     |
#      -------->load_adata_in
# this function creates ONE mudata per row of the CAF file 
# and saves it with sample_id.h5mu in tmp/ 


def create_soft_link(source_path, target_path):
    try:
        os.symlink(source_path, target_path)
        print(f"Soft link created: {target_path}")
    except FileExistsError:
        print(f"Soft link already exists: {target_path}")
    except Exception as e:
        print(f"An error occurred while creating soft link: {e}")

def check_dir_transform(infile_path, transform_file):
    images_path = os.path.join(infile_path, "images")
    if not os.path.exists(images_path):
        os.makedirs(images_path)
        print(f"Directory '{images_path}' created.")
    target_path = os.path.join(images_path, os.path.basename(transform_file))
    tocheck_orig= os.path.join(infile_path,transform_file)
    if os.path.exists(target_path):
        print("File already exists in 'images' directory.")
    elif os.path.islink(tocheck_orig):
        original_path = os.path.realpath(tocheck_orig)
        create_soft_link(original_path, target_path)
    else:
        create_soft_link(os.path.abspath(transform_file), target_path)



if args.spatial_filetype=="vizgen":
    L.info("Reading in Vizgen data with squidpy.read.vizgen() into AnnData from directory " + args.spatial_infile)
    sdata = sd_io.merscope(path = args.spatial_infile)
#    adata = sq.read.vizgen(path = args.spatial_infile, #path, mandatory for squidpy
#                        counts_file=args.spatial_counts, #name of the counts file, mandatory for squidpy
#                        meta_file = args.spatial_metadata, #name of the metadata file, mandatory for squidpy
#                        transformation_file=args.spatial_transformation,
#                        library_id = str(args.sample_id)) #this also has kwargs for read_10x_h5 but keep simple
#    adata.uns["spatial"][str(args.sample_id)]["scalefactors"]["transformation_matrix"].columns = adata.uns["spatial"][str(args.sample_id)]["scalefactors"]["transformation_matrix"].columns.astype(str)
elif args.spatial_filetype =="visium":
    L.info("Reading in Visium data with squidpy.read.visium() into AnnData from directory " + args.spatial_infile)
    sdata = sd_io.visium(path=args.spatial_infile, 
                         dataset_id=str(args.sample_id), 
                         counts_file=args.spatial_counts, 
                         fullres_image_file=args.fullres_image_file,
                         tissue_positions_file=args.tissue_positions_file, 
                         scalefactors_file=args.scalefactors_file)
    #adata = sq.read.visium(path = args.spatial_infile, #path, mandatory for squidpy
    #                    counts_file=args.spatial_counts, #name of the counts file, mandatory for squidpy
    #                    library_id = str(args.sample_id)
    #                    ) #this also has kwargs for read_10x_h5 but keep simple

L.info("Resulting SpatialData is:")
L.info(sdata)
#L.info("Creating MuData with .mod['spatial']")

#mdata = MuData({"spatial": adata})


#---------------
# do some extra processing on the different modalities as needed
#---------------

L.info("Making var names unique")
#make var names unique
#for mm in mdata.mod.keys():
sdata["table"].var_names_make_unique()

L.info("Adding sample_id '%s'to MuData.obs and MuData.mod['spatial'].obs" % args.sample_id)
sdata["table"].obs['sample_id'] = str(args.sample_id)

# copy the sample_id to each modality
#for mm in mdata.mod.keys():
    # mdata[mm].obs['sample_id'] = mdata.obs['sample_id']
sdata["table"].obs['sample_id'] = sdata["table"].obs.loc[sdata["table"].obs_names,:]['sample_id']

#mdata.update()

L.info("Resulting SpatialData is:")
L.info(sdata)

L.info("Saving SpatialData to '%s'" % args.output_file)
L.debug(sdata)
sdata.write(args.output_file)

L.info("Done")


