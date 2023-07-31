#! /usr/bin/env python

from ruffus import *
import sys
import os
import re
from cgatcore import pipeline as P
import pandas as pd
import cgatcore.iotools as IOTools
from panpipes.funcs.io import check_submission_file, gen_load_anndata_jobs
# from scpipelines.funcs.processing import intersection
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning) 

import warnings

# warnings.simplefilter('ignore', category=NumbaDeprecationWarning)
# warnings.simplefilter('ignore', category=NumbaPendingDeprecationWarning)

import yaml


PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "pipeline.yml"])

PARAMS['py_path'] =  os.path.join(os.path.dirname(os.path.dirname(__file__)), 'python_scripts')
PARAMS['r_path'] = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'R_scripts')
PARAMS['resources_path'] = os.path.join(os.path.dirname(os.path.dirname(__file__)), "resources")
job_kwargs = {}

if PARAMS['condaenv'] is not None:
    job_kwargs["job_condaenv"] =PARAMS['condaenv']


#------------------------------------------------------------------------------------------------
## Create a dictionary of modalities
#------------------------------------------------------------------------------------------------
mode_dictionary = PARAMS["modalities"]
#{'spatialT': True}

#------------------------------------------------------------------------------------------------
# setup dirs
#------------------------------------------------------------------------------------------------

@originate("logs/setup_dirs.sentinel")
def set_up_dirs(log_file):
    os.makedirs("logs", exist_ok=True)
    os.makedirs("tmp", exist_ok=True)
    os.makedirs("figures", exist_ok=True)
    IOTools.touch_file(log_file)
    pass

# -----------------------------------------------------------------------------------------------
## Creating h5mu from filtered data files
# -----------------------------------------------------------------------------------------------

def unfilt_file():
    sprefix = PARAMS['sample_prefix']
    unfilt_file = sprefix + "_unfilt.h5mu"
    return unfilt_file



def gen_load_filtered_anndata_jobs():
    caf = pd.read_csv(PARAMS['submission_file'], sep='\t')
    #check with cellranger
    return gen_load_anndata_jobs(caf, load_raw=False, 
                                 mode_dictionary=PARAMS["modalities"])

    

@follows(mkdir("logs"))
@follows(mkdir("tmp"))
@active_if(PARAMS["use_existing_h5mu"] is False)
@files(gen_load_filtered_anndata_jobs)
def load_mudatas(sample_id, outfile, 
                 spatial_path, 
                 spatial_filetype, 
                 spatial_counts, 
                 spatial_metadata, 
                 spatial_transformation):
    
    path_dict = {'spatialT':spatial_path}
                 
    print(path_dict)
    
    modality_dict = {k:True if path_dict[k] is not None else False for k,v in PARAMS['modalities'].items() }
    print(modality_dict)
    
    cmd = """
        python %(py_path)s/make_mudataspatial_from_csv.py 
        --mode_dictionary "%(modality_dict)s"
        --sample_id %(sample_id)s
        --output_file %(outfile)s 
        --spatial_filetype %(spatial_filetype)s
        --spatial_infile %(spatial_infile)s
        --spatial_counts %(spatial_counts)s
        --spatial_metadata %(spatial_metadata)s 
        --spatial_transformation %(spatial_transformation)s

    """
    

    