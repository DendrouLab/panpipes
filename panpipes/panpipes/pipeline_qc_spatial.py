#! /usr/bin/env python

from ruffus import *
import sys
import os
import re
from cgatcore import pipeline as P
import pandas as pd
import cgatcore.iotools as IOTools
from panpipes.funcs.io import gen_load_spatial_jobs
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

print(PARAMS)
#------------------------------------------------------------------------------------------------
## Create a dictionary of modalities
#------------------------------------------------------------------------------------------------
mode_dictionary = PARAMS["modalities"]
#{'spatialT': True}
print(mode_dictionary)
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

def gen_load_spatial_anndata_jobs():
    print(PARAMS['submission_file'])
    caf = pd.read_csv(PARAMS['submission_file'], sep='\t')
    print(gen_load_spatial_jobs(caf,mode_dictionary=PARAMS["modalities"]))
    return gen_load_spatial_jobs(caf,mode_dictionary=PARAMS["modalities"])
   


@follows(mkdir("logs"))
@follows(mkdir("tmp"))
@files(gen_load_spatial_anndata_jobs)
def load_mudatas(spatial_path, outfile, 
                 sample_id,
                 spatial_filetype, 
                 spatial_counts, 
                 spatial_metadata, 
                 spatial_transformation):
    
    path_dict = {'spatialT':spatial_path}
                 
    print(path_dict)
    print('sample_id = %s' % str(sample_id))
    print('outfile = %s' % str(outfile))
    print('spatial_filetype = %s' % str(spatial_filetype))
    print('spatial_counts = %s' % str(spatial_counts))
    print('spatial_metadata = %s' % str(spatial_metadata))
    print('spatial_transformation = %s' % str(spatial_transformation))

    modality_dict = {k:True if path_dict[k] is not None else False for k,v in PARAMS['modalities'].items() }
    print(modality_dict)
    
    cmd = """
        python %(py_path)s/make_mudataspatial_from_csv.py 
        --mode_dictionary "%(modality_dict)s"
        --sample_id %(sample_id)s
        --output_file %(outfile)s 
        --spatial_filetype %(spatial_filetype)s
        --spatial_infile %(spatial_path)s
        --spatial_counts %(spatial_counts)s
        --spatial_metadata %(spatial_metadata)s 
        --spatial_transformation %(spatial_transformation)s

        """
    cmd += " > logs/make_mudatas_%(sample_id)s.log"
    print(cmd)
    job_kwargs["job_threads"] = PARAMS['resources_threads_medium']
    P.run(cmd, **job_kwargs)

# in this workflow we qc each ST independently
# we need to use optional concatenation


@follows(load_mudatas)
@originate(unfilt_file())
def spatialQC(unfilt_file):
    resources_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), "resources")
    sample_prefix = PARAMS["sample_prefix"]
    cmd = """
            python %(py_path)s/run_scanpyQC_spatial.py
              --sampleprefix %(sample_prefix)s
              --input_anndata %(unfilt_file)s
              --outfile %(unfilt_file)s
              --figdir ./figures
              """
    if PARAMS['ccgenes'] is not None:
        if PARAMS['ccgenes'] == "default":
            ccgenesfile = resources_path + "/cell_cycle_genes.tsv"
        else:
            ccgenesfile = PARAMS['ccgenes']
        cmd += " --ccgenes %(ccgenesfile)s"
    if PARAMS['custom_genes_file'] is not None:
        if PARAMS['custom_genes_file'] == "default":
            customgenesfile = resources_path + "/qc_genelist_1.0.csv"
        else:
            customgenesfile = PARAMS['custom_genes_file']
        cmd += " --customgenesfile %(customgenesfile)s"
    if PARAMS['score_genes'] is not None:
        cmd += " --score_genes %(score_genes)s"
    if PARAMS['calc_proportions'] is not None:
        cmd += " --calc_proportions %(calc_proportions)s"


    cmd += " > logs/spatialQC_%(sample_id)s.log"
    job_kwargs["job_threads"] = PARAMS['resources_threads_medium']
    P.run(cmd, **job_kwargs)



    
#----- end

@follows(spatialQC)
def full():
    """
    All cgat pipelines should end with a full() function which updates,
    if needed, all branches of the pipeline.
    The @follows statement should ensure that all functions are covered,
    either directly or as prerequisites.
    """
    pass


def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))

    