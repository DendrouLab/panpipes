from ruffus import *
import sys
import os
from cgatcore import pipeline as P
import pandas as pd
import cgatcore.iotools as IOTools
import re
from itertools import chain
import glob
import warnings
import logging

PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "pipeline.yml"])

PARAMS['py_path'] = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'python_scripts')
PARAMS['r_path'] = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'R_scripts')
PARAMS['resources_path'] = os.path.join(os.path.dirname(os.path.dirname(__file__)), "resources")

job_kwargs = {}

if PARAMS['condaenv'] is not None:
    job_kwargs["job_condaenv"] = PARAMS['condaenv']


@mkdir("logs")
@mkdir("figures")
@active_if(PARAMS['Cell2Location_run'] is True)
#@originate(PARAMS['input_spatial'], PARAMS['input_singlecell'])
def run_cell2location():
    input_spatial = PARAMS["input_spatial"]
    input_singlecell = PARAMS["input_singlecell"]
    cmd = """
        python %(py_path)s/run_cell2location.py
        --input_spatial %(input_spatial)s
        --input_singlecell %(input_singlecell)s
        --figdir ./figures/Cell2Location
              """
    
    # feature selection paramaters
    if PARAMS['Cell2Location_gene_list'] is not None:
        cmd += " --gene_list %(Cell2Location_gene_list)s"
    if PARAMS['Cell2Location_remove_mt'] is not None:
        cmd += " --remove_mt %(Cell2Location_remove_mt)s"
    if PARAMS['Cell2Location_cell_count_cutoff'] is not None:
        cmd += " --cell_count_cutoff %(Cell2Location_cell_count_cutoff)s"
    if PARAMS['Cell2Location_cell_percentage_cutoff2'] is not None:
        cmd += " --cell_percentage_cutoff2 %(Cell2Location_cell_percentage_cutoff2)s"
    if PARAMS['Cell2Location_nonz_mean_cutoff'] is not None:
        cmd += " --nonz_mean_cutoff %(Cell2Location_nonz_mean_cutoff)s"
    
    # parameters for the reference model
    if PARAMS['Cell2Location_labels_key_reference'] is not None:
        cmd += " --labels_key_reference %(Cell2Location_labels_key_reference)s"
    if PARAMS['Cell2Location_batch_key_reference'] is not None:
        cmd += " --batch_key_reference %(Cell2Location_batch_key_reference)s"
    if PARAMS['Cell2Location_layer_reference'] is not None:
        cmd += " --layer_reference %(Cell2Location_layer_reference)s"
    if PARAMS['Cell2Location_categorical_covariate_keys_reference'] is not None:
        cmd += " --categorical_covariate_keys_reference %(Cell2Location_categorical_covariate_keys_reference)s"
    if PARAMS['Cell2Location_continuous_covariate_keys_reference'] is not None:
        cmd += " --continuous_covariate_keys_reference %(Cell2Location_continuous_covariate_keys_reference)s"
    if PARAMS['Cell2Location_max_epochs_reference'] is not None:
        cmd += " --max_epochs_reference %(Cell2Location_max_epochs_reference)s"
    
        # parameters for the spatial model
    if PARAMS['Cell2Location_batch_key_st'] is not None:
        cmd += " --batch_key_st %(Cell2Location_batch_key_st)s"
    if PARAMS['Cell2Location_layer_st'] is not None:
        cmd += " --layer_st %(Cell2Location_layer_st)s"
    if PARAMS['Cell2Location_categorical_covariate_keys_st'] is not None:
        cmd += " --categorical_covariate_keys_st %(Cell2Location_categorical_covariate_keys_st)s"
    if PARAMS['Cell2Location_continuous_covariate_keys_st'] is not None:
        cmd += " --continuous_covariate_keys_st %(Cell2Location_continuous_covariate_keys_st)s"
    if PARAMS['Cell2Location_max_epochs_st'] is not None:
        cmd += " --max_epochs_st %(Cell2Location_max_epochs_st)s"
    if PARAMS['Cell2Location_N_cells_per_location'] is not None:
        cmd += " --N_cells_per_location %(Cell2Location_N_cells_per_location)s"
    if PARAMS['Cell2Location_detection_alpha'] is not None:
        cmd += " --detection_alpha %(Cell2Location_detection_alpha)s"
    
    if PARAMS['Cell2Location_save_models'] is not None:
        cmd += " --save_models %(Cell2Location_save_models)s"
    
    cmd += " > logs/deconvolution.log "
    job_kwargs["job_threads"] = PARAMS['resources_threads_low']
    P.run(cmd, **job_kwargs)


@follows(run_cell2location)
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
