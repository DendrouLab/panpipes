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
import logging

# warnings.simplefilter('ignore', category=NumbaDeprecationWarning)
# warnings.simplefilter('ignore', category=NumbaPendingDeprecationWarning)

import yaml

def get_logger():
    return logging.getLogger("cgatcore.pipeline")

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
mode_dictionary = {'spatial': True} #PARAMS["modalities"]
#{'spatial': True}
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
## Creating spatialData from filtered data files
# -----------------------------------------------------------------------------------------------


def gen_load_spatial_anndata_jobs():
    print(PARAMS['submission_file'])
    caf = pd.read_csv(PARAMS['submission_file'], sep='\t')
    print(gen_load_spatial_jobs(caf,mode_dictionary={'spatial': True}))
    global assays
    assays = {}
    return gen_load_spatial_jobs(caf,mode_dictionary={'spatial': True})
   


@follows(mkdir("logs"))
@follows(mkdir("tmp"))
@files(gen_load_spatial_anndata_jobs)
def load_spatialdatas(spatial_path,  outfile, 
                 sample_id, spatial_filetype, visium_feature_bc_matrix, visium_fullres_image_file, visium_tissue_positions_file, visium_scalefactors_file,
              vpt_cell_by_gene, vpt_cell_metadata, vpt_cell_boundaries):
    
    path_dict = {'spatial':spatial_path}
                 
    print(path_dict)
    print('sample_id = %s' % str(sample_id))
    print('outfile = %s' % str(outfile))
    print('spatial_filetype = %s' % str(spatial_filetype))

    if spatial_filetype == "visium":
        print('visium_feature_bc_matrix = %s' % str(visium_feature_bc_matrix))
        print('visium_fullres_image_file= %s' % str(visium_fullres_image_file))
        print('visium_tissue_positions_file= %s' % str(visium_tissue_positions_file))
        print('visium_scalefactors_file= %s' % str(visium_scalefactors_file))
    if spatial_filetype == "vizgen":
        print('vpt_cell_by_gene = %s' % str(vpt_cell_by_gene))
        print('vpt_cell_metadata= %s' % str(vpt_cell_metadata))
        print('vpt_cell_boundaries= %s' % str(vpt_cell_boundaries))
    modality_dict = {k:True if path_dict[k] is not None else False for k,v in {'spatial': True}.items() }
    print(modality_dict)

    assays[outfile] = spatial_filetype
    
    cmd = """
        python %(py_path)s/make_spatialData_from_csv.py 
        --mode_dictionary "%(modality_dict)s"
        --sample_id %(sample_id)s
        --output_file %(outfile)s 
        --spatial_filetype %(spatial_filetype)s
        --spatial_infile %(spatial_path)s
    """
    if spatial_filetype == "visium":
        cmd += """
        --visium_feature_bc_matrix %(visium_feature_bc_matrix)s
        --scalefactors_file %(visium_scalefactors_file)s 
        --fullres_image_file %(visium_fullres_image_file)s
        --tissue_positions_file %(visium_tissue_positions_file)s
        """
    if spatial_filetype == "vizgen":
        cmd += """
        --vpt_cell_by_gene %(vpt_cell_by_gene)s
        --vpt_cell_metadata %(vpt_cell_metadata)s 
        --vpt_cell_boundaries %(vpt_cell_boundaries)s
        """
    cmd += " > logs/1_make_spatialdatas_%(sample_id)s.log"
    job_kwargs["job_threads"] = PARAMS['resources_threads_medium']
    log_msg = f"TASK: 'load_spatialdatas'" + f" IN CASE OF ERROR, PLEASE REFER TO : 'logs/1_make_spatialdatas_{sample_id}.log' FOR MORE INFORMATION."
    get_logger().info(log_msg)
    P.run(cmd, **job_kwargs)




@follows(load_spatialdatas)
@follows(mkdir("qc.data"))
@follows(mkdir("./figures"))
@transform(load_spatialdatas,
           regex("./tmp/(.*)_raw.zarr"), 
           r"./logs/2_spatialQC_\1.log")
def spatialQC(infile,log_file):
    spatial_filetype = assays[infile]
    resources_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), "resources")
    outfile = infile.replace("_raw","_unfilt")
    outfile = outfile.replace("tmp", "qc.data")
    cmd = """
            python %(py_path)s/run_scanpyQC_spatial.py
              --input_anndata %(infile)s
              --spatial_filetype %(spatial_filetype)s
              --outfile %(outfile)s
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
    cmd += " > %(log_file)s"
    job_kwargs["job_threads"] = PARAMS['resources_threads_medium']
    log_msg = f"TASK: 'spatialQC'" + f" IN CASE OF ERROR, PLEASE REFER TO : '{log_file}' FOR MORE INFORMATION."
    get_logger().info(log_msg)
    P.run(cmd, **job_kwargs)


def run_plotqc_query(pqc_dict):
    # avoid deleting from PARAMS
    pqc = pqc_dict.copy()
    del pqc['grouping_var']
    return any([x != None for x in pqc.values()])


@follows(spatialQC)
@follows(mkdir("./figures/spatial"))
@active_if(run_plotqc_query(PARAMS['plotqc']))
@transform(load_spatialdatas, 
           regex("./tmp/(.*)_raw.zarr"),
           r"./logs/3_qcplot.\1.log")
def plotQC_spatial(unfilt_file,log_file):
    spatial_filetype = assays[unfilt_file]
    unfilt_file = unfilt_file.replace("_raw","_unfilt")
    unfilt_file = unfilt_file.replace("tmp", "qc.data")
    cmd = """
            python %(py_path)s/plot_qc_spatial.py
             --input_spatialdata %(unfilt_file)s
             --spatial_filetype %(spatial_filetype)s
             --figdir ./figures/spatial
            """
    if PARAMS['plotqc']['grouping_var'] is not None:
        cmd += " --grouping_var %(plotqc_grouping_var)s"
    if PARAMS['plotqc']['spatial_metrics'] is not None:
        cmd += " --spatial_qc_metrics %(plotqc_spatial_metrics)s"
    cmd += " > %(log_file)s "
    job_kwargs["job_threads"] = PARAMS['resources_threads_low']
    log_msg = f"TASK: 'plotQC_spatial'" + f" IN CASE OF ERROR, PLEASE REFER TO : '{log_file}' FOR MORE INFORMATION."
    get_logger().info(log_msg)
    P.run(cmd, **job_kwargs)

    
#----- end

@follows(plotQC_spatial)
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

    