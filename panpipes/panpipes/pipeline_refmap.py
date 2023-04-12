from venv import create
from ruffus import *
import sys
import os
from cgatcore import pipeline as P
import pandas as pd
import cgatcore.iotools as IOTools
import re
from itertools import chain
import glob

# __file__="/well/cartography/users/zsj686/non_cart_projects/005-multimodal_scpipelines/src/sc_pipelines_muon_dev/panpipes/pipeline_refmap.py"
PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "pipeline.yml"])

job_kwargs = {}

if PARAMS['condaenv'] is not None:
    job_kwargs["job_condaenv"] =PARAMS['condaenv']


PARAMS['py_path'] =  os.path.join(os.path.dirname(os.path.dirname(__file__)), 'python_scripts')
PARAMS['r_path'] = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'R_scripts')
PARAMS['resources_path'] = os.path.join(os.path.dirname(os.path.dirname(__file__)), "resources")



if PARAMS['totalvi'] is not None:
    totalvi=True
    totalvi_models = PARAMS['totalvi']

if PARAMS['scvi'] is not None:
    if PARAMS['scanvi'] is not None:
        scvi=True


@originate("logs/setup_dirs.sentinel")
def set_up_dirs(log_file):
    os.mkdir("models")
    os.mkdir("refmap")
    os.mkdir("logs")
    os.mkdir("figures")
    IOTools.touch_file(log_file)
    pass

# ------------------------------
# create cluster job for each model in scvi tools: scvi, scanvi, totalvi
# ------------------------------

def gen_refmap_jobs():
    for ref_architecture in ["scvi", "scanvi", "totalvi"]:
        if PARAMS[ref_architecture] is not None:
            mod_list = PARAMS[ref_architecture]
            if type(mod_list) is not list:
                mod_list = [mod_list]
            for infile in mod_list:
                # this is creating output files based on the model folder name
                out_prefix=os.path.basename(os.path.dirname(infile))
                outfile = os.path.join('refmap', out_prefix + ".csv")
                log_file = os.path.join('logs/', out_prefix + ".log")
                yield infile, outfile, log_file, ref_architecture
                


# error if modality is not rna

# ------------------------------
# refmap scvi tools: scvi, scanvi, totalvi
# ------------------------------

@follows(set_up_dirs)
@files(gen_refmap_jobs)
def run_refmap_scvi(infile, outfile, log_file, ref_architecture ):
    cmd = """
    python %(py_path)s/refmap_scvitools.py
    --query_data %(query)s
    --reference_path %(infile)s
    --reference_architecture  %(ref_architecture)s

    --neighbors_n_pcs %(neighbors_npcs)s
    --neighbors_method %(neighbors_method)s
    --neighbors_k %(neighbors_k)s
    --neighbors_metric %(neighbors_metric)s
    """
    if PARAMS['transform_batch']:
        cmd += " --transform_batch %(query_batch)s"
    if PARAMS['reference_data'] is not None:
        cmd += " --adata_reference %(reference_data)s "
    if ref_architecture == "totalvi":
       cmd += "  --impute_proteins %(impute_proteins)s"
    cmd += " > %(log_file)s"
    job_kwargs["job_threads"] = PARAMS['resources_threads_low']
    P.run(cmd, **job_kwargs)



# ------------------------------
# refmap Seurat wnn
# ------------------------------

@follows(set_up_dirs)
@originate("logs/refmap_wnn.log")
def run_refmap_wnn(outfile):
    cmd = """
        Rscript %(r_path)s/refmap_wnn.R
            --reference_path %(wnn)s
            --query_path %(query_path)s 
          """  
    if PARAMS["wnn_mapping_normalization"] is not None:
        cmd +=" --wnn_mapping_normalization %s" % PARAMS["wnn_mapping_normalization"]
    if PARAMS["normalize"]:
        cmd +=" --normalize TRUE"
    else:
        cmd +=" --normalize FALSE"


    job_kwargs["job_threads"] = PARAMS['resources_threads_low']
    
    P.run(cmd, **job_kwargs)
   
#-------------------------------
# collate umaps here (actually not done)
#-------------------------------

# this COLLATE job will become the big final collate job across all possible combinations you may have run

@collate([run_refmap_scvi, run_refmap_wnn], #multimodal
         regex(r"refmap/(.*)"), # this does nothing
          r'refmap/combined_umaps.tsv')
def plot_umaps(infiles,outfile):
    infiles_string = ','.join(infiles)
    cell_mtd_file = "cell_mtd.csv"
    cmd = """python %(py_path)s/refmap_plot_umaps.py 
    --input_mudata %(query)s
    --input_umap_files %(infiles_string)s 
    --output_combined_umaps_tsv %(outfile)s
    """ 
    cmd += " > logs/collate_outputs.log "
    job_kwargs["job_threads"] = PARAMS['resources_threads_low']
    P.run(cmd, **job_kwargs)  



@follows(plot_umaps)
def full():
    """
    All cgat pipelines should end with a full() function which updates,
    if needed, all branches of the pipeline.
    The @follows statement should ensure that all functions are covered,
    either directly or as prerequisites.
    """
    #IOTools.touch_file(file)
    pass

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
