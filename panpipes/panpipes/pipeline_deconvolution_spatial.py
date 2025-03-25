from ruffus import *
import sys
import os
from cgatcore import pipeline as P
import glob
import logging

def get_logger():
    return logging.getLogger("cgatcore.pipeline")


PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "pipeline.yml"])

PARAMS['py_path'] = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'python_scripts')
PARAMS['r_path'] = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'R_scripts')
PARAMS['resources_path'] = os.path.join(os.path.dirname(os.path.dirname(__file__)), "resources")

job_kwargs = {}

if PARAMS['condaenv'] is not None:
    job_kwargs["job_condaenv"] = PARAMS['condaenv']



if (PARAMS["input_spatial"] is not None) and (not os.path.exists(PARAMS["input_spatial"])):
    sys.exit("Spatial input directory doesn't exist or cannot be found")



def gen_filter_jobs():
    input_paths_spatial=glob.glob(os.path.join(PARAMS["input_spatial"],"*.zarr"))
    input_singlecell = PARAMS["input_singlecell"]
    for input_spatial in input_paths_spatial:
        sample_prefix = os.path.basename(input_spatial)
        sample_prefix = sample_prefix.replace(".zarr","")
        yield input_spatial, sample_prefix, input_singlecell    


@mkdir("logs")
@mkdir("figures")
@active_if(PARAMS['Cell2Location_run'] is True)
@mkdir("figures/Cell2Location")
@mkdir("cell2location.output")
@files(gen_filter_jobs)
def run_cell2location(input_spatial, sample_prefix, input_singlecell):

    figdir = "./figures/Cell2Location/" + sample_prefix
    output_dir = "./cell2location.output/" + sample_prefix
    log_file = "1_Cell2Location_" + sample_prefix + ".log"
    cmd = """
        python %(py_path)s/run_cell2location.py
        --input_spatial %(input_spatial)s
        --input_singlecell %(input_singlecell)s
        --figdir %(figdir)s
        --output_dir %(output_dir)s
              """
    # feature selection paramaters
    if PARAMS['Cell2Location_feature_selection']['gene_list'] is not None:
        cmd += f" --gene_list {PARAMS['Cell2Location_feature_selection']['gene_list']}"
    if PARAMS['Cell2Location_feature_selection']['remove_mt'] is not None:
        cmd += f" --remove_mt {PARAMS['Cell2Location_feature_selection']['remove_mt']}"
    if PARAMS['Cell2Location_feature_selection']['cell_count_cutoff'] is not None:
        cmd += f" --cell_count_cutoff {PARAMS['Cell2Location_feature_selection']['cell_count_cutoff']}"
    if PARAMS['Cell2Location_feature_selection']['cell_percentage_cutoff2'] is not None:
        cmd += f" --cell_percentage_cutoff2 {PARAMS['Cell2Location_feature_selection']['cell_percentage_cutoff2']}"
    if PARAMS['Cell2Location_feature_selection']['nonz_mean_cutoff'] is not None:
        cmd += f" --nonz_mean_cutoff {PARAMS['Cell2Location_feature_selection']['nonz_mean_cutoff']}"
    # parameters for the reference model
    if PARAMS['Cell2Location_reference']['labels_key'] is not None:
        cmd += f" --labels_key_reference {PARAMS['Cell2Location_reference']['labels_key']}"
    if PARAMS['Cell2Location_reference']['batch_key'] is not None:
        cmd += f" --batch_key_reference {PARAMS['Cell2Location_reference']['batch_key']}"
    if PARAMS['Cell2Location_reference']['layer'] is not None:
        cmd += f" --layer_reference {PARAMS['Cell2Location_reference']['layer']}"
    if PARAMS['Cell2Location_reference']['categorical_covariate_keys'] is not None:
        cmd += f" --categorical_covariate_keys_reference {PARAMS['Cell2Location_reference']['categorical_covariate_keys']}"
    if PARAMS['Cell2Location_reference']['continuous_covariate_keys'] is not None:
        cmd += f" --continuous_covariate_keys_reference {PARAMS['Cell2Location_reference']['continuous_covariate_keys']}"
    if PARAMS['Cell2Location_reference']['max_epochs'] is not None:
        cmd += f" --max_epochs_reference {PARAMS['Cell2Location_reference']['max_epochs']}"
    if PARAMS['Cell2Location_reference']['accelerator'] is not None:
        cmd += f" --accelerator_reference {PARAMS['Cell2Location_reference']['accelerator']}" 
    #if PARAMS['Cell2Location_reference']['use_gpu'] is not None:
    #    cmd += f" --use_gpu_reference {PARAMS['Cell2Location_reference']['use_gpu']}" 
    # parameters for the spatial model
    if PARAMS['Cell2Location_spatial']['batch_key'] is not None:
        cmd += f" --batch_key_st {PARAMS['Cell2Location_spatial']['batch_key']}"
    if PARAMS['Cell2Location_spatial']['layer'] is not None:
        cmd += f" --layer_st {PARAMS['Cell2Location_spatial']['layer']}"
    if PARAMS['Cell2Location_spatial']['categorical_covariate_keys'] is not None:
        cmd += f" --categorical_covariate_keys_st {PARAMS['Cell2Location_spatial']['categorical_covariate_keys']}"
    if PARAMS['Cell2Location_spatial']['continuous_covariate_keys'] is not None:
        cmd += f" --continuous_covariate_keys_st {PARAMS['Cell2Location_spatial']['continuous_covariate_keys']}"
    if PARAMS['Cell2Location_spatial']['max_epochs'] is not None:
        cmd += f" --max_epochs_st {PARAMS['Cell2Location_spatial']['max_epochs']}"
    if PARAMS['Cell2Location_spatial']['N_cells_per_location'] is not None:
        cmd += f" --N_cells_per_location {PARAMS['Cell2Location_spatial']['N_cells_per_location']}"
    if PARAMS['Cell2Location_spatial']['detection_alpha'] is not None:
        cmd += f" --detection_alpha {PARAMS['Cell2Location_spatial']['detection_alpha']}"
    if PARAMS['Cell2Location_spatial']['accelerator'] is not None:
        cmd += f" --accelerator_spatial {PARAMS['Cell2Location_spatial']['accelerator']}"
    
    if PARAMS['Cell2Location_save_models'] is not None:
        cmd += " --save_models %(Cell2Location_save_models)s"   
    if PARAMS['Cell2Location_export_gene_by_spot'] is not None:
        cmd += " --export_gene_by_spot %(Cell2Location_export_gene_by_spot)s"

    cmd += " > logs/%(log_file)s "
    job_kwargs["job_threads"] = PARAMS['resources_threads_low']
    log_msg = f"TASK: 'run_cell2location'" + f" IN CASE OF ERROR, PLEASE REFER TO : 'logs/{log_file}' FOR MORE INFORMATION."
    get_logger().info(log_msg)
    P.run(cmd, **job_kwargs)



@active_if(PARAMS['Tangram_run'] is True)
@mkdir("figures/Tangram")
@mkdir("tangram.output")
@files(gen_filter_jobs)
def run_tangram(input_spatial, sample_prefix, input_singlecell):

    figdir = "./figures/Tangram/" + sample_prefix
    output_dir = "./tangram.output/" + sample_prefix
    log_file = "2_Tangram_" + sample_prefix + ".log"
    cmd = """
        python %(py_path)s/run_tangram.py
        --input_spatial %(input_spatial)s
        --input_singlecell %(input_singlecell)s
        --figdir %(figdir)s
        --output_dir %(output_dir)s
              """
    # feature selection paramaters
    if PARAMS['Tangram_feature_selection']['gene_list'] is not None:
        cmd += f' --gene_list {PARAMS["Tangram_feature_selection"]["gene_list"]}'
    if PARAMS['Tangram_feature_selection']['rank_genes']['labels_key'] is not None:
        cmd += f' --labels_key_rank_genes {PARAMS["Tangram_feature_selection"]["rank_genes"]["labels_key"]}'
    if PARAMS['Tangram_feature_selection']['rank_genes']['layer'] is not None:
        cmd += f" --layer_rank_genes {PARAMS['Tangram_feature_selection']['rank_genes']['layer']}"
    if PARAMS['Tangram_feature_selection']['rank_genes']['n_genes'] is not None:
        cmd += f" --n_genes_rank {PARAMS['Tangram_feature_selection']['rank_genes']['n_genes']}"
    if PARAMS['Tangram_feature_selection']['rank_genes']['test_method'] is not None:
        cmd += f" --method_rank_genes {PARAMS['Tangram_feature_selection']['rank_genes']['test_method']}"
    if PARAMS['Tangram_feature_selection']['rank_genes']['correction_method'] is not None:
        cmd += f" --corr_method_rank_genes {PARAMS['Tangram_feature_selection']['rank_genes']['correction_method']}"    
    
    # model parameters 
    if PARAMS['Tangram_model']['labels_key'] is not None:
        cmd += f" --labels_key_model {PARAMS['Tangram_model']['labels_key']}"
    if PARAMS['Tangram_model']['num_epochs'] is not None:
        cmd += f" --num_epochs {PARAMS['Tangram_model']['num_epochs']}"
    if PARAMS['Tangram_model']['device'] is not None:
        cmd += f" --device {PARAMS['Tangram_model']['device']}"
    if PARAMS['Tangram_model']['kwargs'] is not None:
        kwargs = PARAMS['Tangram_model']['kwargs'].__str__().replace("'", '"')
        cmd += f" --kwargs '{kwargs}'"

    cmd += " > logs/%(log_file)s "
    job_kwargs["job_threads"] = PARAMS['resources_threads_low']
    log_msg = f"TASK: 'run_tangram'" + f" IN CASE OF ERROR, PLEASE REFER TO : 'logs/{log_file}' FOR MORE INFORMATION."
    get_logger().info(log_msg)
    P.run(cmd, **job_kwargs)


@follows(run_cell2location, run_tangram)
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
