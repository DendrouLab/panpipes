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


PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "pipeline.yml"])


PARAMS['py_path'] =  os.path.join(os.path.dirname(os.path.dirname(__file__)), 'python_scripts')
PARAMS['r_path'] = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'R_scripts')
PARAMS['resources_path'] = os.path.join(os.path.dirname(os.path.dirname(__file__)), "resources")

sprefix = PARAMS['sample_prefix']

job_kwargs={}
if PARAMS['condaenv'] is not None:
    job_kwargs["job_condaenv"] =PARAMS['condaenv']


preprocessed_file = PARAMS['preprocessed_file']

@originate("logs/setup_dirs.sentinel")
def set_up_dirs(log_file):
    os.makedirs("logs", exist_ok=True)
    os.makedirs("tmp", exist_ok=True)
    os.makedirs("batch_correction", exist_ok=True)
    if PARAMS['rna']['run']:
        os.makedirs("figures/rna")
    if PARAMS['prot']['run']:
        os.makedirs("figures/prot")
    if PARAMS['atac']['run']:
        os.makedirs("figures/atac")
    if PARAMS['multimodal']['run']:
        os.makedirs("figures/multimodal")
    os.makedirs("figures/rep")
    IOTools.touch_file(log_file)
    pass
# ------------------------------------------------------------------------
# Unimodal Integration
# ------------------------------------------------------------------------
#
#
# ------------------------------------------------------------------------
# Rna unimodal integration:
# - nobatch
# - BBKNN
# - COMBAT
# - Harmony
# - Scanorama
# - scvi
# ------------------------------------------------------------------------
#rna No batch correction
@follows(set_up_dirs)
@originate("batch_correction/umap_rna_none.csv")
def run_no_batch_umap(outfile):
    # print(PARAMS)
    rna_params = PARAMS['rna']

    cmd = """python %(py_path)s/batch_correct_none.py 
     --input_anndata %(preprocessed_obj)s
     --output_csv %(outfile)s
     --integration_col %(rna_column)s
     """
    # cannot use the normal method for importing params from yaml, because it only works up to depth 2
    neighbor_params = PARAMS['rna']['neighbors']
    if neighbor_params['method'] is not None:
        cmd += " --neighbors_method %s" % neighbor_params['method']
    if neighbor_params['metric'] is not None:
        cmd += " --neighbors_metric %s" % neighbor_params['metric']
    if neighbor_params['npcs'] is not None:
        cmd += " --neighbors_n_pcs %s"  % neighbor_params['npcs']
    if neighbor_params['k'] is not None:
        cmd += " --neighbors_k %s" % neighbor_params['k']
    cmd += " > logs/rna_no_correct.log"
    
    if PARAMS['queues_long'] is not None:
        job_kwargs["job_queue"] = job_queue=PARAMS['queues_long']
    job_kwargs["job_threads"] = PARAMS['resources_threads_high']    
    P.run(cmd, **job_kwargs) 



# rna BBKNN
@follows(set_up_dirs)
@active_if(PARAMS['rna_run'])
@active_if(PARAMS['rna_tools'] is not None and 'bbknn' in PARAMS['rna_tools'])
# @follows(normalise_log_hvg_regress_scale)
@originate("batch_correction/umap_rna_bbknn.csv")
def run_bbknn_rna(outfile):
    cmd = """python %(py_path)s/batch_correct_bbknn.py 
     --input_anndata %(preprocessed_obj)s
     --output_csv %(outfile)s
     --integration_col %(rna_column)s
     --modality rna
     """
    if PARAMS['rna']['bbknn']['neighbors_within_batch'] is not None:
        cmd += " --neighbors_within_batch %i" % PARAMS['rna']['bbknn']['neighbors_within_batch']
    if PARAMS['rna']['neighbors']['npcs'] is not None:
        cmd += " --neighbors_n_pcs %s" % PARAMS['rna']['neighbors']['npcs']
    cmd += " > logs/rna_bbknn.log "
    
    if PARAMS['queues_long'] is not None:
        job_kwargs["job_queue"] = job_queue=PARAMS['queues_long']
    job_kwargs["job_threads"] = PARAMS['resources_threads_high']
    P.run(cmd, **job_kwargs) 


# rna COMBAT
@follows(set_up_dirs)
@active_if(PARAMS['rna_run'])
@active_if(PARAMS['rna_tools'] is not None and 'combat' in PARAMS['rna_tools'])
@originate("batch_correction/umap_rna_combat.csv")
def run_combat(outfile):
    cmd = """python %(py_path)s/batch_correct_combat.py 
     --input_anndata %(preprocessed_obj)s
     --output_csv %(outfile)s
     --integration_col %(rna_column)s
     --n_threads %(resources_threads_high)s
     --modality rna
     """
    # cannot use the normal method for importing params from yaml, because it only works up to depth 2
    neighbor_params = PARAMS['rna']['neighbors']
    if neighbor_params['method'] is not None:
        cmd += " --neighbors_method %s" % neighbor_params['method']
    if neighbor_params['metric'] is not None:
        cmd += " --neighbors_metric %s" % neighbor_params['metric']
    if neighbor_params['npcs'] is not None:
        cmd += " --neighbors_n_pcs %s"  % neighbor_params['npcs']
    if neighbor_params['k'] is not None:
        cmd += " --neighbors_k %s" % neighbor_params['k']
    cmd += " > logs/rna_combat.log "
    
    if PARAMS['queues_long'] is not None:
        job_kwargs["job_queue"] = job_queue=PARAMS['queues_long']

    job_kwargs["job_threads"] = PARAMS['resources_threads_high']
    P.run(cmd, **job_kwargs) 

# rna HARMONY
@follows(set_up_dirs)
@active_if(PARAMS['rna_run'])
@active_if(PARAMS['rna_tools'] is not None and 'harmony' in PARAMS['rna_tools'])
@originate("batch_correction/umap_rna_harmony.csv")
def run_harmony(outfile):
    cmd = """python %(py_path)s/batch_correct_harmony.py 
     --input_anndata %(preprocessed_obj)s
     --output_csv %(outfile)s 
     --integration_col %(rna_column)s
     --n_threads %(resources_threads_high)s
     --modality rna
     """
    # cannot use the normal method for importing params from yaml, because it only works up to depth 2
    harmony_params = PARAMS['rna']['harmony']
    if harmony_params['npcs'] is not None:
        cmd += " --harmony_npcs %s" % harmony_params['npcs']
    if harmony_params['sigma'] is not None:
        cmd += " --sigma_val %s" % harmony_params['sigma'] 
    if harmony_params['theta'] is not None:
        cmd += " --theta_val %s" % harmony_params['theta']       
    neighbor_params = PARAMS['rna']['neighbors']
    if neighbor_params['method'] is not None:
        cmd += " --neighbors_method %s" % neighbor_params['method']
    if neighbor_params['metric'] is not None:
        cmd += " --neighbors_metric %s" % neighbor_params['metric']
    if neighbor_params['npcs'] is not None:
        cmd += " --neighbors_n_pcs %s"  % neighbor_params['npcs']
    if neighbor_params['k'] is not None:
        cmd += " --neighbors_k %s" % neighbor_params['k']
    cmd += " > logs/rna_harmony.log " 
     #job arguments
    
    if PARAMS['queues_long'] is not None:
        job_kwargs["job_queue"] = job_queue=PARAMS['queues_long']
    job_kwargs["job_threads"] = PARAMS['resources_threads_high']
    P.run(cmd, **job_kwargs)

# rna SCANORAMA
@follows(set_up_dirs)
@active_if(PARAMS['rna_run'])
@active_if(PARAMS['rna_tools'] is not None and 'scanorama' in PARAMS['rna_tools'])
@originate("batch_correction/umap_rna_scanorama.csv")
def run_scanorama(outfile):
    cmd = """python %(py_path)s/batch_correct_scanorama.py 
     --input_anndata %(preprocessed_obj)s
     --output_csv %(outfile)s 
     --integration_col %(rna_column)s
     --n_threads %(resources_threads_high)s
     --modality rna
     """
    # cannot use the normal method for importing params from yaml, because it only works up to depth 2
    neighbor_params = PARAMS['rna']['neighbors']
    if neighbor_params['method'] is not None:
        cmd += " --neighbors_method %s" % neighbor_params['method']
    if neighbor_params['metric'] is not None:
        cmd += " --neighbors_metric %s" % neighbor_params['metric']
    if neighbor_params['npcs'] is not None:
        cmd += " --neighbors_n_pcs %s"  % neighbor_params['npcs']
    if neighbor_params['k'] is not None:
        cmd += " --neighbors_k %s" % neighbor_params['k']
    cmd += " > logs/rna_scanorama.log " 
    #job arguments
    
    if PARAMS['queues_long'] is not None:
        job_kwargs["job_queue"] = job_queue=PARAMS['queues_long']
    job_kwargs["job_threads"] = PARAMS['resources_threads_high']
    P.run(cmd, **job_kwargs)

# rna scvi
@follows(set_up_dirs)
@active_if(PARAMS['rna_run'])
@active_if(PARAMS['rna_tools'] is not None and 'scvi' in PARAMS['rna_tools'])
@originate("batch_correction/umap_rna_scvi.csv")
def run_scvi(outfile):
    cmd = """python %(py_path)s/batch_correct_scvi.py 
     --scaled_anndata %(preprocessed_obj)s
     --output_csv %(outfile)s 
     --integration_col %(rna_column)s
     --figdir figures/rna/
     """
    # cannot use the normal method for importing params from yaml, because it only works up to depth 2
    neighbor_params = PARAMS['rna']['neighbors']
    if neighbor_params['method'] is not None:
        cmd += " --neighbors_method %s" % neighbor_params['method']
    if neighbor_params['metric'] is not None:
        cmd += " --neighbors_metric %s" % neighbor_params['metric']
    if neighbor_params['npcs'] is not None:
        cmd += " --neighbors_n_pcs %s"  % neighbor_params['npcs']
    if neighbor_params['k'] is not None:
        cmd += " --neighbors_k %s" % neighbor_params['k']
    cmd += " > logs/rna_scvi.log "
    job_kwargs = {}
    if PARAMS['queues_gpu'] is not None:
        job_kwargs["job_queue"] = PARAMS['queues_gpu']
        job_kwargs["job_threads"] = int(PARAMS['resources_threads_gpu'])
    elif PARAMS['queues_long'] is not None:
            job_kwargs["job_queue"] = job_queue=PARAMS['queues_long']
            job_kwargs["job_threads"] = int(PARAMS['resources_threads_high'])
    else:
        job_kwargs["job_threads"] = int(PARAMS['resources_threads_high'])
    P.run(cmd, **job_kwargs)

# rna scib
# @transform([run_scvi, run_scanorama, run_bbknn, run_harmony, run_combat],
#             regex(r"batch_correction/umap_rna_(.*).csv"),
#             r'batch_correction/scib_metrics_\1.csv',r'\1', 
#             scaled_file)
# def run_scib(infiles, outfile, method,  uncorrected_anndata):
#     print(infiles, outfile)
#     log_file = "logs/rna_scib_" + method + ".log"
#     # this is going to be apain since I want the temp h5ad files not the csvs
#     # get batch method from infile
#     corrected_anndata = "tmp/" + method + "_scaled_adata.h5ad"
#     cmd = """python %(py_path)s/run_scib_metrics.py 
#     --uncorrected_anndata %(uncorrected_anndata)s \
#     --batch_corrected_anndata %(corrected_anndata)s \
#     --integration_method %(method)s \
#     --integration_col %(rna_column)s \
#     --outfile %(outfile)s  
#     """
#     if PARAMS['use_muon']:
#         cmd += " --use_muon True"
#     if PARAMS['reference_col'] is not None:
#         cmd += " --reference_col %(scib_reference_col)s "
#     if PARAMS['scib_rough_clust'] is not None:
#         cmd += " --rough_consensus_clustering  %(scib_rough_clust)s"

#     cmd += " > %(log_file)s "
#     print(cmd)
#     P.run(cmd, job_threads=PARAMS['resources_threads_high'])
@active_if(PARAMS['rna_run'])
@follows(run_harmony, run_bbknn_rna, run_combat, run_scanorama, run_scvi, run_no_batch_umap)
def run_unimodal_integration_rna():
    pass

# ------------------------------------------------------------------------
# Prot unimodal Integration:
# - nobatch
# - Harmony
# ------------------------------------------------------------------------

@follows(set_up_dirs)
@active_if(PARAMS['prot_run'])
@originate("batch_correction/umap_prot_none.csv")
def run_no_batch_umap_prot(outfile):
    cmd = """python %(py_path)s/batch_correct_none.py 
     --input_anndata %(preprocessed_obj)s
     --output_csv %(outfile)s
     --integration_col %(prot_column)s
     """
    cmd += " --modality prot"
    # cannot use the normal method for importing params from yaml, because it only works up to depth 2
    neighbor_params = PARAMS['prot']['neighbors']
    if neighbor_params['method'] is not None:
        cmd += " --neighbors_method %s" % neighbor_params['method']
    if neighbor_params['metric'] is not None:
        cmd += " --neighbors_metric %s" % neighbor_params['metric']
    if neighbor_params['npcs'] is not None:
        cmd += " --neighbors_n_pcs %s"  % neighbor_params['npcs']
    if neighbor_params['k'] is not None:
        cmd += " --neighbors_k %s" % neighbor_params['k']
    cmd += " > logs/prot_no_correct.log"
    
    if PARAMS['queues_long'] is not None:
        job_kwargs["job_queue"] = job_queue=PARAMS['queues_long']
    job_kwargs["job_threads"] = PARAMS['resources_threads_high']
    P.run(cmd, **job_kwargs) 

# prot HARMONY
@follows(set_up_dirs)
@active_if(PARAMS['prot_run'])
@active_if(PARAMS['prot_tools'] is not None and 'harmony' in PARAMS['prot_tools'])
@originate("batch_correction/umap_prot_harmony.csv")
def run_harmony_prot( outfile):
    cmd = """python %(py_path)s/batch_correct_harmony.py 
     --input_anndata %(preprocessed_obj)s
     --output_csv %(outfile)s 
     --integration_col %(prot_column)s
     --n_threads %(resources_threads_high)s
     """
    cmd += " --modality prot"
    harmony_params = PARAMS['prot']['harmony']
    if harmony_params['npcs'] is not None:
        cmd += " --harmony_npcs %s" % harmony_params['npcs']
    if harmony_params['sigma'] is not None:
        cmd += " --sigma_val %s" % harmony_params['sigma'] 
    if harmony_params['theta'] is not None:
        cmd += " --theta_val %s" % harmony_params['theta']   
    neighbor_params = PARAMS['prot']['neighbors']
    if neighbor_params['method'] is not None:
        cmd += " --neighbors_method %s" % neighbor_params['method']
    if neighbor_params['metric'] is not None:
        cmd += " --neighbors_metric %s" % neighbor_params['metric']
    if neighbor_params['npcs'] is not None:
        cmd += " --neighbors_n_pcs %s"  % neighbor_params['npcs']
    if neighbor_params['k'] is not None:
        cmd += " --neighbors_k %s" % neighbor_params['k']
    cmd += " > logs/prot_harmony.log " 
     #job arguments
    
    if PARAMS['queues_long'] is not None:
        job_kwargs["job_queue"] = job_queue=PARAMS['queues_long']
    job_kwargs["job_threads"] = PARAMS['resources_threads_high']
    P.run(cmd, **job_kwargs)



# prot BBKNN
@follows(set_up_dirs)
@active_if(PARAMS['prot_run'])
@active_if(PARAMS['prot_tools'] is not None and 'bbknn' in PARAMS['prot_tools'])
# @follows(normalise_log_hvg_regress_scale)
@originate("batch_correction/umap_prot_bbknn.csv")
def run_bbknn_prot(outfile):
    cmd = """python %(py_path)s/batch_correct_bbknn.py 
     --input_anndata %(preprocessed_obj)s
     --output_csv %(outfile)s
     --integration_col %(prot_column)s
     --modality prot
     """
    if PARAMS['prot']['bbknn']['neighbors_within_batch'] is not None:
        cmd += " --neighbors_within_batch %i" % PARAMS['prot']['bbknn']['neighbors_within_batch']
    if PARAMS['prot']['neighbors']['npcs'] is not None:
        cmd += " --neighbors_n_pcs %s" % PARAMS['prot']['neighbors']['npcs']
    cmd += " > logs/prot_bbknn.log "
    
    if PARAMS['queues_long'] is not None:
        job_kwargs["job_queue"] = job_queue=PARAMS['queues_long']
    job_kwargs["job_threads"] = PARAMS['resources_threads_high']
    P.run(cmd, **job_kwargs) 


# prot COMBAT
@follows(set_up_dirs)
@active_if(PARAMS['prot_run'])
@active_if(PARAMS['prot_tools'] is not None and 'combat' in PARAMS['prot_tools'])
@originate("batch_correction/umap_prot_combat.csv")
def run_combat_prot(outfile):
    cmd = """python %(py_path)s/batch_correct_combat.py 
     --input_anndata %(preprocessed_obj)s
     --output_csv %(outfile)s
     --integration_col %(prot_column)s
     --n_threads %(resources_threads_high)s
     --modality prot
     """
    # cannot use the normal method for importing params from yaml, because it only works up to depth 2
    neighbor_params = PARAMS['prot']['neighbors']
    if neighbor_params['method'] is not None:
        cmd += " --neighbors_method %s" % neighbor_params['method']
    if neighbor_params['metric'] is not None:
        cmd += " --neighbors_metric %s" % neighbor_params['metric']
    if neighbor_params['npcs'] is not None:
        cmd += " --neighbors_n_pcs %s"  % neighbor_params['npcs']
    if neighbor_params['k'] is not None:
        cmd += " --neighbors_k %s" % neighbor_params['k']
    cmd += " > logs/prot_combat.log "
    
    if PARAMS['queues_long'] is not None:
        job_kwargs["job_queue"] = PARAMS['queues_long']

    job_kwargs["job_threads"] = PARAMS['resources_threads_high']
    P.run(cmd, **job_kwargs) 
    

@active_if(PARAMS['prot_run'])
@follows(run_harmony_prot, run_bbknn_prot,run_combat_prot, run_no_batch_umap_prot)
def run_unimodal_integration_prot():
    pass


# ------------------------------------------------------------------------
# ATAC unimodal Integration:
# - nobatch
# - Harmony
# - BBKNN
# ------------------------------------------------------------------------
@follows(set_up_dirs)
@active_if(PARAMS['atac_run'])
@originate("batch_correction/umap_atac_none.csv")
def run_no_batch_umap_atac(outfile):
    cmd = """python %(py_path)s/batch_correct_none.py 
     --input_anndata %(preprocessed_obj)s
     --output_csv %(outfile)s
     --integration_col %(atac_column)s
     """
    cmd += " --modality atac"
    neighbor_params = PARAMS['atac']['neighbors']
    if neighbor_params['method'] is not None:
        cmd += " --neighbors_method %s" % neighbor_params['method']
    if neighbor_params['metric'] is not None:
        cmd += " --neighbors_metric %s" % neighbor_params['metric']
    if neighbor_params['npcs'] is not None:
        cmd += " --neighbors_n_pcs %s"  % neighbor_params['npcs']
    if neighbor_params['k'] is not None:
        cmd += " --neighbors_k %s" % neighbor_params['k']
    atac_dimred = PARAMS['atac']['dimred']
    if PARAMS['atac']['dimred'] is not None:
        cmd += " --dimred %s" % atac_dimred
    else:
        cmd += " --dimred PCA"
    cmd += " > logs/atac_no_correct.log"
    
    if PARAMS['queues_long'] is not None:
        job_kwargs["job_queue"] = job_queue=PARAMS['queues_long']
    job_kwargs["job_threads"] = PARAMS['resources_threads_high']
    P.run(cmd, **job_kwargs) 

# atac HARMONY
@follows(set_up_dirs)
@active_if(PARAMS['atac_run'])
@active_if(PARAMS['atac_tools'] is not None and 'harmony' in PARAMS['atac_tools'])
@originate("batch_correction/umap_atac_harmony.csv")
def run_harmony_atac( outfile):
    cmd = """python %(py_path)s/batch_correct_harmony.py 
     --input_anndata %(preprocessed_obj)s
     --output_csv %(outfile)s 
     --integration_col %(atac_column)s
     --n_threads %(resources_threads_high)s
     --modality atac
     """

    harmony_params = PARAMS['atac']['harmony']
    if harmony_params['npcs'] is not None:
        cmd += " --harmony_npcs %s" % harmony_params['npcs']
    if harmony_params['sigma'] is not None:
        cmd += " --sigma_val %s" % harmony_params['sigma'] 
    if harmony_params['theta'] is not None:
        cmd += " --theta_val %s" % harmony_params['theta']   
    neighbor_params = PARAMS['atac']['neighbors']
    if neighbor_params['method'] is not None:
        cmd += " --neighbors_method %s" % neighbor_params['method']
    if neighbor_params['metric'] is not None:
        cmd += " --neighbors_metric %s" % neighbor_params['metric']
    if neighbor_params['npcs'] is not None:
        cmd += " --neighbors_n_pcs %s"  % neighbor_params['npcs']
    if neighbor_params['k'] is not None:
        cmd += " --neighbors_k %s" % neighbor_params['k']
    atac_dimred = PARAMS['atac']['dimred']
    if PARAMS['atac']['dimred'] is not None:
        cmd += " --dimred %s" % atac_dimred
    else:
        cmd += " --dimred PCA"
    cmd += " > logs/atac_harmony.log " 
     #job arguments
    
    if PARAMS['queues_long'] is not None:
        job_kwargs["job_queue"] = job_queue=PARAMS['queues_long']
    job_kwargs["job_threads"] = PARAMS['resources_threads_high']
    P.run(cmd, **job_kwargs)

# atac BBKNN
@follows(set_up_dirs)
@active_if(PARAMS['atac_run'] )
@active_if(PARAMS['atac_tools'] is not None and 'bbknn' in PARAMS['atac_tools'])
# @follows(normalise_log_hvg_regress_scale)
@originate("batch_correction/umap_atac_bbknn.csv")
def run_bbknn_atac(outfile):
    cmd = """python %(py_path)s/batch_correct_bbknn.py 
     --input_anndata %(preprocessed_obj)s
     --output_csv %(outfile)s
     --integration_col %(atac_column)s
     --modality atac
     """
    if PARAMS['atac']['bbknn']['neighbors_within_batch'] is not None:
        cmd += " --neighbors_within_batch %i" % PARAMS['atac']['bbknn']['neighbors_within_batch']
    if PARAMS['atac']['neighbors']['npcs'] is not None:
        cmd += " --neighbors_n_pcs %s" % PARAMS['atac']['neighbors']['npcs']
    #Forcing bbknn to run on PCA in case of atac
    cmd += " --dimred PCA"
    cmd += " > logs/atac_bbknn.log "
    if PARAMS['queues_long'] is not None:
        job_kwargs["job_queue"] = PARAMS['queues_long']
    job_kwargs["job_threads"] = PARAMS['resources_threads_high']
    P.run(cmd, **job_kwargs) 


@active_if(PARAMS['atac_run'])
@follows(run_harmony_atac, run_no_batch_umap_atac, run_bbknn_atac)
def run_unimodal_integration_atac():
    pass

@follows(run_unimodal_integration_rna,
         run_unimodal_integration_prot,
         run_unimodal_integration_atac)
def run_unimodal_integration():
    pass
# ------------------------------------------------------------------------
# Multimodal  Integration
# ------------------------------------------------------------------------
# RUN totalvi
@follows(set_up_dirs)
@active_if(PARAMS['multimodal_run'])
@active_if(PARAMS['multimodal_tools'] is not None)
@active_if('totalvi' in [x.lower() for x in PARAMS['multimodal_tools']])
@originate("batch_correction/umap_multimodal_totalvi.csv")
def run_totalvi(outfile):
    if os.path.exists("figures/multimodal") is False:
        os.makedirs("figures/multimodal")
    cmd = """python %(py_path)s/batch_correct_totalvi.py 
     --scaled_anndata %(preprocessed_obj)s
     --output_csv %(outfile)s 
     --figdir figures/
     """
    if PARAMS['multimodal_column_categorical'] is not None:
        cmd += "--integration_col_categorical %(multimodal_column_categorical)s "
    neighbor_params = PARAMS['multimodal']['neighbors']
    if neighbor_params['method'] is not None:
        cmd += " --neighbors_method %s" % neighbor_params['method']
    if neighbor_params['metric'] is not None:
        cmd += " --neighbors_metric %s" % neighbor_params['metric']
    if neighbor_params['npcs'] is not None:
        cmd += " --neighbors_n_pcs %s"  % neighbor_params['npcs']
    if neighbor_params['k'] is not None:
        cmd += " --neighbors_k %s" % neighbor_params['k']
    cmd += " > logs/multimodal_totalvi.log "
    
    if PARAMS['queues_gpu'] is not None:
        job_kwargs["job_queue"] = PARAMS['queues_gpu']
        job_kwargs["job_threads"] = int(PARAMS['resources_threads_gpu'])
    else:
        if PARAMS['queues_long'] is not None:
            job_kwargs["job_queue"] = job_queue=PARAMS['queues_long']
        job_kwargs["job_threads"] = int(PARAMS['resources_threads_high'])
    P.run(cmd, **job_kwargs)

# Run MultiVI

@follows(set_up_dirs)
@active_if(PARAMS['multimodal_run'])
@active_if( PARAMS['multimodal_tools'] is not None)
@active_if('multivi' in [x.lower() for x in PARAMS['multimodal_tools']])
@originate("batch_correction/umap_multimodal_multivi.csv")
def run_multivi(outfile):
    if os.path.exists("figures/multimodal") is False:
        os.makedirs("figures/multimodal")
    cmd = """python %(py_path)s/batch_correct_multivi.py 
     --scaled_anndata %(preprocessed_obj)s
     --output_csv %(outfile)s 
     --figdir figures/
     """
    
    if PARAMS['multimodal_column_categorical'] is not None:
        cmd += "--integration_col_categorical %(multimodal_column_categorical)s "

    if PARAMS['multimodal_column_continuous'] is not None:
        cmd += "--integration_col_continuous %(multimodal_column_continuous)s "

    neighbor_params = PARAMS['multimodal']['neighbors']
    if neighbor_params['method'] is not None:
        cmd += " --neighbors_method %s" % neighbor_params['method']
    if neighbor_params['metric'] is not None:
        cmd += " --neighbors_metric %s" % neighbor_params['metric']
    if neighbor_params['npcs'] is not None:
        cmd += " --neighbors_n_pcs %s"  % neighbor_params['npcs']
    if neighbor_params['k'] is not None:
        cmd += " --neighbors_k %s" % neighbor_params['k']
    cmd += " > logs/multimodal_multivi.log "
    
    if PARAMS['queues_gpu'] is not None:
        job_kwargs["job_queue"] = PARAMS['queues_gpu']
        job_kwargs["job_threads"] = int(PARAMS['resources_threads_gpu'])
    else:
        if PARAMS['queues_long'] is not None:
            job_kwargs["job_queue"] = job_queue=PARAMS['queues_long']
        job_kwargs["job_threads"] = int(PARAMS['resources_threads_high'])
    P.run(cmd, **job_kwargs)



# Run Mofa
@follows(set_up_dirs)
@active_if(PARAMS['multimodal_run'])
@active_if( PARAMS['multimodal_tools'] is not None and 'mofa' in PARAMS['multimodal_tools'])
@originate("batch_correction/umap_multimodal_mofa.csv")
def run_mofa(outfile):
    if os.path.exists("figures/multimodal") is False:
        os.makedirs("figures/multimodal")
    cmd = """python %(py_path)s/batch_correct_mofa.py 
     --scaled_anndata %(preprocessed_obj)s
     --output_csv %(outfile)s 
     --figdir figures/
     """
    if PARAMS['multimodal_column_categorical'] is not None:
        cmd += " --integration_col_categorical %(multimodal_column_categorical)s"
    mofa_params = PARAMS['multimodal']['mofa']
    if mofa_params['n_factors'] is not None:
        cmd += " --n_factors %s" % mofa_params['n_factors']
    if mofa_params['n_iterations'] is not None:
        cmd += " --n_iterations %s" % mofa_params['n_iterations']
    if mofa_params['convergence_mode'] is not None:
        cmd += " --convergence_mode %s" % mofa_params['convergence_mode']
    if mofa_params['save_parameters']:
        cmd += " --save_parameters %s" % mofa_params['save_parameters']
        cmd += " --outfile_model %s" % mofa_params['outfile']
    neighbor_params = PARAMS['multimodal']['neighbors']
    if neighbor_params['method'] is not None:
        cmd += " --neighbors_method %s" % neighbor_params['method']
    if neighbor_params['metric'] is not None:
        cmd += " --neighbors_metric %s" % neighbor_params['metric']
    if neighbor_params['npcs'] is not None:
        cmd += " --neighbors_n_pcs %s"  % neighbor_params['npcs']
    if neighbor_params['k'] is not None:
        cmd += " --neighbors_k %s" % neighbor_params['k']
    
    if PARAMS['queues_gpu'] is not None:
        cmd += " --use_gpu True"
        job_kwargs["job_queue"] = PARAMS['queues_gpu']
        job_kwargs["job_threads"] = int(PARAMS['resources_threads_gpu'])
    else:
        if PARAMS['queues_long'] is not None:
            job_kwargs["job_queue"] = job_queue=PARAMS['queues_long']
        job_kwargs["job_threads"] = int(PARAMS['resources_threads_high'])
    cmd += " > logs/multimodal_mofa.log "
    P.run(cmd, **job_kwargs)

# Run WNN

@follows(set_up_dirs)
@follows(run_unimodal_integration)
@active_if(PARAMS['multimodal_run'])
@active_if( PARAMS['multimodal_tools'] is not None)
@active_if('wnn' in [x.lower() for x in PARAMS['multimodal_tools']])
@originate("batch_correction/umap_multimodal_wnn.csv")
def run_wnn(outfile):
    if os.path.exists("figures/multimodal") is False:
        os.makedirs("figures/multimodal")
    cmd = """python %(py_path)s/batch_correct_wnn.py 
     --scaled_anndata %(preprocessed_obj)s
     --output_csv %(outfile)s 
     --figdir figures/
     """
    wnn_params = PARAMS['multimodal']['WNN']
    if wnn_params['n_neighbors'] is not None:
       cmd += " --n_neighbors %s" % wnn_params['n_neighbors']
    if wnn_params['n_bandwidth_neighbors'] is not None:
       cmd += " --n_bandwidth_neighbors %s" % wnn_params['n_bandwidth_neighbors']
    if wnn_params['n_multineighbors'] is not None:
       cmd += " --n_multineighbors %s" % wnn_params['n_multineighbors']
    if wnn_params['metric'] is not None:
       cmd += " --metric %s" % wnn_params['metric']
    if wnn_params['low_memory'] is not None:
       cmd += " --low_memory %s" % wnn_params['low_memory']
    cmd += " > logs/multimodal_wnn.log"

    
    if PARAMS['queues_gpu'] is not None:
        job_kwargs["job_queue"] = PARAMS['queues_gpu']
        job_kwargs["job_threads"] = int(PARAMS['resources_threads_gpu'])
    else:
        if PARAMS['queues_long'] is not None:
            job_kwargs["job_queue"] = job_queue=PARAMS['queues_long']
        job_kwargs["job_threads"] = int(PARAMS['resources_threads_high'])
        
    P.run(cmd, **job_kwargs)

# end of multimodal
@follows(run_mofa, run_wnn, run_totalvi, run_multivi)
def run_multimodal_integration():
    pass
# ------------------------------------------------------------------------
# Evaluation
# ------------------------------------------------------------------------

@collate([run_no_batch_umap, run_scanorama, run_bbknn_rna,run_harmony, run_scvi, run_combat, 
          run_no_batch_umap_prot, run_harmony_prot, run_bbknn_prot, run_combat_prot,
          run_no_batch_umap_atac, run_bbknn_atac,run_harmony_atac, 
          run_totalvi, run_wnn, run_multivi, run_mofa], 
         regex(r"(.*)/(.*)"), 
          r'batch_correction/combined_umaps.tsv')
def collate_integration_outputs(infiles,outfile):
    #infiles_string = ','.join(infiles)
    contents = glob.glob("batch_correction/*.csv")
    union = list(set(infiles) | set(contents))
    infiles_string = ','.join(union)
    batch_dict =  "batch_correction/batch_dict.yml"
    cell_mtd_file = sprefix + "_cell_mtd.csv"
    cmd = """python %(py_path)s/run_collate_mtd_files.py 
    --input_mudata %(preprocessed_obj)s
    --input_umap_files %(infiles_string)s 
    --output_batch_yml %(batch_dict)s 
    --output_cell_metadata_csv %(cell_mtd_file)s
    --output_combined_umaps_tsv %(outfile)s
    """ 
        # add the integration cols
    if PARAMS['rna_run']:
        cmd += " --rna_integration_col %(rna_column)s"
    if PARAMS['prot_run']:
        cmd += " --prot_integration_col %(prot_column)s"
    if PARAMS['atac_run']:
        cmd += " --atac_integration_col %(atac_column)s"
    if PARAMS['multimodal_run']:
        cmd += " --multimodal_integration_col %(multimodal_column_categorical)s"
    cmd += " > logs/collate_mtd.log "
    job_kwargs["job_threads"] = PARAMS['resources_threads_medium']
    P.run(cmd, **job_kwargs)  


# this COLLATE job will become the big final collate job across all possible combinations you may have run
@follows(set_up_dirs)
@transform(collate_integration_outputs, formatter(), "logs/plot_batch_corrected_umaps.log")
def plot_umaps(infile, outfile):
    print(infile, outfile)
    cell_mtd_file = sprefix + "_cell_mtd.csv"
    cmd = """python %(py_path)s/plot_umaps_batch_correct.py 
    --fig_dir figures/
    --combined_umaps_tsv %(infile)s
    --batch_dict batch_correction/batch_dict.yml
    --cell_meta_df %(cell_mtd_file)s
    """
    # 
    # Parse the whole dict of metrics
    # this contains grouping var and qc_metrics per modality
    cmd += ' --qc_dict "%(plotqc)s"'
    cmd += " > %(outfile)s" 
    job_kwargs["job_threads"] = PARAMS['resources_threads_medium']
    P.run(cmd,**job_kwargs)


#this can follow now any mtd generation, but it will collate only RNA jobs for lisi
@follows(collate_integration_outputs)
@transform(collate_integration_outputs, 
           formatter(),  'logs/lisi.log')
def run_lisi(infile, outfile):
    cell_mtd_file = sprefix + "_cell_mtd.csv"  
    cmd = """python %(py_path)s/run_lisi.py 
    --combined_umaps_df %(infile)s 
    --cell_meta_df %(cell_mtd_file)s
    --integration_dict batch_correction/batch_dict.yml
    --fig_dir figures/  > %(outfile)s 
    """ 
    job_kwargs["job_threads"] = PARAMS['resources_threads_low']
    P.run(cmd,**job_kwargs)



@follows(run_unimodal_integration, run_multimodal_integration,run_lisi, plot_umaps)
@originate("logs/batch_correction_complete.log")
def batch_correction(outfile):
    IOTools.touch_file(outfile)





@follows( batch_correction)
def full():
    """
    All cgat pipelines should end with a full() function which updates,
    if needed, all branches of the pipeline.
    The @follows statement should ensure that all functions are covered,
    either directly or as prerequisites.
    """
    pass

# Final choices: Leave blank until you have reviewed the results from running
# sc_pipelines integration make full
# choose the parameters to integrate into the final anndata object
# then run
# `panpipes integration make merge_integration`
@originate(sprefix + "_corrected.h5mu")
def merge_integration(output_mudata):
    cmd = """
    python %(py_path)s/batch_correct_merge.py 
    --preprocessed_mudata %(preprocessed_obj)s
    --output_mudata %(output_mudata)s
    """
    choices = PARAMS['final_obj']
    if choices['rna']['include'] and choices['rna']['bc_choice'] is not None :
        rna_choice = choices['rna']['bc_choice']
        cmd += " --rna_correction_choice %(rna_choice)s"
    if choices['prot']['include'] and choices['prot']['bc_choice'] is not None:
        prot_choice = choices['prot']['bc_choice']
        cmd += " --prot_correction_choice %(prot_choice)s"
    if choices['atac']['include'] and choices['atac']['bc_choice'] is not None:
        atac_choice = choices['atac']['bc_choice']
        cmd += " --atac_correction_choice %(atac_choice)s"
    if choices['multimodal']['include'] and choices['multimodal']['bc_choice'] is not None:
        multimodal_choice = choices['multimodal']['bc_choice']
        cmd += " --multimodal_correction_choice %(multimodal_choice)s"
    cmd += " > logs/merge_final_obj.log"
    job_kwargs["job_threads"] = PARAMS['resources_threads_high']
    P.run(cmd,**job_kwargs)
    # clear up tmp since it is no longer required
    # P.run("rm -r tmp")


def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
