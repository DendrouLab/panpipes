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
#{'rna': True, 'prot': False, 'Repertoire': False, 'atac': False}

# -----------------------------------------------------------------------------------------------
## 10X metrics plotting
# -----------------------------------------------------------------------------------------------

# TODO: update this to check for each modality csv file
# TODO: check this works with current submission file
@active_if(PARAMS['plot_10X_metrics'])
@follows(mkdir('figures'))
@follows(mkdir('figures/tenx_metrics'))
@follows(mkdir("logs"))
@originate("logs/tenx_metrics_multi_aggregate.log")
def aggregate_tenx_metrics_multi(outfile):
    """this is to aggregate all the cellranger multi metric_summary files
    it also does some plotting
    """    
    cmd = """
        python %(py_path)s/aggregate_cellranger_summary_metrics.py
            --pipe_df %(submission_file)s
            --figdir figures/tenx_metrics/
            --cellranger_column_conversion_df %(resources_path)s/metrics_summary_col_conversion.tsv
            --output_file 10x_metrics.csv > %(outfile)s
            """
    job_kwargs["job_threads"] = PARAMS['resources_threads_low']
    P.run(cmd, **job_kwargs)


@follows(aggregate_tenx_metrics_multi)
def process_all_tenx_metrics():
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
    return gen_load_anndata_jobs(caf, load_raw=False, mode_dictionary=PARAMS["modalities"], load_prot_from_raw=PARAMS['load_prot_from_raw'])

    

@follows(mkdir("logs"))
@follows(mkdir("tmp"))
@active_if(PARAMS["use_existing_h5mu"] is False)
@files(gen_load_filtered_anndata_jobs)
def load_mudatas(rna_path, outfile, 
              sample_id, 
              rna_filetype,  
              prot_path, prot_filetype, 
              tcr_path, tcr_filetype,  
              bcr_path, bcr_filetype, 
              atac_path, atac_filetype, 
              fragments_file, per_barcode_metrics_file, peak_annotation_file, 
              cell_mtd_path):
    
    path_dict = {'rna':rna_path,
                 'prot': prot_path,
                 'bcr': bcr_path,
                 'tcr' : tcr_path,
                 'atac':atac_path}
    print(path_dict)
    
    modality_dict = {k:True if path_dict[k] is not None else False for k,v in PARAMS['modalities'].items() }
    print(modality_dict)
    
    cmd = """
        python %(py_path)s/make_adata_from_csv.py 
        --mode_dictionary "%(modality_dict)s"
        --sample_id %(sample_id)s
        --output_file %(outfile)s 
    """
    
    cmd += " --use_muon True" #this is still here?
    if rna_path is not None and pd.notna(rna_path):
        cmd += " --rna_infile %(rna_path)s"
        cmd += " --rna_filetype %(rna_filetype)s"
    if prot_path is not None and pd.notna(prot_path):
        cmd += " --prot_infile %(prot_path)s"
        cmd += " --prot_filetype %(prot_filetype)s"
        cmd += " --subset_prot_barcodes_to_rna %(subset_prot_barcodes_to_rna)s"
    if atac_path is not None and pd.notna(atac_path):
        cmd += " --atac_infile %(atac_path)s"
        cmd += " --atac_filetype %(atac_filetype)s"
    if per_barcode_metrics_file is not None and pd.notna(per_barcode_metrics_file):
        cmd += " --per_barcode_metrics_file %(per_barcode_metrics_file)s"
    if fragments_file is not None and pd.notna(fragments_file):
        cmd += " --fragments_file %(fragments_file)s"
    if peak_annotation_file is not None and pd.notna(peak_annotation_file):
        cmd += " --peak_annotation_file %(peak_annotation_file)s"
      # ~ means this tests "is not nan"
    if tcr_path is not None and pd.notna(tcr_path):
        cmd += " --tcr_filtered_contigs %(tcr_path)s"
        cmd += " --tcr_filetype %(tcr_filetype)s"
    if bcr_path is not None and pd.notna(bcr_path):
        cmd += " --bcr_filtered_contigs %(bcr_path)s"
        cmd += " --bcr_filetype %(bcr_filetype)s"
    cmd += " > logs/load_mudatas_%(sample_id)s.log"
    # print(cmd)
    job_kwargs["job_threads"] = PARAMS['resources_threads_medium']
    P.run(cmd, **job_kwargs)



@active_if(PARAMS["use_existing_h5mu"] is False)
@collate(load_mudatas,
         formatter(""),
         unfilt_file())
def concat_filtered_mudatas(infiles, outfile):
    # print(infiles)
    # print(outfile)
    infiles_str = ','.join(infiles)
    cmd = """
    python %(py_path)s/concat_adata.py 
        --input_files_str %(infiles_str)s 
        --output_file %(outfile)s 
        --submissionfile %(submission_file)s
        --sampleprefix %(sample_prefix)s
        --join_type %(concat_join_type)s
        """
    if PARAMS['metadatacols'] is not None and PARAMS['metadatacols'] != "":
        cmd += " --metadatacols  %(metadatacols)s"
    if PARAMS["barcode_mtd_include"] is True:
        cmd += " --barcode_mtd_df %(barcode_mtd_path)s"
        cmd += " --barcode_mtd_metadatacols %(barcode_mtd_metadatacols)s"
    if PARAMS['protein_metadata_table'] is not None:
        cmd += " --protein_var_table %(protein_metadata_table)s"
    if PARAMS['index_col_choice'] is not None:
        cmd += " --protein_new_index_col %(index_col_choice)s"
    cmd += " > logs/concat_filtered_mudatas.log"
    job_kwargs["job_threads"] = PARAMS['resources_threads_high']
    P.run(cmd, **job_kwargs)
    # P.run("rm tmp/*", job_threads=PARAMS['resources_threads_low'])


# -----------------------------------------------------------------------------------------------
## Creating h5mu from bg data files
# -----------------------------------------------------------------------------------------------

if PARAMS['assess_background'] or (PARAMS["modalities"]["prot"] and "dsb" in PARAMS['normalisation_methods']):
    PARAMS['bg_required'] = True
else:
    PARAMS['bg_required'] = False
    

def bg_file():
    sprefix = PARAMS['sample_prefix']
    use_muon = True
    if use_muon:
        bg_file = sprefix + "_bg.h5mu"
    else:
        bg_file = sprefix + "_bg.h5ad"
    return bg_file


def gen_load_bg_anndata_jobs():
    caf = pd.read_csv(PARAMS['submission_file'], sep='\t')
    return gen_load_anndata_jobs(caf, load_raw=True, mode_dictionary=PARAMS["modalities"], load_prot_from_raw=True)
    

@active_if(PARAMS["bg_required"])
@follows(mkdir("logs"))
@follows(mkdir("tmp"))
@active_if(PARAMS["use_existing_h5mu"] is False)
@files(gen_load_bg_anndata_jobs)
def load_bg_mudatas(rna_path, outfile, 
              sample_id, 
              rna_filetype,  
              prot_path, prot_filetype, 
              tcr_path, tcr_filetype,  
              bcr_path, bcr_filetype, 
              atac_path, atac_filetype, 
              fragments_file, per_barcode_metrics_file, peak_annotation_file, 
              cell_mtd_path):
    path_dict = {'rna':rna_path,
                 'prot': prot_path,
                 'bcr': bcr_path,
                 'tcr' : tcr_path,
                 'atac':atac_path}
    print(path_dict)
    
    # need to remove unnessacry modality
    modality_dict = {k:True if path_dict[k] is not None else False for k,v in PARAMS['modalities'].items() }

    modality_dict['tcr'] = False
    modality_dict['bcr'] = False
    cmd = """
    python %(py_path)s/make_adata_from_csv.py 
    --sample_id %(sample_id)s
    --mode_dictionary "%(modality_dict)s"
    --rna_infile %(rna_path)s
    --rna_filetype %(rna_filetype)s
    --output_file %(outfile)s 
    """
    if prot_path is not None and pd.notna(prot_path):
        cmd += " --prot_infile %(prot_path)s"
        cmd += " --prot_filetype %(prot_filetype)s"
    if atac_path is not None and pd.notna(atac_path):
        cmd += " --atac_infile %(atac_path)s"
        cmd += " --atac_filetype %(atac_filetype)s"
    if per_barcode_metrics_file is not None and pd.notna(per_barcode_metrics_file):
        cmd += " --per_barcode_metrics_file %(per_barcode_metrics_file)s"
    if fragments_file is not None and pd.notna(fragments_file):
        cmd += " --fragments_file %(fragments_file)s"
    if peak_annotation_file is not None and pd.notna(peak_annotation_file):
        cmd += " --peak_annotation_file %(peak_annotation_file)s"
    if PARAMS['protein_metadata_table'] is not None:
        cmd += " --protein_var_table %(protein_metadata_table)s"  #check which of these 2 needs to stay!!!
    if PARAMS['index_col_choice'] is not None:
        cmd += " --protein_new_index_col %(index_col_choice)s"
    cmd += " > logs/load_bg_mudatas_%(sample_id)s.log"
    job_kwargs["job_threads"] = PARAMS['resources_threads_medium']
    P.run(cmd, **job_kwargs)


@active_if(PARAMS["bg_required"])
@active_if(PARAMS['downsample_background'])
@follows(load_bg_mudatas)
@transform(load_bg_mudatas, 
           regex("./tmp/(.*)_raw.h5(.*)"), 
           r"./logs/\1_bg_downsampled.log")
def downsample_bg_mudatas(infile, outfile):
    cmd = """
    python %(py_path)s/downsample.py 
    --input_mudata %(infile)s
    --output_mudata "%(infile)s" 
    --downsample_value 20000  > %(outfile)s
    """
    job_kwargs["job_threads"] = PARAMS['resources_threads_medium']
    P.run(cmd, **job_kwargs)


@active_if(PARAMS["bg_required"])
@active_if(PARAMS["assess_background"])
@active_if(PARAMS["use_existing_h5mu"] is False)
@follows(load_bg_mudatas)
@follows(downsample_bg_mudatas)
@collate(load_bg_mudatas,
         formatter(""),
         bg_file())
def concat_bg_mudatas(infiles, outfile):
    # print(infiles)
    # print(outfile)
    infiles_str = ','.join(infiles)
    cmd = """
    python %(py_path)s/concat_adata.py 
        --input_files_str %(infiles_str)s 
        --output_file %(outfile)s 
        --submissionfile %(submission_file)s
        --sampleprefix %(sample_prefix)s
        --join_type %(concat_join_type)s
        """
    if PARAMS['metadatacols'] is not None and PARAMS['metadatacols'] != "":
        cmd += " --metadatacols  %(metadatacols)s"
  #  if PARAMS["barcode_mtd_include"] is True:
   #     cmd += " --barcode_mtd_df %(barcode_mtd_path)s"
    #    cmd += " --barcode_mtd_metadatacols %(barcode_mtd_metadatacols)s"
    cmd += " > logs/concat_bg_mudatas.log"
    job_kwargs["job_threads"] = PARAMS['resources_threads_high']
    P.run(cmd, **job_kwargs)
    # P.run("rm tmp/*", job_threads=PARAMS['resources_threads_low'])

# -----------------------------------------------------------------------------------------------
## rna QC 
# -----------------------------------------------------------------------------------------------


def orfile():
    return PARAMS['sample_prefix'] + "_cell_metadata.tsv"
@active_if(PARAMS['modalities_rna'])
@active_if(PARAMS["use_existing_h5mu"] is False)
@follows(mkdir("scrublet"))
@transform(load_mudatas, regex(r"./tmp/(.*).h5(.*)"),
           r"scrublet/\1_scrublet_scores.txt", r"\1")
def run_scrublet(infile, outfile, sample_id):
    outdir = "./scrublet"
    cmd = """python %(py_path)s/run_scrublet_scores.py
        --sample_id %(sample_id)s
        --inputpath %(infile)s
        --outdir %(outdir)s
        """
    if PARAMS['scr_expected_doublet_rate'] is not None:
        cmd += " --expected_doublet_rate %(scr_expected_doublet_rate)s"
    if PARAMS['scr_sim_doublet_ratio'] is not None:
        cmd += " --sim_doublet_ratio %(scr_sim_doublet_ratio)s"
    if PARAMS['scr_n_neighbours'] is not None:
        cmd += " --n_neighbors %(scr_n_neighbours)s"
    if PARAMS['scr_min_counts'] is not None:
        cmd += " --min_counts %(scr_min_counts)s"
    if PARAMS['scr_min_cells'] is not None:
        cmd += " --min_cells %(scr_min_cells)s"
    if PARAMS['scr_min_gene_variability_pctl'] is not None:
        cmd += " --min_gene_variability_pctl %(scr_min_gene_variability_pctl)s"
    if PARAMS['scr_n_prin_comps'] is not None:
        cmd += " --n_prin_comps %(scr_n_prin_comps)s"
    if PARAMS['scr_use_thr'] is not None:
        cmd += " --use_thr %(scr_use_thr)s"
    if PARAMS['scr_call_doublets_thr'] is not None:
        cmd += " --call_doublets_thr %(scr_call_doublets_thr)s"
    cmd += " > logs/run_scrublet_" + sample_id + ".log"
    job_kwargs["job_threads"] = PARAMS['resources_threads_medium']
    P.run(cmd,**job_kwargs)
    IOTools.touch_file(outfile)


@active_if(PARAMS['modalities_rna'])
@follows(mkdir("figures"))
@follows(mkdir("figures/rna"))
@follows(concat_filtered_mudatas, run_scrublet)
@originate("logs/run_scanpy_qc_rna.log", orfile(), unfilt_file())
def run_rna_qc(log_file, outfile, unfilt_file):
    # infile = submission file
    resources_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), "resources")
    customgenesfile = PARAMS["custom_genes_file"]
    calcproportions= PARAMS["calc_proportions"]
    cmd = """
        python %(py_path)s/run_scanpyQC_rna.py
          --sampleprefix %(sample_prefix)s
          --input_anndata %(unfilt_file)s
          --outfile %(unfilt_file)s
          --customgenesfile %(customgenesfile)s
          --calc_proportions %(calcproportions)s
          --figdir figures
          """
    if len(os.listdir('scrublet')) != 0:
        cmd += " --scrubletdir scrublet"
    # if params is specified otherwise default values will be used.
    # if params is specified otherwise default values will be used.
    if PARAMS['score_genes'] is not None:
        scoregenes = PARAMS['score_genes']
        cmd += " --score_genes %(scoregenes)s"
    if PARAMS['ccgenes'] is not None:
        if PARAMS['ccgenes'] == "default":
            ccgenesfile = resources_path + "/cell_cycle_genes.tsv"
        else:
            ccgenesfile = PARAMS['ccgenes']
        cmd += " --ccgenes %(ccgenesfile)s"
    cmd += " --channel_col %(channel_col)s" # is this needed?

    # if PARAMS['isotype_upper_quantile'] is not None:
    #     cmd += " --isotype_upper_quantile %(isotype_upper_quantile)s"
    # if PARAMS['isotype_n_pass'] is not None:
    #     cmd += " --isotype_n_pass %(isotype_n_pass)s"
    # add log file
    cmd += " > %(log_file)s"
    job_kwargs["job_threads"] = PARAMS['resources_threads_high']
    P.run(cmd, **job_kwargs)
    if os.path.exists("cache"):
        P.run("rm -r cache")
    #IOTools.touch_file(metadata)

#TO DO ADD the other modalities in this order:

# prot run_scanpy_qc_prot.log


# -----------------------------------------------------------------------------------------------
## prot QC 
# -----------------------------------------------------------------------------------------------
    
# this is because we don't want to try and load the file during the `config` stage

# def prot_cell_mtd():
#     return PARAMS['sample_prefix'] + "_prot_cell_metadata.tsv"
@active_if(PARAMS['modalities_prot'])
@follows(mkdir('figures/prot'))
@follows(concat_filtered_mudatas, run_rna_qc)
@originate("logs/run_scanpy_qc_prot.log", orfile(), unfilt_file())
def run_scanpy_prot_qc(log_file, outfile, unfilt_file):
    # infile = submission file
    cmd = """
        python %(py_path)s/run_scanpyQC_prot.py
          --sampleprefix %(sample_prefix)s
          --input_anndata %(unfilt_file)s
          --outfile %(unfilt_file)s
          --figdir figures/prot/
          """
    if PARAMS['channel_col'] is not None:
        cmd += " --channel_col %(channel_col)s"
    else:
        cmd += " --channel_col sample_id"
    if PARAMS['prot_plotqc_metrics']:
        cmd += " --per_cell_metrics %(prot_plotqc_metrics)s"
    if PARAMS['prot_metrics_per_adt']:
        cmd += " --per_adt_metrics %(prot_metrics_per_adt)s"
    if PARAMS['identify_isotype_outliers']:
        cmd += " --identify_isotype_outliers %(identify_isotype_outliers)s"
        cmd += " --isotype_upper_quantile %(isotype_upper_quantile)s"
        cmd += " --isotype_n_pass %(isotype_n_pass)s"
    # add log file
    cmd += " > %(log_file)s"
    job_kwargs["job_threads"] = PARAMS['resources_threads_high']
    P.run(cmd, **job_kwargs)
    pass

@active_if(PARAMS['modalities_prot'])
@follows(run_scanpy_prot_qc, concat_filtered_mudatas, concat_bg_mudatas)
@originate("logs/run_dsb_clr.log", unfilt_file(), bg_file())
def run_dsb_clr(outfile, unfilt_file, bg_file):
    print(unfilt_file)
    # infile = mdata_unfilt
    # outfile sampleprefix + "prot_qc_per_cell.tsv"
    cmd = """
        python %(py_path)s/run_preprocess_prot.py
        --filtered_mudata %(unfilt_file)s
        --figpath ./figures/prot
        """
    if PARAMS['channel_col'] is not None:
        cmd += " --channel_col %(channel_col)s"
    if PARAMS['normalisation_methods'] is not None:
        cmd += " --normalisation_methods %(normalisation_methods)s"
    if PARAMS['quantile_clipping'] is not None:
        cmd += " --quantile_clipping %(quantile_clipping)s"
    if PARAMS['clr_margin'] is not None:
        cmd += " --clr_margin %(clr_margin)s"
    if PARAMS['save_norm_prot_mtx'] is True:
        cmd += " --save_mtx True"
    # find the bg data if available
    if os.path.exists(bg_file):
        cmd += " --bg_mudata %(bg_file)s"
    cmd += " > %(outfile)s"
    job_kwargs["job_threads"] = PARAMS['resources_threads_high']
    P.run(cmd, **job_kwargs)

@follows(run_scanpy_prot_qc, run_dsb_clr)
def run_prot_qc():
    pass


# -----------------------------------------------------------------------------------------------
## Repertoire QC  # Repertoire run_scanpy_qc_rep.log
# -----------------------------------------------------------------------------------------------
# toggle_rep_qc = ()
@active_if(PARAMS['modalities_bcr'] or PARAMS['modalities_tcr']  )
@follows(mkdir('figures/rep'))
@follows(run_rna_qc, run_prot_qc)
@originate("logs/run_scanpy_qc_rep.log", unfilt_file())
def run_repertoire_qc(logfile, unfilt_file):
    cmd = """python %(py_path)s/run_scanpyQC_rep.py
          --sampleprefix %(sample_prefix)s
          --input_mudata %(unfilt_file)s
          --output_mudata %(unfilt_file)s
          --figdir figures/rep/
          --distance_metrics "%(ir_dist)s"
          --clonotype_metrics "%(clonotype_definition)s"
          """
    cmd += " > %(logfile)s"
    job_kwargs["job_threads"] = PARAMS['resources_threads_low']
    P.run(cmd, **job_kwargs)

# -----------------------------------------------------------------------------------------------
## atac QC # run_scanpy_qc_atac.log
# -----------------------------------------------------------------------------------------------
@active_if(PARAMS['modalities_atac'])
@mkdir('figures/atac')
@follows(run_rna_qc, run_prot_qc, run_repertoire_qc)
@originate("logs/run_scanpy_qc_atac.log", orfile(), unfilt_file())
def run_atac_qc(log_file, outfile, unfilt_file):
    # if this is a multiple samples project
    # they should be run together upfront 
    # for independent atac run cellranger-atac count, then cellranger-atac aggr
    # https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/using/aggr
    # for multiome run cellranger-arc count then cellranger-arc aggr
    # https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/pipelines/latest/using/aggr
    # 
    

    plot_qc_atac_metrics = PARAMS["plotqc_atac_metrics"]
    cmd = """
        python %(py_path)s/run_scanpyQC_atac.py
          --sampleprefix %(sample_prefix)s
          --input_anndata %(unfilt_file)s
          --outfile %(unfilt_file)s
          --figdir figures/atac/  
          --atac_qc_metrics %(plot_qc_atac_metrics)s
          """
    if PARAMS["is_paired"]:
        cmd += " --is_paired True"
    else:
        cmd += " --is_paired False"
    if PARAMS["features_tss"] is not None:
        cmd += " --is_paired False"
        feat_file=PARAMS["features_tss"]
        cmd += " --tss_coordinates %(feat_file)s"
    if PARAMS["partner_rna"] is not None:
        cmd += " --is_paired False"
        prna_file = PARAMS["partner_rna"]
        cmd += " --paired_rna %(prna_file)s"
    
    cmd += " --use_muon True"

    cmd += " > %(log_file)s"
    job_kwargs["job_threads"] = PARAMS['resources_threads_low']
    P.run(cmd, **job_kwargs)

@follows(run_rna_qc, run_prot_qc, run_repertoire_qc, run_atac_qc)
def run_qc():
    pass


# -----------------------------------------------------------------------------------------------
## QC plotting
# -----------------------------------------------------------------------------------------------

@follows(run_qc)
# @transform(run_rna_qc,
#             regex(r"(.*)_cell_metadata.tsv"),
#             r"logs/plot_qc.log", )
@originate("logs/plot_qc.log", orfile())
def plot_qc(log_file, cell_file):
    qcmetrics = PARAMS['plotqc_rna_metrics']
    cmd = """
    Rscript %(r_path)s/plotQC.R 
    --prefilter TRUE
    --cell_metadata %(cell_file)s 
    --sampleprefix %(sample_prefix)s
    --groupingvar %(plotqc_grouping_var)s
    """
    cmd += " --scanpy_or_muon muon"
    if PARAMS["modalities"]['rna'] and PARAMS['plotqc_rna_metrics'] is not None:
        cmd += " --rna_qc_metrics %(plotqc_rna_metrics)s"
    if PARAMS["modalities"]['prot'] and PARAMS['plotqc_prot_metrics'] is not None:
        cmd += " --prot_qc_metrics %(plotqc_prot_metrics)s"
    if PARAMS["modalities"]['atac'] and PARAMS['plotqc_atac_metrics'] is not None:
        cmd += " --atac_qc_metrics %(plotqc_atac_metrics)s"
    if PARAMS['modalities_bcr'] or PARAMS['modalities_tcr']:
        if PARAMS['plotqc_rep_metrics'] is not None:
            pqrm = ','.join(['plotqc_rep_metrics'])
            cmd += " --rep_qc_metrics %(pqrm)s"
    cmd += " > %(log_file)s"
    job_kwargs["job_threads"] = PARAMS['resources_threads_low']
    P.run(cmd, **job_kwargs)


# ------
# assess background
# -------

@active_if(PARAMS['assess_background'])
@follows(mkdir("figures/background"))
@follows(run_qc)
@follows(concat_bg_mudatas)
@originate("logs/assess_background.log", unfilt_file(), bg_file())
def run_assess_background(log_file, unfilt_file, bg_file):
    cmd = """
    python %(py_path)s/assess_background.py
        --filtered_mudata %(unfilt_file)s
        --bg_mudata %(bg_file)s
        --channel_col %(channel_col)s
        --figpath "./figures/background/"
    """
    cmd += " > %(log_file)s"
    job_kwargs["job_threads"] = PARAMS['resources_threads_high']
    P.run(cmd, **job_kwargs)



# # ------------
# TO DO change the decorator to follow all qc plots?
@follows(plot_qc)
@follows(run_rna_qc, run_assess_background)
def all_rna_qc():
    pass

@follows(run_dsb_clr, plot_qc)
def all_prot_qc():
    pass

@follows(process_all_tenx_metrics)
@follows(all_prot_qc)
@follows(all_rna_qc)
#@originate("ruffus.check")
def full():
    """
    All cgat pipelines should end with a full() function which updates,
    if needed, all branches of the pipeline.
    The @follows statement should ensure that all functions are covered,
    either directly or as prerequisites.
    """
    #IOTools.touch_file(file)
    pass

@originate("cleanup_done.txt")
def cleanup(file):
    # remove any ctmp fails
    P.run("rm ctmp*", without_cluster=True)
    # remove tmp dir
    P.run("rm -r tmp", without_cluster=True)
    # delete empty dirs
    P.run("find ./ -empty -type d -delete", without_cluster=True)

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
