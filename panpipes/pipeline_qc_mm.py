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
import yaml


PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "../pipeline.yml",
     "pipeline.yml"])

PARAMS['py_path'] = os.path.join(os.path.dirname(os.path.dirname(__file__)), "python")
PARAMS['R_path'] = os.path.join(os.path.dirname(os.path.dirname(__file__)), "R")
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
@follows(mkdir("logs"))
@follows(mkdir("figures"))
@follows(mkdir("figures/tenx_metrics"))
@originate("logs/plot_tenx_metrics.log")
def plot_tenx_metrics(outfile):
    r_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), "R")
    cmd = """
        Rscript %(r_path)s/produce_barplot_10xmetric.v3.R
            --csvpaths %(submission_file)s
            --outdir ./
            --figdir ./figures/tenx_metrics/
            --project %(project)s
            --kneeplot %(kneeplot)s > %(outfile)s
            """
    P.run(cmd, job_threads=PARAMS['resources_threads_low'], **job_kwargs)

# -----------------------------------------------------------------------------------------------
## Creating h5mu from filtered data files
# -----------------------------------------------------------------------------------------------

def unfilt_file():
    sprefix = PARAMS['sample_prefix']
    if PARAMS['use_muon']:
        unfilt_file = sprefix + "_unfilt.h5mu"
    else:
        unfilt_file = sprefix + "_unfilt.h5ad"
    return unfilt_file



def gen_load_filtered_anndata_jobs():
    caf = pd.read_csv(PARAMS['submission_file'], sep='\t')
    return gen_load_anndata_jobs(caf, load_raw=False, mode_dictionary=PARAMS["modalities"])

    

@follows(mkdir("logs"))
@follows(mkdir("tmp"))
@active_if(PARAMS["use_existing_h5ad"] is False)
@files(gen_load_filtered_anndata_jobs)
def load_mudatas(gex_path, outfile, 
              sample_id, 
              gex_filetype,  
              adt_path, adt_filetype, 
              tcr_path, tcr_filetype,  
              bcr_path, bcr_filetype, 
              atac_path, atac_filetype, 
              fragments_file, per_barcode_metrics_file, peak_annotation_file, 
              cell_mtd_path):
    print(gex_path, outfile, \
              sample_id, \
              gex_filetype,  \
              adt_path, adt_filetype, \
              tcr_path, tcr_filetype,  \
              bcr_path, bcr_filetype, \
              atac_path, atac_filetype, \
              fragments_file, per_barcode_metrics_file, peak_annotation_file, \
              cell_mtd_path)
    cmd = """
    python %(py_path)s/make_adata_from_csv.py 
    --mode_dictionary "%(modalities)s"
    --sample_id %(sample_id)s
    --output_file %(outfile)s 
    """
    if gex_path is not None and pd.notna(gex_path):
        cmd += " --gex_infile %(gex_path)s"
        cmd += " --gex_filetype %(gex_filetype)s"
    if adt_path is not None and pd.notna(adt_path):
        cmd += " --adt_infile %(adt_path)s"
        cmd += " --adt_filetype %(adt_filetype)s"
        cmd += " --subset_adt_barcodes_to_gex %(subset_adt_barcodes_to_gex)s"
    if atac_path is not None and pd.notna(atac_path):
        cmd += " --atac_infile %(atac_path)s"
        cmd += " --atac_filetype %(atac_filetype)s"
    if per_barcode_metrics_file is not None and pd.notna(per_barcode_metrics_file):
        cmd += " --per_barcode_metrics_file %(per_barcode_metrics_file)s"
    if fragments_file is not None and pd.notna(fragments_file):
        cmd += " --fragments_file %(fragments_file)s"
    if peak_annotation_file is not None and pd.notna(peak_annotation_file):
        cmd += " --peak_annotation_file %(peak_annotation_file)s"
    if PARAMS['use_muon']:
        cmd += " --use_muon True"
    if PARAMS['protein_metadata_table'] is not None:
        cmd += " --protein_var_table %(protein_metadata_table)s"
    if PARAMS['index_col_choice'] is not None:
        cmd += " --protein_new_index_col %(index_col_choice)s"
      # ~ means this tests "is not nan"
    if tcr_path is not None and pd.notna(tcr_path):
        cmd += " --tcr_filtered_contigs %(tcr_path)s"
        cmd += " --tcr_filetype %(tcr_filetype)s"
    if bcr_path is not None and pd.notna(bcr_path):
        cmd += " --bcr_filtered_contigs %(bcr_path)s"
        cmd += " --bcr_filetype %(bcr_filetype)s"
    if PARAMS["barcode_mtd_include"] is True:
        cmd += " --barcode_mtd_df %(barcode_mtd_path)s"
        cmd += " --barcode_mtd_metadatacols %(barcode_mtd_metadatacols)s"
    cmd += " > logs/load_mudatas_%(sample_id)s.log"
    # print(cmd)
    P.run(cmd, job_threads=PARAMS['resources_threads_medium'], **job_kwargs)


@active_if(PARAMS["use_existing_h5ad"] is False)
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
    # if PARAMS['demultiplex_include'] is not False:
    #     cmd += " --demultiplexing %(demultiplex_include)s"
    #     cmd += " --demultiplexing_metadatacols %(demultiplex_metadatacols)s"
    # if PARAMS['use_muon']:
    #     cmd += " --use_muon True"
    cmd += " > logs/concat_filtered_mudatas.log"
    P.run(cmd, job_threads=PARAMS['resources_threads_high'], **job_kwargs)
    # P.run("rm tmp/*", job_threads=PARAMS['resources_threads_low'])


# -----------------------------------------------------------------------------------------------
## Creating h5mu from raw data files
# -----------------------------------------------------------------------------------------------

if PARAMS['assess_background'] or (PARAMS["modalities"]["prot"] and "dsb" in PARAMS['normalisation_methods']):
    PARAMS['raw_required'] = True
else:
    PARAMS['raw_required'] = False
    

def raw_file():
    sprefix = PARAMS['sample_prefix']
    if PARAMS['use_muon']:
        raw_file = sprefix + "_raw.h5mu"
    else:
        raw_file = sprefix + "_raw.h5ad"
    return raw_file


def gen_load_raw_anndata_jobs():
    caf = pd.read_csv(PARAMS['submission_file'], sep='\t')
    return gen_load_anndata_jobs(caf, load_raw=True, mode_dictionary=PARAMS["modalities"])
    

@active_if(PARAMS["raw_required"])
@follows(mkdir("logs"))
@follows(mkdir("tmp"))
@active_if(PARAMS["use_existing_h5ad"] is False)
@files(gen_load_raw_anndata_jobs)
def load_raw_mudatas(gex_path, outfile, 
              sample_id, 
              gex_filetype,  
              adt_path, adt_filetype, 
              tcr_path, tcr_filetype,  
              bcr_path, bcr_filetype, 
              atac_path, atac_filetype, 
              fragments_file, per_barcode_metrics_file, peak_annotation_file, 
              cell_mtd_path):
    # need to remove unnessacry modality
    mod_dict = PARAMS['modalities'].copy()
    mod_dict['tcr'] = False
    mod_dict['bcr'] = False
    cmd = """
    python %(py_path)s/make_adata_from_csv.py 
    --sample_id %(sample_id)s
    --mode_dictionary "%(mod_dict)s"
    --gex_infile %(gex_path)s
    --gex_filetype %(gex_filetype)s
    --output_file %(outfile)s 
    """
    if adt_path is not None and pd.notna(adt_path):
        cmd += " --adt_infile %(adt_path)s"
        cmd += " --adt_filetype %(adt_filetype)s"
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
    cmd += " > logs/load_raw_mudatas_%(sample_id)s.log"
    
    P.run(cmd, job_threads=PARAMS['resources_threads_medium'], **job_kwargs)

@active_if(PARAMS["raw_required"])
@active_if(PARAMS['downsample_raw'])
@follows(load_raw_mudatas)
@transform(load_raw_mudatas, regex("(.*)_raw.h5mu"), r"\1_raw_downsampled.log")
def downsample_raw_mudatas(infile, outfile):
    cmd = """
    python %(py_path)s/downsample.py 
    --input_anndata %(infile)s
    --output_anndata "%(infile)s" 
    --downsample_value 20000
    --use_muon %(use_muon)s  > %(outfile)s
    """
    P.run(cmd, job_threads=PARAMS['resources_threads_medium'], **job_kwargs)

@active_if(PARAMS["raw_required"])
@active_if(PARAMS["assess_background"])
@active_if(PARAMS["use_existing_h5ad"] is False)
@collate(load_raw_mudatas,
         formatter(""),
         raw_file())
def concat_raw_mudatas(infiles, outfile):
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
    # if PARAMS['demultiplex_include'] is not False:
    #     cmd += " --demultiplexing %(demultiplex_include)s"
    #     cmd += " --demultiplexing_metadatacols %(demultiplex_metadatacols)s"
    # if PARAMS['use_muon']:
    #     cmd += " --use_muon True"
    cmd += " > logs/concat_raw_mudatas.log"
    P.run(cmd, job_threads=PARAMS['resources_threads_high'], **job_kwargs)
    # P.run("rm tmp/*", job_threads=PARAMS['resources_threads_low'])

# -----------------------------------------------------------------------------------------------
## rna QC 
# -----------------------------------------------------------------------------------------------


def orfile():
    return PARAMS['sample_prefix'] + "_cell_metadata.tsv"
@active_if(PARAMS['modalities_rna'])
@active_if(PARAMS["use_existing_h5ad"] is False)
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
    P.run(cmd, job_threads=PARAMS['resources_threads_medium'],**job_kwargs)
    IOTools.touch_file(outfile)


@active_if(PARAMS['modalities_rna'])
@follows(mkdir("figures"))
@follows(mkdir("figures/gex"))
@follows(concat_filtered_mudatas, run_scrublet)
@originate("logs/run_scanpy_qc_rna.log", orfile(), unfilt_file())
# @transform(concat_mudata, suffix("_unfilt.h5ad"), "_cell_metadata.tsv")
def run_gex_qc(log_file, outfile, unfilt_file):
    # infile = submission file
    resources_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), "resources")
    customgenesfile = PARAMS["custom_genes_file"]
    calcproportions= PARAMS["calc_proportions"]
    cmd = """
        python %(py_path)s/run_scanpyQC_GEX.py
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
    P.run(cmd, job_threads=PARAMS['resources_threads_high'], **job_kwargs)
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
@follows(mkdir('figures/adt'))
@follows(concat_filtered_mudatas, run_gex_qc)
@originate("logs/run_scanpy_qc_prot.log", orfile(), unfilt_file())
def run_scanpy_adt_qc(log_file, outfile, unfilt_file):
    # infile = submission file
    cmd = """
        python %(py_path)s/run_scanpyQC_ADT.py
          --sampleprefix %(sample_prefix)s
          --input_anndata %(unfilt_file)s
          --outfile %(unfilt_file)s
          --figdir figures/adt/
          """
    if PARAMS['channel_col'] is not None:
        cmd += " --channel_col %(channel_col)s"
    else:
        cmd += " --channel_col sample_id"
    if PARAMS['adt_metrics_per_cell']:
        cmd += " --per_cell_metrics %(adt_metrics_per_cell)s"
    if PARAMS['adt_metrics_per_adt']:
        cmd += " --per_adt_metrics %(adt_metrics_per_adt)s"
    if PARAMS['identify_isotype_outliers']:
        cmd += " --identify_isotype_outliers %(identify_isotype_outliers)s"
        cmd += " --isotype_upper_quantile %(isotype_upper_quantile)s"
        cmd += " --isotype_n_pass %(isotype_n_pass)s"
    # add log file
    cmd += " > %(log_file)s"
    P.run(cmd, job_threads=PARAMS['resources_threads_high'], **job_kwargs)
    pass

@active_if(PARAMS['modalities_prot'])
@follows(run_scanpy_adt_qc, concat_filtered_mudatas, concat_raw_mudatas)
@originate("logs/run_dsb_clr.log", unfilt_file(), raw_file())
def run_dsb_clr(outfile, unfilt_file, raw_file):
    print(unfilt_file)
    # infile = mdata_unfilt
    # outfile sampleprefix + "adt_qc_per_cell.tsv"
    cmd = """
        python %(py_path)s/run_preprocess_prot.py
        --filtered_mudata %(unfilt_file)s
        --channel_col %(channel_col)s
        --figpath ./figures/adt
        """
    if PARAMS['normalisation_methods'] is not None:
        cmd += " --normalisation_methods %(normalisation_methods)s"
    if PARAMS['quantile_clipping'] is not None:
        cmd += " --quantile_clipping %(quantile_clipping)s"
    if PARAMS['clr_margin'] is not None:
        cmd += " --clr_margin %(clr_margin)s"
    # find the raw data if available
    if os.path.exists(raw_file):
        cmd += " --raw_mudata %(raw_file)s"
    cmd += " > %(outfile)s"
    P.run(cmd, job_threads=PARAMS['resources_threads_low'], **job_kwargs)

@follows(run_scanpy_adt_qc, run_dsb_clr)
def run_adt_qc():
    pass


# -----------------------------------------------------------------------------------------------
## Repertoire QC  # Repertoire run_scanpy_qc_rep.log
# -----------------------------------------------------------------------------------------------
# toggle_rep_qc = ()
@active_if(PARAMS['modalities_bcr'] or PARAMS['modalities_tcr']  )
@follows(mkdir('figures/rep'))
@follows(run_gex_qc, run_adt_qc)
@originate("logs/run_scanpy_qc_rep.log", unfilt_file())
def run_repertoire_qc(logfile, unfilt_file):
    cmd = """python %(py_path)s/run_scanpyQC_REP.py
          --sampleprefix %(sample_prefix)s
          --input_mudata %(unfilt_file)s
          --output_mudata %(unfilt_file)s
          --figdir figures/rep/
          --distance_metrics "%(ir_dist)s"
          --clonotype_metrics "%(clonotype_definition)s"
          """
    cmd += " > %(logfile)s"
    P.run(cmd, job_threads=PARAMS['resources_threads_low'], **job_kwargs)

# -----------------------------------------------------------------------------------------------
## atac QC # run_scanpy_qc_atac.log
# -----------------------------------------------------------------------------------------------
@active_if(PARAMS['modalities_atac'])
@mkdir('figures/atac')
@follows(run_gex_qc, run_adt_qc, run_repertoire_qc)
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
        python %(py_path)s/run_scanpyQC_ATAC.py
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
    if PARAMS['use_muon']:
        cmd += " --use_muon True"

    cmd += " > %(log_file)s"
    P.run(cmd, job_threads=PARAMS['resources_threads_low'], **job_kwargs)

@follows(run_gex_qc, run_adt_qc, run_repertoire_qc, run_atac_qc)
def run_qc():
    pass


# -----------------------------------------------------------------------------------------------
## QC plotting
# -----------------------------------------------------------------------------------------------

@follows(run_qc)
# @transform(run_gex_qc,
#             regex(r"(.*)_cell_metadata.tsv"),
#             r"logs/plot_qc.log", )
@originate("logs/plot_qc.log", orfile())
def plot_qc(log_file, cell_file):
    R_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), "R")
    qcmetrics = PARAMS['plotqc_gex_metrics']
    cmd = """
    Rscript %(R_path)s/plotQC.R 
    --prefilter TRUE
    --cell_metadata %(cell_file)s 
    --sampleprefix %(sample_prefix)s
    --qc_metrics %(qcmetrics)s
    --groupingvar %(plotqc_grouping_var)s
    """
    if PARAMS['use_muon']:
        cmd += " --scanpy_or_muon muon"
    else:
        cmd += " --scanpy_or_muon scanpy"
    cmd += " > %(log_file)s"
    P.run(cmd, job_threads=PARAMS['resources_threads_low'], **job_kwargs)


# ------
# assess background
# -------

@active_if(PARAMS['assess_background'])
@follows(mkdir("figures/background"))
@follows(run_qc)
@follows(concat_raw_mudatas)
@originate("logs/assess_background.log", unfilt_file(), raw_file())
def run_assess_background(log_file, unfilt_file, raw_file):
    cmd = """
    python %(py_path)s/assess_background.py
        --filtered_mudata %(unfilt_file)s
        --raw_mudata %(raw_file)s
        --channel_col %(channel_col)s
        --figpath "./figures/background/"
    """
    cmd += " > %(log_file)s"
    P.run(cmd, job_threads=PARAMS['resources_threads_medium'], **job_kwargs)



# # ------------
# TO DO change the decorator to follow all qc plots?
@follows(plot_qc)
@follows(run_gex_qc, plot_tenx_metrics, run_assess_background)
def all_gex_qc():
    pass

@follows(run_dsb_clr, plot_qc)
def all_adt_qc():
    pass


@follows(all_adt_qc)
@follows(all_gex_qc)
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
