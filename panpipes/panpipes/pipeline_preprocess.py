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
from panpipes.funcs.io import dictionary_stripper
# from itertools import chain
# import glob

# import pandas as pd

PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "pipeline.yml"])


PARAMS['py_path'] =  os.path.join(os.path.dirname(os.path.dirname(__file__)), 'python_scripts')
PARAMS['r_path'] = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'R_scripts')
PARAMS['resources_path'] = os.path.join(os.path.dirname(os.path.dirname(__file__)), "resources")

if  PARAMS['sample_prefix'] is not None:
    PARAMS['mudata_file'] = PARAMS['sample_prefix'] + ".h5mu"
else:
    logging.warning("sample prefix is None, please use a string")
    
job_kwargs = {}

if PARAMS['condaenv'] is not None:
    job_kwargs["job_condaenv"] =PARAMS['condaenv']


mode_dictionary = PARAMS["modalities"]
@mkdir("logs")
@mkdir("figures")
@originate(PARAMS['mudata_file'])
def filter_mudata(outfile):
    if PARAMS['filtering_run']:
        filter_dict = dictionary_stripper(PARAMS['filtering'])
        cmd = """
        python %(py_path)s/run_filter.py
        --input_mudata %(unfiltered_obj)s
        --output_mudata %(outfile)s
        --filter_dict "%(filter_dict)s"
        """
        if PARAMS['filtering_keep_barcodes'] is not None:
            cmd += " --keep_barcodes %(filtering_keep_barcodes)s"
        if PARAMS['intersect_mods'] is not None:
            cmd += " --intersect_mods %(intersect_mods)s"
        cmd += " > logs/filtering.log "
        job_kwargs["job_threads"] = PARAMS['resources_threads_low']
        P.run(cmd, **job_kwargs)
    else:
        try:
            f = open(outfile)
            f.close(outfile)
        except IOError:
            print("filter is not set to True, but there is no %s file" % outfile)
            sys.exit(1)


def run_plotqc_query(pqc_dict):
    # avoid deleting from PARAMS
    pqc = pqc_dict.copy()
    del pqc['grouping_var']
    return any([x != None for x in pqc.values()])


@active_if(PARAMS['filtering_run'])
@active_if(run_plotqc_query(PARAMS['plotqc']))
@follows(filter_mudata)
@originate("logs/postfilterplot.log")
def postfilterplot(log_file):
    cell_mtd_file = PARAMS['sample_prefix'] + "_filtered_cell_metadata.tsv"
    cmd = """
    Rscript %(r_path)s/plotQC.R 
    --prefilter FALSE
    --cell_metadata %(cell_mtd_file)s 
    --sampleprefix %(sample_prefix)s
    --groupingvar %(plotqc_grouping_var)s
    """
    cmd += " --scanpy_or_muon muon"
    if PARAMS['plotqc']['rna_metrics'] is not None:
        cmd += " --rna_qc_metrics %(plotqc_rna_metrics)s"
    if PARAMS['plotqc']['prot_metrics'] is not None:
        cmd += " --prot_qc_metrics %(plotqc_prot_metrics)s"
    if PARAMS['plotqc']['atac_metrics'] is not None:
        cmd += " --atac_qc_metrics %(plotqc_atac_metrics)s"
    if PARAMS['plotqc']['rep_metrics'] is not None:
        cmd += " --rep_qc_metrics %(plotqc_rep_metrics)s"
    cmd += " > %(log_file)s "
    job_kwargs["job_threads"] = PARAMS['resources_threads_low']
    P.run(cmd, **job_kwargs)




@active_if(PARAMS['downsample_n'] is not None)
@follows(filter_mudata)
@originate("logs/downsample.log", PARAMS['mudata_file'])
def downsample(log_file, filt_obj):
    if PARAMS['downsample'] is not None:
        cmd="""
        python %(py_path)s/downsample.py
         --input_mudata %(filt_obj)s \
         --output_mudata %(filt_obj)s \
         --sampleprefix %(sample_prefix)s \
         --downsample_value %(downsample_n)s \
         --downsample_col %(downsample_col)s
         --intersect_mods %(downsample_mods)s
        """
        cmd += " > %(log_file)s  "
        job_kwargs["job_threads"] = PARAMS['resources_threads_low']
        P.run(cmd, **job_kwargs)
    IOTools.touch_file(log_file)
    pass


# setting these follows means that it still works if filter is False
@active_if(mode_dictionary['rna'] is True)
@follows(downsample)
@transform(filter_mudata, formatter(), "logs/preprocess_rna.log")
def rna_preprocess(adata_obj, log_file):
    cmd = """python %(py_path)s/run_preprocess_rna.py 
            --input_mudata %(adata_obj)s 
            --fig_dir figures/ 
            --output_scaled_mudata %(adata_obj)s
            """
    # add in the options if specified in pipeline.yml
    if PARAMS['hvg_exclude'] is not None:
        cmd += " --exclude_file %(hvg_exclude_file)s"
        if PARAMS['hvg_exclude'] =="default":
            cmd += " --exclude exclude"
        else:    
            cmd += " --exclude %(hvg_exclude)s"
    if PARAMS['hvg_flavor'] is not None:
        cmd += " --flavor %(hvg_flavor)s"
    if PARAMS['hvg_batch_key'] is not None:
        cmd += " --hvg_batch_key %(hvg_batch_key)s"
    if PARAMS['hvg_n_top_genes'] is not None:
        cmd += " --n_top_genes %(hvg_n_top_genes)s"
    if PARAMS['hvg_min_mean'] is not None:
        cmd += " --min_mean %(hvg_min_mean)s"
    if PARAMS['hvg_max_mean'] is not None:
        cmd += " --max_mean %(hvg_max_mean)s"
    if PARAMS['hvg_min_disp'] is not None:
        cmd += " --min_disp %(hvg_min_disp)s"
    if PARAMS['hvg_filter'] is not None:
        cmd += " --filter_by_hvg %(hvg_filter)s"
    if PARAMS['regress_variables'] is not None:
        cmd += " --regress_out %(regress_variables)s"
    if PARAMS['run_scale'] is True:
        cmd += " --scale True"
        if PARAMS['scale_max_value'] is not None:
            cmd += " --scale_max_value %(scale_max_value)s"
    else:
        cmd += " --scale False"
    if PARAMS['pca_scree_n_pcs'] is not None:
        cmd += " --n_pcs %(pca_scree_n_pcs)s"
    if PARAMS['pca_color_by'] is not None:
        cmd += " --color_by %(pca_color_by)s"
    if PARAMS['output_logged_mudata'] is not None:
        cmd += " --output_logged_mudata %(unscaled_outfile)s"
    cmd += " > %(log_file)s "
    job_kwargs["job_threads"] = PARAMS['resources_threads_high']
    P.run(cmd, **job_kwargs)


@active_if(mode_dictionary['prot'] is True)
# @transform(rna_preprocess,formatter(),"logs/preprocess_prot.log")
@follows(rna_preprocess)
@originate("logs/preprocess_prot.log", PARAMS['mudata_file'])
def prot_preprocess( log_file, scaled_file, ):
    if os.path.exists("figures/prot") is False:
        os.mkdir("figures/prot" )
    cmd = """
        python %(py_path)s/run_preprocess_prot.py
        --filtered_mudata %(scaled_file)s
        --figpath ./figures/prot
        --save_mudata_path %(scaled_file)s
        """
    if PARAMS['normalisation_methods'] is not None:
        cmd += " --normalisation_methods %(prot_normalisation_methods)s"
    if PARAMS['quantile_clipping'] is not None:
        cmd += " --quantile_clipping %(prot_quantile_clipping)s"
    if PARAMS['clr_margin'] is not None:
        cmd += " --clr_margin %(prot_clr_margin)s"
    if PARAMS['background_obj'] is not None:
        cmd += " --bg_mudata %(prot_background_obj)s"
    if PARAMS['prot']['store_as_X'] is not None:
        cmd += " --store_as_x %(prot_store_as_X)s"
    if PARAMS['prot_save_norm_prot_mtx'] is True:
        cmd += " --save_mtx True"
    cmd += " > %(log_file)s"
    job_kwargs["job_threads"] = PARAMS['resources_threads_high']
    P.run(cmd, **job_kwargs)


@active_if(mode_dictionary['atac'] is True)
@follows(prot_preprocess)
# @transform(rna_preprocess,formatter(),"logs/preprocess_atac.log")
@originate("logs/preprocess_atac.log", PARAMS['mudata_file'])
def atac_preprocess(log_file, scaled_file):
    if os.path.exists("figures/atac") is False:
        os.mkdir("figures/atac" )
    cmd = """
        python %(py_path)s/run_preprocess_atac.py
        --input_mudata %(scaled_file)s
        --output_mudata %(scaled_file)s
        --figdir ./figures/atac
        """
    if PARAMS['atac_binarize'] is True:
        cmd += " --binarize True" 
    else:
        cmd += " --binarize False" 
    if PARAMS['atac_normalize'] is not None:
        cmd += " --normalize %(atac_normalize)s"
    if PARAMS['atac_TFIDF_flavour'] is not None:
        cmd += " --TFIDF_flavour %(atac_TFIDF_flavour)s"
    if PARAMS['atac_min_mean'] is not None:
        cmd += " --min_mean %(atac_min_mean)s"    
    if PARAMS['atac_max_mean'] is not None:
        cmd += " --max_mean %(atac_max_mean)s"    
    if PARAMS['atac_min_disp'] is not None:
        cmd += " --min_disp %(atac_min_disp)s"    
    if PARAMS['atac_dimred'] is not None:
        cmd += " --dimred %(atac_dimred)s"
    else:    
        cmd += " --dimred PCA"
    if PARAMS['atac_dim_remove'] is not None:
        cmd += " --dim_remove %(atac_dim_remove)s"

    if PARAMS['atac_feature_selection_flavour'] is not None:
        cmd += " --feature_selection_flavour %(atac_feature_selection_flavour)s"
    if PARAMS['atac_min_cutoff'] is not None:
        cmd += " --min_cutoff %(atac_min_cutoff)s"
    
    cmd += " > %(log_file)s"
    job_kwargs["job_threads"] = PARAMS['resources_threads_high']
    P.run(cmd, **job_kwargs)


'''
@active_if(mode_dictionary['spatialT'] is True)
@follows(atac_preprocess)
@originate("logs/preprocess_spatialT.log", PARAMS['mudata_file'])
def spatialT_preprocess(log_file, scaled_file):
    if os.path.exists("figures/spatialT") is False:
        os.mkdir("figures/spatialT")
    cmd = """
        python %(py_path)s/run_preprocess_spatialT.py
        --input_mudata %(scaled_file)s
        --output_mudata %(scaled_file)s
        --figdir ./figures/spatialT
        """
    if PARAMS['spatialT_norm_hvg_flavour'] is not None:
        cmd += " --norm_hvg_flavour %(spatialT_norm_hvg_flavour)s"
    if PARAMS['spatialT_n_top_genes'] is not None:
        cmd += " --n_top_genes %(spatialT_n_top_genes)s"
    if PARAMS['spatialT_filter_by_hvg'] is True:
        cmd += " --filter_by_hvg True"
    else:
        cmd += " --filter_by_hvg False"
    if PARAMS['spatialT_hvg_batch_key'] is not None:
        cmd += " --hvg_batch_key %(spatialT_hvg_batch_key)s"
    if PARAMS['spatialT_squidpy_hvg_flavour'] is not None:
        cmd += " --squidpy_hvg_flavour %(spatialT_squidpy_hvg_flavour)s"
    if PARAMS['spatialT_min_mean'] is not None:
        cmd += " --min_mean %(spatialT_min_mean)s"
    if PARAMS['spatialT_max_mean'] is not None:
        cmd += " --max_mean %(spatialT_max_mean)s"
    if PARAMS['spatialT_min_disp'] is not None:
        cmd += " --min_disp %(spatialT_min_disp)s"
    if PARAMS['spatialT_theta'] is not None:
        cmd += " --theta %(spatialT_theta)s"
    if PARAMS['spatialT_clip'] is not None:
        cmd += " --clip %(spatialT_clip)s"
    if PARAMS['spatialT_n_pcs'] is not None:
        cmd += " --n_pcs %(spatialT_n_pcs)s"

    cmd += " > %(log_file)s"
    job_kwargs["job_threads"] = PARAMS['resources_threads_high']
    P.run(cmd, **job_kwargs)
'''



# @active_if(mode_dictionary['rep'] is True)
# @follows(atac_preprocess)
# # @transform(rna_preprocess,formatter(),"logs/preprocess_rep.log")
# @originate("logs/preprocess_rep.log", PARAMS['mudata_file'])
# def rep_preprocess( log_file, scaled_file):
#     pass

# ---- end stub
@follows(postfilterplot,rna_preprocess,atac_preprocess,prot_preprocess)
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
