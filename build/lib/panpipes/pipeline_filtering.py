from ruffus import *
import sys
import os
from cgatcore import pipeline as P
import pandas as pd
import cgatcore.iotools as IOTools
import re
from itertools import chain
import glob
# from itertools import chain
# import glob

# import pandas as pd

PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "../pipeline.yml",
     "pipeline.yml"])


PARAMS['py_path'] = os.path.join(os.path.dirname(os.path.dirname(__file__)), "python")
PARAMS['R_path'] = os.path.join(os.path.dirname(os.path.dirname(__file__)), "R")
PARAMS['resources_path'] = os.path.join(os.path.dirname(os.path.dirname(__file__)), "resources")

PARAMS['filt_file'] = PARAMS['sample_prefix'] + "_filt.h5mu"
PARAMS['scaled_file'] = PARAMS['sample_prefix'] + "_scaled.h5mu"
PARAMS['mudata_file'] = PARAMS['sample_prefix'] + ".h5mu"

mode_dictionary = PARAMS["modalities"]

@mkdir("logs")
@mkdir("figures")
@originate(PARAMS['mudata_file'])
def filter_mudata(outfile):
    if PARAMS['filtering_run']:
        cmd = """
        python %(py_path)s/run_filter.py
        --input_mudata %(unfiltered_obj)s
        --output_mudata %(outfile)s
        --filter_dict "%(filtering)s"
        """
        if PARAMS['filtering_keep_barcodes'] is not None:
            cmd += " --keep_barcodes %(filtering_keep_barcodes)s"
        if PARAMS['intersect_mods'] is not None:
            cmd += " --intersect_mods %(intersect_mods)s"
        cmd += " > logs/filtering.log "
        P.run(cmd,  job_threads=PARAMS['resources_threads_low'])
    else:
        try:
            f = open(outfile)
            f.close(outfile)
        except IOError:
            print("filter is not set to True, but there is no %s file" % outfile)
            sys.exit(1)


@active_if(PARAMS['plotqc_metrics'] is not None)
@follows(filter_mudata)
@originate("logs/postfilterplot.log")
def postfilterplot(log_file):
    R_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), "R")
    cell_mtd_file = PARAMS['sample_prefix'] + "_filtered_cell_metadata.tsv"
    cmd = """
    Rscript %(R_path)s/plotQC.R 
    --prefilter FALSE
    --cell_metadata %(cell_mtd_file)s 
    --sampleprefix %(sample_prefix)s
    --groupingvar %(plotqc_grouping_var)s
    """
    cmd += " --scanpy_or_muon muon"
    cmd += " > %(log_file)s "
    P.run(cmd, job_threads=PARAMS['resources_threads_low'])




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
         --batch_col %(downsample_col)s
         --intersect_mods %(downsample_mods)s
        """
        cmd += " > %(log_file)s  "
        P.run(cmd, job_threads=PARAMS['resources_threads_low'])
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
        if PARAMS['hvg_exclude'] == "default":
            hvg_exclude = PARAMS['resources_path'] + "/exclude_genes_HLAIGTR_v1.txt"
        cmd += " --exclude_file %(hvg_exclude)s"
    if PARAMS['hvg_flavor'] is not None:
        cmd += " --flavor %(hvg_flavor)s"
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
    if PARAMS['scale_max_value'] is not None:
        cmd += " --scale_max_value %(scale_max_values)s"
    if PARAMS['pca_scree_n_pcs'] is not None:
        cmd += " --n_pcs %(pca_scree_n_pcs)s"
    if PARAMS['pca_color_by'] is not None:
        cmd += " --color_by %(pca_color_by)s"
    if PARAMS['output_logged_mudata'] is not None:
        cmd += " --output_logged_mudata %(unscaled_outfile)s"
    cmd += " > %(log_file)s "
    P.run(cmd, job_threads=PARAMS['resources_threads_high'])


@active_if(mode_dictionary['prot'] is True)
# @transform(rna_preprocess,formatter(),"logs/preprocess_prot.log")
@follows(rna_preprocess)
@originate("logs/preprocess_prot.log", PARAMS['mudata_file'])
def prot_preprocess( log_file, scaled_file, ):
    if os.path.exists("figures/adt") is False:
        os.mkdir("figures/adt" )
    cmd = """
        python %(py_path)s/run_preprocess_prot.py
        --filtered_mudata %(scaled_file)s
        --figpath ./figures/adt
        --save_mudata_path %(scaled_file)s
        """
    if PARAMS['normalisation_methods'] is not None:
        cmd += " --normalisation_methods %(adt_normalisation_methods)s"
        cmd += " --clr_margin %(adt_clr_margin)s"
        cmd += " --quantile_clipping %(adt_quantile_clipping)s"
    if PARAMS['raw_obj'] is not None:
        cmd += " --raw_mudata %(raw_obj)s"
    cmd += " > %(log_file)s"
    P.run(cmd, job_threads=PARAMS['resources_threads_high'])


@active_if(mode_dictionary['atac'] is True)
# @follows(prot_preprocess)
# @transform(rna_preprocess,formatter(),"logs/preprocess_atac.log")
@originate("logs/preprocess_atac.log", add_inputs(PARAMS['mudata_file']))
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
        cmd +- " --normalize %(atac_normalize)s"
    if PARAMS['atac_TFIDF_flavour'] is not None:
        cmd +- " --TFIDF_flavour %(atac_TFIDF_flavour)s"
    if PARAMS['atac_min_mean'] is not None:
        cmd +- " --min_mean %(atac_min_mean)s"    
    if PARAMS['atac_max_mean'] is not None:
        cmd +- " --max_mean %(atac_max_mean)s"    
    if PARAMS['atac_min_disp'] is not None:
        cmd +- " --min_disp %(atac_min_disp)s"    
    cmd += " > %(log_file)s"
    P.run(cmd, job_threads=PARAMS['resources_threads_high'])

@active_if(mode_dictionary['rep'] is True)
@follows(atac_preprocess)
# @transform(rna_preprocess,formatter(),"logs/preprocess_rep.log")
@originate("logs/preprocess_rep.log", add_inputs(PARAMS['mudata_file']))
def rep_preprocess( log_file, scaled_file):
    pass

# ---- end stub
@follows(postfilterplot,rep_preprocess)
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
