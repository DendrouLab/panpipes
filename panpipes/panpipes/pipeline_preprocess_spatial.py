from ruffus import *
import sys
import os
from cgatcore import pipeline as P
import logging
from panpipes.funcs.io import dictionary_stripper



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
        python %(py_path)s/run_filter_spatialT.py
        --input_mudata %(unfiltered_obj)s
        --output_mudata %(outfile)s
        --filter_dict "%(filter_dict)s"
        """
        if PARAMS['filtering_keep_barcodes'] is not None:
            cmd += " --keep_barcodes %(filtering_keep_barcodes)s"
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

@active_if(mode_dictionary['spatialT'] is True)
@active_if(PARAMS['filtering_run'])
@active_if(run_plotqc_query(PARAMS['plotqc']))
@follows(filter_mudata)
@originate("logs/postfilterplot_spatialT.log" , PARAMS['mudata_file'])
def postfilterplot_spatialT(log_file, filt_file):
    cmd = """
            python %(py_path)s/plot_qc_spatialT.py
             --input_mudata %(filt_file)s
             --output_mudata %(filt_file)s
             --figdir ./figures/spatialT
            """

    if PARAMS['plotqc']['grouping_var'] is not None:
        cmd += " --grouping_var %(plotqc_grouping_var)s"
    if PARAMS['plotqc']['spatialT_metrics'] is not None:
        cmd += " --spatialT_qc_metrics %(plotqc_spatialT_metrics)s"
    cmd += " > %(log_file)s "
    job_kwargs["job_threads"] = PARAMS['resources_threads_low']
    P.run(cmd, **job_kwargs)





@active_if(mode_dictionary['spatialT'] is True)
@follows(filter_mudata)
@originate("logs/preprocess_spatialT.log", PARAMS['mudata_file'])
def spatialT_preprocess(log_file, filt_file):
    if os.path.exists("figures/spatialT") is False:
        os.mkdir("figures/spatialT")
    cmd = """
        python %(py_path)s/run_preprocess_spatialT.py
        --input_mudata %(filt_file)s
        --output_mudata %(filt_file)s
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




@follows(filter_mudata, postfilterplot_spatialT, spatialT_preprocess)
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
