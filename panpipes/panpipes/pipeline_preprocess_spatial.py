from ruffus import *
import sys
import os
from cgatcore import pipeline as P
import glob
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

if PARAMS["input_dir"] is not None:
    input_dir= PARAMS["input_dir"]
    if not os.path.exists(input_dir):
        logging.warning("input directory doesn't exist, checking if files are in the default ../qc.data/ dir")
        input_dir = "../qc.data"
        if not os.path.exists(input_dir):
            sys.exit("can't find input data")
else:
    input_dir = "../qc.data"
    if not os.path.exists(input_dir):
            sys.exit("can't find input data")



def gen_filter_jobs():
    input_paths=glob.glob(os.path.join(input_dir,"*unfilt.h5mu"))
    for infile_path in input_paths:
        file_name = os.path.basename(infile_path)
        outfile = file_name.replace("unfilt","filtered")
        yield infile_path, outfile
    

@mkdir("logs")
@mkdir("figures")
@mkdir("tables")
@mkdir("filtered.data")
@files(gen_filter_jobs)
def filter_mudata(infile_path,outfile):
    print('processing file = %s' % str(infile_path))
    log_file = os.path.basename(outfile)
    log_file= "filtering."+log_file.replace("filtered.h5mu","") + ".log"
    if PARAMS['filtering_run']:
        filter_dict = dictionary_stripper(PARAMS['filtering'])
        cmd = """
        python %(py_path)s/run_filter_spatial.py
        --input_mudata %(infile_path)s
        --output_mudata filtered.data/%(outfile)s
        --filter_dict "%(filter_dict)s"
        """
        if PARAMS['filtering_keep_barcodes'] is not None:
            cmd += " --keep_barcodes %(filtering_keep_barcodes)s"
        cmd += " > logs/%(log_file)s "
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

@active_if(mode_dictionary['spatial'] is True)
@active_if(PARAMS['filtering_run'])
@active_if(run_plotqc_query(PARAMS['plotqc']))
@follows(filter_mudata)
@originate("logs/postfilterplot_spatial.log" , PARAMS['mudata_file'])
def postfilterplot_spatial(log_file, filt_file):
    cmd = """
            python %(py_path)s/plot_qc_spatial.py
             --input_mudata %(filt_file)s
             --output_mudata %(filt_file)s
             --figdir ./figures/spatial
            """

    if PARAMS['plotqc']['grouping_var'] is not None:
        cmd += " --grouping_var %(plotqc_grouping_var)s"
    if PARAMS['plotqc']['spatial_metrics'] is not None:
        cmd += " --spatial_qc_metrics %(plotqc_spatial_metrics)s"
    cmd += " > %(log_file)s "
    job_kwargs["job_threads"] = PARAMS['resources_threads_low']
    P.run(cmd, **job_kwargs)





@active_if(mode_dictionary['spatial'] is True)
@follows(filter_mudata)
@originate("logs/preprocess_spatial.log", PARAMS['mudata_file'])
def spatial_preprocess(log_file, filt_file):
    if os.path.exists("figures/spatial") is False:
        os.mkdir("figures/spatial")
    cmd = """
        python %(py_path)s/run_preprocess_spatial.py
        --input_mudata %(filt_file)s
        --output_mudata %(filt_file)s
        --figdir ./figures/spatial
        """
    if PARAMS['spatial_norm_hvg_flavour'] is not None:
        cmd += " --norm_hvg_flavour %(spatial_norm_hvg_flavour)s"
    if PARAMS['spatial_n_top_genes'] is not None:
        cmd += " --n_top_genes %(spatial_n_top_genes)s"
    if PARAMS['spatial_filter_by_hvg'] is True:
        cmd += " --filter_by_hvg True"
    else:
        cmd += " --filter_by_hvg False"
    if PARAMS['spatial_hvg_batch_key'] is not None:
        cmd += " --hvg_batch_key %(spatial_hvg_batch_key)s"
    if PARAMS['spatial_squidpy_hvg_flavour'] is not None:
        cmd += " --squidpy_hvg_flavour %(spatial_squidpy_hvg_flavour)s"
    if PARAMS['spatial_min_mean'] is not None:
        cmd += " --min_mean %(spatial_min_mean)s"
    if PARAMS['spatial_max_mean'] is not None:
        cmd += " --max_mean %(spatial_max_mean)s"
    if PARAMS['spatial_min_disp'] is not None:
        cmd += " --min_disp %(spatial_min_disp)s"
    if PARAMS['spatial_theta'] is not None:
        cmd += " --theta %(spatial_theta)s"
    if PARAMS['spatial_clip'] is not None:
        cmd += " --clip %(spatial_clip)s"
    if PARAMS['spatial_n_pcs'] is not None:
        cmd += " --n_pcs %(spatial_n_pcs)s"

    cmd += " > %(log_file)s"
    job_kwargs["job_threads"] = PARAMS['resources_threads_high']
    P.run(cmd, **job_kwargs)




@follows(filter_mudata, postfilterplot_spatial, spatial_preprocess)
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
