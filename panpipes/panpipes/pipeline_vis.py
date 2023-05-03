
from ruffus import *
import sys
import os
from cgatcore import pipeline as P
import pandas as pd
import cgatcore.iotools as IOTools
import re
from itertools import chain, product
import glob
from panpipes.funcs.io import dictionary_stripper 

PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "pipeline.yml"])

PARAMS['py_path'] =  os.path.join(os.path.dirname(os.path.dirname(__file__)), 'python_scripts')
PARAMS['r_path'] = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'R_scripts')

job_kwargs = {}

if PARAMS['condaenv'] is not None:
    job_kwargs["job_condaenv"] =PARAMS['condaenv']


@originate("logs/setup_dirs.sentinel")
def set_up_dirs(log_file):
    if os.path.exists("logs") is False:
        os.mkdir("logs")
    mods =  [key for key, value in PARAMS['modalities'].items() if value is True]
    for mod in mods:
        if os.path.exists(mod) is False:
            os.makedirs(os.path.join(mod))
    IOTools.touch_file(log_file)
    pass

#--------------------------------------------------
# Plot markers
#--------------------------------------------------



marker_files = PARAMS['custom_markers']['files']
marker_files  = list(chain(*[marker_files[x]  for x in ['full', 'minimal'] if marker_files[x]  is not None]))


do_plot_features = (PARAMS['custom_markers_files'] is not None or PARAMS['custom_markers_minimal'] is not None) and \
                (PARAMS['do_plots_marker_dotplots'] or PARAMS['do_plots_marker_matrixplots'])


@follows(set_up_dirs)
@active_if(do_plot_features)
@transform(marker_files, formatter(), "logs/markers_by_group__{basename[0]}.log")
def plot_custom_markers_per_group(marker_file, log_file):
    print(marker_file)
    print(log_file)
    group_vars = " ".join(PARAMS["grouping_vars"])
    # pull out the modalities to use, not including multimodal.
    modalities = [key for key, val in PARAMS['modalities'].items() if val and key not in ["multimodal", "rep"]]
    modalities = ",".join(modalities)
    # pull out the (optional) layer choices
    layer_vars = {key:val  for key, val in PARAMS['custom_markers_layers'].items() if PARAMS['modalities'][key]}
    layer_vars = dictionary_stripper(layer_vars)
    cmd = """
        python %(py_path)s/plot_custom_markers.py \
            --infile %(mudata_obj)s \
            --modalities %(modalities)s \
            --layers "%(layer_vars)s" \
            --marker_file %(marker_file)s \
            --group_cols %(group_vars)s \
            --base_figure_dir ./
     """
    cmd += " > %(log_file)s "
    job_kwargs["job_threads"] = PARAMS['resources_threads_high']
    P.run(cmd, **job_kwargs)


plot_embeddings = any([True if val['run'] is True else False for val in PARAMS["embedding"].values() ])
@active_if(plot_embeddings )
@follows(set_up_dirs)
@active_if(PARAMS['custom_markers_minimal'] is not None)
@transform(PARAMS['custom_markers']['files']['minimal'], formatter(), "logs/markers_umap_{basename[0]}.log")
def plot_custom_markers_umap(marker_file, log_file):
    embedding_dict = PARAMS["embedding"]
    embedding_dict =  {mod:val['basis'] for mod, val in embedding_dict.items() if val['run'] is True}
    embedding_dict = dictionary_stripper(embedding_dict)
    layer_vars = {key:val  for key, val in PARAMS['custom_markers_layers'].items() if PARAMS['modalities'][key]}
    layer_vars = dictionary_stripper(layer_vars)
    modalities = ",".join([key for key, val in PARAMS['modalities'].items() if val])
    cmd = """
        python %(py_path)s/plot_custom_markers_umap.py \
            --infile %(mudata_obj)s \
            --modalities %(modalities)s \
            --layers "%(layer_vars)s" \
            --basis_dict "%(embedding_dict)s" \
            --marker_file %(marker_file)s \
            --base_figure_dir ./
    """
    cmd += " > %(log_file)s "
    job_kwargs["job_threads"] = PARAMS['resources_threads_high']
    P.run(cmd, **job_kwargs)




# @transform(PARAMS['custom_markers']['files']['minimal'], formatter(), "logs/markers_umap__{basename[0]}.log")
@active_if(plot_embeddings )
@follows(set_up_dirs)
@originate("logs/categorical_variables_umap.log")
def plot_categorical_umaps(log_file):
    embedding_dict = PARAMS["embedding"]
    embedding_dict =  {mod:val['basis'] for mod, val in embedding_dict.items() if val['run'] is True}
    cat_vars = PARAMS['categorical_vars']
    if "all" in cat_vars.keys() and cat_vars["all"] is not None:    
        for k, v in cat_vars.items():
            if v is None:
                cat_vars[k] = cat_vars['all']
            else:
                cat_vars[k] = cat_vars[k] + cat_vars['all']
    del cat_vars['all']
    # subset by the modalities we are using according to PARAMS['modalities']
    cat_vars = {key:val  for key, val in cat_vars.items() if PARAMS['modalities'][key] }
    layer_vars = {key:val  for key, val in PARAMS['custom_markers_layers'].items() if PARAMS['modalities'][key]}
    cmd = """
        python %(py_path)s/plot_variables_umaps.py \
            --infile %(mudata_obj)s \
            --layers "%(layer_vars)s" \
            --basis_dict "%(embedding_dict)s" \
            --categorical_variables "%(cat_vars)s" \
            --base_figure_dir ./
            --fig_suffix categorical_vars.png
    """
    cmd += " > %(log_file)s "
    job_kwargs["job_threads"] = PARAMS['resources_threads_high']
    P.run(cmd, **job_kwargs)


@active_if(plot_embeddings)
@follows(set_up_dirs)
@originate("logs/continuous_variables_umap.log")
def plot_continuous_umaps(log_file):
    embedding_dict = PARAMS["embedding"]
    embedding_dict =  {mod:val['basis'] for mod, val in embedding_dict.items() if val['run'] is True}
    cont_vars = PARAMS['continuous_vars']
    if "all" in cont_vars.keys() and cont_vars["all"] is not None:    
        for k, v in cont_vars.items():
            if v is None:
                cont_vars[k] = cont_vars['all']
            else:
                cont_vars[k] = cont_vars[k] + cont_vars['all']
    del cont_vars['all']
    # subset by the modalities we are using according to PARAMS['modalities']
    cont_vars = {key:val  for key, val in cont_vars.items() if PARAMS['modalities'][key] }
    layer_vars = {key:val  for key, val in PARAMS['custom_markers_layers'].items() if PARAMS['modalities'][key]}
  
    cmd = """
        python %(py_path)s/plot_variables_umaps.py 
            --infile %(mudata_obj)s \
            --layers "%(layer_vars)s" 
            --basis_dict "%(embedding_dict)s" 
            --continuous_variables "%(cont_vars)s" 
            --base_figure_dir ./
            --fig_suffix continuous_vars.png
    """
    cmd += " > %(log_file)s "
    job_kwargs["job_threads"] = PARAMS['resources_threads_high']
    P.run(cmd, **job_kwargs)


do_plot_metrics = any([PARAMS['do_plots_categorical_barplots'],
                PARAMS['do_plots_categorical_stacked_barplots'],
                PARAMS['do_plots_continuous_violin']])
# write metadata out
@active_if(do_plot_metrics)
@originate(PARAMS['sample_prefix'] + "_cell_metadata.tsv")
def write_obs(cmtd):
    cmd = """
        python %(py_path)s/write_metadata.py
        --infile %(mudata_obj)s
        --outfile %(cmtd)s > logs/write_obs.log
        """
    job_kwargs["job_threads"] = PARAMS['resources_threads_high']
    P.run(cmd, **job_kwargs)
  

# Plot cluster metrics
@follows(set_up_dirs)
@active_if(do_plot_metrics)
@transform(write_obs, formatter(), "logs/plot_metrics.log")
def plot_metrics(mtd, log_file):
    cmd = """
    Rscript %(r_path)s/plot_metrics.R \
        --mtd_object %(mtd)s \
        --params_yaml pipeline.yml > %(log_file)s
    """
    job_kwargs["job_threads"] = PARAMS['resources_threads_low']
    P.run(cmd, **job_kwargs)


# parse the scatter files
scatter_files = [PARAMS['paired_scatters'], PARAMS['custom_markers_files']['paired_scatters']]
scatter_files = [sc if type(sc) is list else [sc] for sc in scatter_files if sc is not None]
scatter_files = list(set(chain(*scatter_files)))
@follows(set_up_dirs)
@active_if(PARAMS['do_plots_paired_scatters'])
@active_if(len(scatter_files) != 0)
@transform(scatter_files,
            formatter(), "logs/scatters__{basename[0]}.log")
def plot_scatters(infile, log_file):
    print(infile)
    cmd = """
        python %(py_path)s/plot_features_scatter.py
        --mdata_object %(mudata_obj)s
        --layers_dict "%(custom_markers_layers)s" 
        --scatters_csv %(infile)s
        > %(log_file)s
        """
    job_kwargs["job_threads"] = PARAMS['resources_threads_high']
    P.run(cmd, **job_kwargs)    


@follows(plot_custom_markers_per_group, plot_custom_markers_umap,
 plot_continuous_umaps, plot_categorical_umaps, plot_metrics, plot_scatters)
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
