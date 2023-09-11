"""
CGAT pipeline for clustering single cell data with Scanpy.
# ASSUMED INPUT: mudata with normalized and scaled data in the relevant layers
# This pipeline is designed to follow on from pipeline_integration.py

"""

from ruffus import *
import sys
import os
from cgatcore import pipeline as P
import pandas as pd
import cgatcore.iotools as IOTools
import re
from itertools import chain, product
import glob

from panpipes.funcs.processing import extract_parameter_from_fname
PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "pipeline.yml"])

PARAMS['py_path'] =  os.path.join(os.path.dirname(os.path.dirname(__file__)), 'python_scripts')
PARAMS['r_path'] = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'R_scripts')
PARAMS['mudata_with_knn'] = 'mudata_w_neighbors.h5mu'

job_kwargs={}
if PARAMS['condaenv'] is not None:
    job_kwargs["job_condaenv"] =PARAMS['condaenv']


@originate("logs/setup_dirs.sentinel")
def set_up_dirs(log_file):
    os.mkdir("logs")
    mods =  [key for key, value in PARAMS['modalities'].items() if value is True]
    if PARAMS["multimodal"]["run_clustering"] is True :
        mods.append("multimodal")
    for mod in mods:
        os.makedirs(os.path.join(mod, "figures"))
    IOTools.touch_file(log_file)
    pass
## ------------------------------------
## Single modality scripts
## ------------------------------------

# -----------------------------------=
# neighbors
# --------------------------------------
@follows(set_up_dirs)
@originate(PARAMS['mudata_with_knn'])
def run_neighbors(outfile):
    mods =  [key for key, value in PARAMS['modalities'].items() if value is True]
    if  any([PARAMS['neighbors'][mod]['use_existing'] is False for mod in mods]):
        # this means we want to rerun neighbors for at least one assay
        #we want to replace thhe scaled obj with the new neighbors
        log_file="logs/run_single_mod_neighbors.log"
        cmd="""
        python %(py_path)s/rerun_find_neighbors_for_clustering.py \
            --infile %(scaled_obj)s \
            --outfile %(outfile)s  \
            --neighbor_dict '%(neighbors)s' \
            --n_threads %(resources_threads_high)s
            """
        cmd += " > %(log_file)s"
        job_kwargs["job_threads"] = PARAMS['resources_threads_high']
        P.run(cmd, **job_kwargs)
    else:
        P.run('ln -s %(scaled_obj)s %(outfile)s', without_cluster=True)




# ------------------------------------
# UMAP
# ------------------------------------

def gen_umap_jobs():
    """
    Generate find neighbor jobs with all parameter combinations.
    """
    # same infile for all jobs
    # define files based on jobs
    infile = PARAMS['mudata_with_knn']
    mods =  [key for key, value in PARAMS['modalities'].items() if value is True]
    if PARAMS["multimodal"]["run_clustering"] is True :
        mods.append("multimodal")
    for mod in mods:
        for md in PARAMS['umap'][mod]['mindist']:
            if PARAMS['umap'][mod]['mindist'] is not None:
                output_file = os.path.join(mod,  'md' + str(md) +"_umap.txt.gz")
                log_file = os.path.join("logs","_".join([mod,  'md' + str(md) +"_umap.log"]))
                yield [infile, output_file, mod, md, log_file]

@follows(run_neighbors)
@follows(set_up_dirs)
@files(gen_umap_jobs)
def calc_sm_umaps(infile, outfile, mod, mindist, log_file):
    prefix = os.path.split(infile)[0]
    cmd = """
        python %(py_path)s/run_umap.py \
            --infile %(infile)s \
            --outfile %(outfile)s \
            --min_dist %(mindist)s 
            """
    if mod is not None and mod != "multimodal":
        # if this is not specified it will use gex default
        cmd += " --modality %(mod)s"
    elif mod=="multimodal":
        if PARAMS['multimodal_integration_method'].lower() == "wnn":
            cmd += " --neighbors_key wnn"
    cmd += " > %(log_file)s"
    job_kwargs["job_threads"] = PARAMS['resources_threads_high']
    P.run(cmd, **job_kwargs)


# ------------------------------------
# Clustering
# ------------------------------------

def gen_cluster_jobs():
    """
    Generate find neighbor jobs with all parameter combinations.
    """
    # same infile for all jobs
    # define files based on jobs
    infile = PARAMS['mudata_with_knn']
    mods =  [key for key, value in PARAMS['modalities'].items() if value is True]
    if PARAMS['multimodal']['run_clustering'] is True:
        mods.append("multimodal")
    for mod in mods:
        for res in PARAMS['clusterspecs'][mod]['resolutions']:
            if PARAMS['clusterspecs'][mod]['resolutions'] is not None:
                alg = PARAMS['clusterspecs'][mod]['algorithm']
                output_file = os.path.join(mod, 'alg' + alg + '_res' + str(res), "clusters.txt.gz")
                log_file = os.path.join("logs", "_".join([mod, 'alg' + alg + '_res' + str(res), "clusters.log"]))
                yield [infile, output_file, mod, res, alg, log_file]

@follows(set_up_dirs)
@files(gen_cluster_jobs)
@follows(run_neighbors)
def calc_cluster(infile, outfile,  mod, res, alg, log_file):
    cmd = """python %(py_path)s/run_clustering.py 
            --infile %(infile)s 
            --outfile %(outfile)s 
            --resolution %(res)s 
            --algorithm %(alg)s
    """ 
    if mod is not None and mod != "multimodal":
        # if this is not specified it will use gex default
        cmd += " --modality %(mod)s"
    elif mod=="multimodal":
        if PARAMS['multimodal_integration_method'].lower() == "wnn":
            cmd += " --neighbors_key wnn"
    cmd += " > %(log_file)s"
    job_kwargs["job_threads"] = PARAMS['resources_threads_medium']
    P.run(cmd, **job_kwargs)


@collate(calc_cluster,
         regex("(.*)/(.*)/clusters.txt.gz"),
         r"\1/all_res_clusters_list.txt.gz")
def aggregate_clusters(infiles, outfile):
    print(infiles)
    print(outfile)
    infiles_str = ','.join(infiles)
    cmd = "python %(py_path)s/aggregate_csvs.py \
               --input_files_str %(infiles_str)s \
               --output_file %(outfile)s \
               --clusters_or_markers clusters > logs/aggregate_clusters.log"
    job_kwargs["job_threads"] = PARAMS['resources_threads_low']           
    P.run(cmd, **job_kwargs)


@collate([[calc_cluster], [calc_sm_umaps]], formatter(), PARAMS['sample_prefix'] + "_clustered.h5mu")
def collate_mdata(infiles,outfile):
    cluster_files = infiles[0]
    cluster_mods = [os.path.dirname(os.path.dirname(x)) for x in cluster_files]
    cluster_col = [os.path.basename(os.path.dirname(x)) for x in cluster_files]
    cluster_col = [re.sub("alg", "", x) for x in cluster_col]
    cluster_df = pd.DataFrame(zip(cluster_mods,cluster_col,cluster_files), columns = ['mod', 'new_key', 'fpath'])
    cluster_df.to_csv("cluster_file_paths.csv", index=None)
    umap_files = infiles[1]
    print(umap_files)
    umap_mods = [os.path.dirname(x) for x in umap_files]
    umap_key = [extract_parameter_from_fname(x, "md", prefix="") for x in umap_files]
    umap_key = ["X_umap_mindist_" + str(x) for x in umap_key]
    umap_df = pd.DataFrame(zip(umap_mods,umap_key,umap_files), columns = ['mod', 'new_key', 'fpath'])
    umap_df.to_csv("umap_file_paths.csv", index=None)
    
    cmd = """python %(py_path)s/collate_mdata.py
        --clusters_files_csv cluster_file_paths.csv
        --umap_files_csv umap_file_paths.csv
        --output_mudata %(outfile)s 
    """
    if PARAMS['full_obj'] is None:
        mdata_in = PARAMS['mudata_with_knn']
        cmd += "--input_mudata %(mdata_in)s"
    else:
        cmd += "--input_mudata  %(full_obj)s"
    cmd += " > logs/collate_data.log"
    job_kwargs["job_threads"] = PARAMS['resources_threads_medium']
    P.run(cmd, **job_kwargs)


@transform(collate_mdata, 
            formatter(""),
            'logs/plot_clusters_umaps.log')
def plot_cluster_umaps(infile, log_file,):
    # get associated umap
    mods =  [key for key, value in PARAMS['modalities'].items() if value is True]
    if PARAMS['multimodal']['run_clustering']:
        mods.append("multimodal")
    mods = ','.join(mods)
    cmd = """python %(py_path)s/plot_cluster_umaps.py \
    --infile %(infile)s 
    --modalities '%(mods)s'
    """
    cmd += " >> %(log_file)s"
    job_kwargs["job_threads"] = PARAMS['resources_threads_medium']
    P.run(cmd, jobs_limit=1, **job_kwargs)

@transform(aggregate_clusters, regex("(.*)/all_res_clusters_list.txt.gz"),
            r'logs/\1_clustree.log',
            r'\1/figures/clustree.png', ) 
def plot_clustree(infile, log_file, outfile):
    # convert infiles to comma sep. string
    prefix = re.sub('_dir', '', os.path.dirname(infile))
    # call R
    cmd = "Rscript %(r_path)s/plotclustree.R \
        --infile %(infile)s  \
        --plot_title %(prefix)s \
        --outfile %(outfile)s > %(log_file)s"
    
    job_kwargs["job_threads"] = PARAMS['resources_threads_low']
    P.run(cmd,  **job_kwargs)



# # all the defs
# # def clustering_umbrella
@follows(calc_cluster, collate_mdata, plot_cluster_umaps, plot_clustree)
@originate("logs/cluster_analysis.sentinel")
def cluster_analysis(fname):
    IOTools.touch_file(fname)
    pass


# ------------------------------------
# Markers
# ------------------------------------

@subdivide(calc_cluster, regex("(.*)/(.*)/clusters.txt.gz"),
           r"\1/\2/*_markers.txt",
           r"\1/\2")
def gen_marker_jobs(infile, outfile, base_dir):
    mods =  [key for key, value in PARAMS['modalities'].items() if value is True]
    for md in mods:
        IOTools.touch_file(os.path.join(base_dir, md + "_markers.txt"))

@follows(collate_mdata)
@transform(gen_marker_jobs,
           regex("(.*)/(.*)/(.*)_markers.txt"),
           r"logs/\1_\2_\3_markers.log",
           r"\1/\2/\3_markers",
           r"\1",
           r"\1/\2",
           r"\3"
           )
def find_markers(infile, log_file, outfile_prefix, base_mod, cluster_dir, data_mod):
    """
    Runs scanpy.tl.rank_gene_groups in parallel for each cluster
    """
    anndata_file=PARAMS['sample_prefix'] + "_clustered.h5mu" + "/" + data_mod
    cluster_file = os.path.join(cluster_dir, "clusters.txt.gz")
    cmd = """
    python %(py_path)s/run_find_markers_multi.py 
    --infile %(anndata_file)s
    --cluster_file %(cluster_file)s 
    --output_file_prefix %(outfile_prefix)s 
    --mod %(data_mod)s
    """
    min_cells = PARAMS["markerspecs"][data_mod]["mincells"]
    cmd += " --mincells %(min_cells)s"
    layer_choice = PARAMS["markerspecs"][data_mod]["layer"]
    cmd += " --layer '%(layer_choice)s'"
    testuse = PARAMS["markerspecs"][data_mod]["method"]
    cmd += " --testuse '%(testuse)s'"
    if PARAMS['markerspecs'][data_mod]['pseudo_seurat'] is True:
        min_pct =  PARAMS["markerspecs"][data_mod]["minpct"]
        threshuse =  PARAMS["markerspecs"][data_mod]["threshuse"]
        cmd += """ --pseudo_seurat True 
        --minpct %(min_pct)s 
        --threshuse %(threshuse)s 
        """
    cmd += " > %(log_file)s "
    job_kwargs["job_threads"] = PARAMS['resources_threads_high']
    P.run(cmd, **job_kwargs)




# # limit jobs because reading the h5 file simultaneouly is problematic
# @collate([ [collate_mdata],[find_markers]],
#         formatter(),
#         os.path.join("logs", "marker_heatmap.log"))
# def heatmap_markers(infiles, log_file):
#     """
#     Plots a standard Seurat Heatmap, or a Complex Heatmap
#     """
#     # print(infiles, log_file)
#     marker_files = infiles[1]
#     scaled_obj = infiles[0][0]
#     cmd = """
#     Rscript %(r_path)s/marker_heatmaps.R \
#     --anndata %(scaled_obj)s \
#     --marker_files %(marker_file)s \
#     --docomplex %(plotspecs_docomplex)s \
#     --subgroup %(plotspecs_discrete_variables)s > %(log_file)s"""
#     P.run(cmd, job_threads=PARAMS['resources_threads_high'])

# this transforms gen_marker_jobs instead of find_markers because gen_marker_jobs creates empty marker files
# which are then filled by find_markers.
@follows(collate_mdata)
@follows(find_markers)
@transform(gen_marker_jobs,
            regex(r"(.*)/alg(.*)/(.*)_markers.txt"),
            r"logs/cluster\1_alg\2_exprs_\3_dotplots.log",
            r"\1/alg\2/figures/dotplot_top_markers_\3.png",
            r"\1/alg\2/figures",
             r"\1", r"\2", r"\3")
def plot_marker_dotplots(marker_file, log_file, outfile, 
                         fig_path, cluster_mod, cluster_col, expr_mod):
    """
    Plots some additional marker plots
    Read in the h5mu file, pull the cluster col from the h5mu obs.
    then pull the epxression values from the expr_mod.
    
    """
    data_obj=PARAMS['sample_prefix'] + "_clustered.h5mu"
    # check there is a figures directory
    if not os.path.exists(fig_path):
        os.makedirs(fig_path)
    cmd = """
        python %(py_path)s/plot_scanpy_markers.py \
            --infile %(data_obj)s \
            --modality %(expr_mod)s
            --marker_file %(marker_file)s \
            --figure_prefix %(fig_path)s \
            --n %(plotspecs_top_n_markers)s
            
     """
    if cluster_mod != "multimodal":
        cluster_col = cluster_mod + ":" + cluster_col
    cmd += " --group_col %(cluster_col)s"
    layer_choice = PARAMS["markerspecs"][expr_mod]["layer"]
    cmd += " --layer %(layer_choice)s"
    cmd += " > %(log_file)s "
    job_kwargs["job_threads"] = PARAMS['resources_threads_medium']
    P.run(cmd, **job_kwargs)





@follows(find_markers, plot_marker_dotplots)
@originate("marker_analysis.sentinel")
def marker_analysis(fname):
    IOTools.touch_file(fname)
    pass





# #--------
# # Generate cellxgene
# #--------

# cellxgene_file = PARAMS['sample_prefix'] + "_" + PARAMS['modality'] + "_cellxgene.h5ad"
# # @originate(cellxgene_file, add_inputs(gen_neighbor_jobs))
# @merge(calc_neighbors, cellxgene_file)
# def cellxgene(infiles, cellxgene_file):
#     print(infiles)
#     print(cellxgene_file)
#     desired_match = "nneigh" + str(PARAMS['cellxgene_nneighbours']) + "_pcs" + str(PARAMS['cellxgene_npcs'])
#     print(desired_match)
#     desired_match = re.compile(".*" + desired_match + ".*")
#     infiles = list(filter(desired_match.match, infiles))
#     print(infiles)
#     if len(infiles) > 1:
#         sys.exit("too many matches, ask charlotte to debug")
#     elif len(infiles) == 0:
#         sys.exit("no matches, did you fill in the cellxgene section of the yml?")
#     else:
#         infile = infiles[0]
    
#     umaps = glob.glob(os.path.dirname(infile) +"/*_umap.txt.gz")
#     umaps_str = ",".join(umaps)

#     clusters_file = os.path.dirname(infile) + "/all_res_clusters_list.txt.gz"
#     best_cluster_res = PARAMS['clusterspecs_algorithm'] + "_res_" + str(PARAMS['cellxgene_cluster_res'])
#     cmd="""
#     python  %(py_path)s/generate_cellxgene.py
#     --input_anndata %(infile)s \
#     --output_anndata %(cellxgene_file)s \
#     --sample_prefix %(sample_prefix)s
#     --clusters %(clusters_file)s \
#     --umaps %(umaps_str)s \
#     --best_md %(cellxgene_umap_md)s \
#     --best_cluster_col %(best_cluster_res)s
#     """
#     if PARAMS['use_muon']:
#         cmd += " --use_muon True"
#     if PARAMS['modality']:
#         # if this is not specified it will use gex default
#         cmd += " --modality %(modality)s"
#     print(cmd)
#     P.run(cmd, job_threads=PARAMS['resources_threads_medium'])

# #--------


@follows(cluster_analysis, marker_analysis )
@originate("cleanup_done.txt")
def cleanup(file):
    # remove any ctmp fails
    P.run("rm ctmp*", without_cluster=True)
    # delete empty dirs
    P.run("find ./ -empty -type d -delete", without_cluster=True)


@follows(cleanup)
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
