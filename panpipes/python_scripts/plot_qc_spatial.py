'''
Plot spatial transcriptomics data
'''
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import muon as mu
import os
import argparse
import sys
import logging
import re 
L = logging.getLogger()
L.setLevel(logging.INFO)
log_handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s: %(levelname)s - %(message)s')
log_handler.setFormatter(formatter)
L.addHandler(log_handler)
L.debug("testing logger works")



sc.settings.verbosity = 3

parser = argparse.ArgumentParser()

parser.add_argument("--input_mudata",
                    default="mudata_unfilt.h5mu",
                    help="")
parser.add_argument("--figdir",
                    default="./figures/",
                    help="path to save the figures to")
parser.add_argument("--spatial_filetype",
                    default="",
                    help="")
parser.add_argument("--spatial_qc_metrics",
                    default="None",
                    help="metrics to plot")
parser.add_argument("--grouping_var",
                    default="None",
                    help="variables to group the cells by")



args, opt = parser.parse_known_args()

L.info("running with args:")
L.info(args)

figdir = args.figdir

if not os.path.exists(figdir):
    os.mkdir(figdir)

sc.settings.figdir = figdir
sc.set_figure_params(scanpy=True, fontsize=14, dpi=300, facecolor='white', figsize=(5,5))

mdata = mu.read(args.input_mudata)
spatial = mdata.mod['spatial']

input_data = os.path.basename(args.input_mudata)
pattern = r"_filtered.h5(.*)"
match = re.search(pattern, input_data)
if match is None:
    match = re.search(r"_unfilt.h5(.*)", input_data)
sprefix = input_data[:match.start()]

# convert string to list of strings
qc_metrics = list(args.spatial_qc_metrics.split(","))
group_var = list(args.grouping_var.split(","))


# check if metrics in adata.obs or adata.var
qc_metrics = [metric if metric in spatial.obs.columns or metric in spatial.var.columns else L.warning("Variable '%s' not found in adata.var or adata.obs, will not be plotted" % metric) for metric in qc_metrics]
qc_metrics = [metric for metric in qc_metrics if metric is not None]

# check that group_vars are in adata.obs
group_var = [group if group in spatial.obs.columns else L.warning("group_var '%s' not found in adata.obs, will be ignored" % group) for group in group_var]
group_var = [group for group in group_var if group is not None]
# make sure that it's saved as categorical 
for group in group_var: 
    spatial.obs[group] = spatial.obs[group].astype("category") 
    
if group_var == []:
    group_var = None


# plot the metrics

for metric in qc_metrics: 
    
    # check if in adata.obs:
    if metric in spatial.obs.columns: 
        # check that it's a numeric column, so that it can be plotted: 
        if metric not in spatial.obs._get_numeric_data().columns:
            L.warning("Variable '%s' not numerical in adata.obs, will not be plotted" % metric)
        else:
            if group_var is None: 
                sc.pl.violin(spatial, keys = metric, xlabel = metric+ " in .obs",
                            save =  "_obs_" + metric+ "_" + "."+sprefix + ".png", show = False)
            
            else: #plot violin for each group
                for group in group_var: 
                    sc.pl.violin(spatial, keys = metric,groupby = group, xlabel = group + ", "+ metric+ " in .obs",
                            save = "_obs_" + metric+ "_" + group+ "."+sprefix +".png", show = False)
            #plot spatial 
            sc.pl.embedding(spatial,basis="spatial", color = metric, save = "_spatial_" + metric + "."+sprefix +".png", show = False)

    #check if in adata.var: 
    if metric in spatial.var.columns:
        
        if metric not in spatial.var._get_numeric_data().columns:
            L.warning("Variable '%s' not numerical in adata.var, will not be plotted" % metric)
        else:
            # plot violins 
            ax = sns.violinplot(
                    data=spatial.var[[metric]],
                    orient='vertical', 
                )
            ax.set(xlabel=metric+ " in .var" )
            ax.figure.savefig(figdir + "/" +"violin_var_" + metric + "."+sprefix +".png")



if args.spatial_filetype == "vizgen":
    fig, axs = plt.subplots(1, 4, figsize=(15, 4))

    axs[0].set_title("Total transcripts per cell")
    sns.histplot(
        spatial.obs["total_counts"],
        kde=False,
        ax=axs[0],
    )

    axs[1].set_title("Unique transcripts per cell")
    sns.histplot(
        spatial.obs["n_genes_by_counts"],
        kde=False,
        ax=axs[1],
    )

    axs[2].set_title("Transcripts per FOV")
    sns.histplot(
        spatial.obs.groupby("fov").sum()["total_counts"],
        kde=False,
        ax=axs[2],
    )

    axs[3].set_title("Volume of segmented cells")
    sns.histplot(
        spatial.obs["volume"],
        kde=False,
        ax=axs[3],
    )

    plt.tight_layout()  # Ensures proper spacing between subplots
    #plt.savefig("merfish_histo.png", dpi=300)
    plt.savefig(figdir + "/histograms."+sprefix +".png", dpi=300)  # Adjust dpi as needed
    plt.close()  # Close the figure to free up memory

            

L.info("Done")

