'''
scanpy QC script PROT
order of QC:
- RNA
- PROT
- Repertoire
- ATAC
'''
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import scipy.io
import matplotlib.pyplot as plt
plt.ioff()
import numpy as np
import pandas as pd
import os
import argparse
import scanpy as sc
import seaborn as sns
import muon as mu

from panpipes.funcs.io import  write_obs
from panpipes.funcs.processing import check_for_bool
from panpipes.funcs.scmethods import identify_isotype_outliers
from panpipes.funcs.plotting import adjust_x_axis


import sys
import logging
L = logging.getLogger()
L.setLevel(logging.INFO)
log_handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s: %(levelname)s - %(message)s')
log_handler.setFormatter(formatter)
L.addHandler(log_handler)

parser = argparse.ArgumentParser()
# required option
parser.add_argument("--sampleprefix",
                    default="",
                    help="prefix to prepend when saving the metadata file")
parser.add_argument("--input_anndata",
                    default="adata_unfilt.h5ad",
                    help="")
parser.add_argument("--outfile",
                    default="adata_unfilt.h5ad",
                    help="")
parser.add_argument("--figdir",
                    default="./figures/prot/",
                    help="path to save the figures to")
parser.add_argument('--identify_isotype_outliers',
                    default=False, type=check_for_bool,
                    help='')
parser.add_argument('--isotype_upper_quantile',
                    default=0.9,
                    help='')
parser.add_argument('--isotype_n_pass',
                    default=2,
                    help='')
parser.add_argument('--channel_col',
                    default=None,
                    help='')
parser.add_argument("--per_cell_metrics",
                    default=None,
                    help="")
parser.add_argument("--per_prot_metrics",
                    default=None,
                    help="")


args, opt = parser.parse_known_args()
L.info("Running with params: %s", args)

sc.settings.verbosity = 3
figdir = args.figdir

if not os.path.exists(figdir):
    os.mkdir(figdir)

sc.settings.figdir = figdir
sc.set_figure_params(scanpy=True, fontsize=14, dpi=300, facecolor='white', figsize=(5,5))


L.info("Reading in MuData from '%s'" % args.input_anndata)
mdata = mu.read(args.input_anndata)
prot = mdata['prot']


# determine what QC metrics we want
if args.per_cell_metrics is not None:
    per_cell_metrics = args.per_cell_metrics.split(",")
    per_cell_metrics = [a.strip() for a in per_cell_metrics]



# work out if we already have isotype column, if not try to infer from index.
L.info("Looking for 'isotype' column in .var")
if 'isotype' not in prot.var.columns:
    # this means that isotype column was not included in the protein conversion table 
    # so we are going to have a wwhack at identifying them
    L.info("No 'isotype' column found in .var, adding 'isotype' column.")
    prot.var['isotype'] = prot.var_names.str.contains("isotype")
else:
    # there is already an isotype column in var
    L.info("Isotype column found in var")
    pass

# are there any actual isotypes though?
isotypes = list(prot.var_names[prot.var.isotype].unique())
L.info("Number of isotypes found in data %i" % len(isotypes))
if len(isotypes) > 0:
    L.info("Including isotypes in QC metrics call")
    qc_vars=["isotype"]
else:
    L.info("Not including isotypes in QC metrics call")
    qc_vars=[]
# if isotypes are found then include them in the caluclate qc metrics call.
# this calculates standard metrics 

L.info("Calculating QC metrics with scanpy.pp.calculate_qc_metrics()")
sc.pp.calculate_qc_metrics(prot,
                        var_type="prot",  
                        qc_vars=qc_vars, 
                        percent_top=None,log1p=True, inplace=True)

## let's assess the isotype outlier cells. 
#(Cells with an excessive amount of isotype indicating stickiness)
if (len(isotypes) > 0) & check_for_bool(args.identify_isotype_outliers):
    L.info("Identifying isotype outliers")
    # this measn we found some isotypes earlier
    # we identify outliers on a per channel basis 
    # using the groupby argument
    identify_isotype_outliers(prot, isotypes, 
            quantile_val=float(args.isotype_upper_quantile),
            n_isotypes_pass=int(args.isotype_n_pass), 
            groupby=args.channel_col)

mdata.update()

# write out the cell_metadata
L.info("Saving updated obs in a metadata tsv file to ./" + args.sampleprefix + "_cell_metadata.tsv")
write_obs(mdata, output_prefix=args.sampleprefix, 
        output_suffix="_cell_metadata.tsv")
# write out whole object data.
L.info("Saving updated MuData to '%s'" % args.outfile)
mdata.write(args.outfile)




# now run some per channel metrics - this is not stored in the final object.
# the reason these are not run on the whole object is to get the "var_df" 
# values on a per channel basis 
# IS this cause you need the pandas per channel to plot? YES

# the rest of the script is concerned with Per PROT metrics instead of per cell. 
# per cell metrics are plotted in plotqc.R


channel_col = args.channel_col # TO DO why is this just not simply sample_id 

L.info("Calculating per channel metrics for channel '%s'" % channel_col)

out = {}
for si in prot.obs[channel_col].unique():
    obs_df, var_df = sc.pp.calculate_qc_metrics(prot[prot.obs[channel_col]==si],
                            var_type="prot", 
                            qc_vars=qc_vars, 
                            percent_top=None,log1p=True, inplace=False)
    # store the var df for the metrics we care about in the dictionary.
    out[si] = var_df.copy()
    # add an extra column t the var_df to indicate which channel the data came from
    out[si][channel_col] = si

# concatenate all the var_dfs, so we can summarise the per channel data in boxplots.
var_dat = pd.concat(out).reset_index().rename(columns={'level_0':channel_col, 'level_1': "prot_id"})
prot_id_col=var_dat.columns[1]

if args.per_prot_metrics is not None:
    per_prot_metrics = args.per_prot_metrics.split(",")
    per_prot_metrics = [a.strip() for a in per_prot_metrics]

    # let's make some boxplot figures.
    L.info("Plotting per channel metrics in boxplots")
    for qcp in per_prot_metrics:
        fig, ax = plt.subplots(nrows=1,figsize=(24,6))
        sns.boxplot(data=var_dat, x=prot_id_col, y=qcp, ax=ax)
        adjust_x_axis(ax)
        fig.tight_layout()
        plt.savefig(os.path.join(figdir, "boxplot_" + channel_col + "_" + qcp + ".png"))
    plt.close()

    # write out the per channel metrics in a separate csv.
    L.info("Saving per channel metrics to csv file '"+ args.sampleprefix + "_prot_qc_metrics_per_" + channel_col + ".csv'")
    var_dat.to_csv(args.sampleprefix + "_prot_qc_metrics_per_" + channel_col + ".csv")

L.info("Done")

