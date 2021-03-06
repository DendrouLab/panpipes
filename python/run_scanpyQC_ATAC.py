'''
scanpy QC script ATAC
'''
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import scipy.io
import matplotlib.pyplot as plt
plt.ioff()
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import os
import argparse
import seaborn as sns
from muon import atac as ac
import muon as mu
import sys
import logging
L = logging.getLogger()
L.setLevel(logging.INFO)
log_handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s: %(levelname)s - %(message)s')
log_handler.setFormatter(formatter)
L.addHandler(log_handler)
L.debug("testing logger works")

from panpipes.funcs.io import read_anndata, write_anndata, write_obs
from panpipes.funcs.processing import check_for_bool
from panpipes.funcs.scmethods import identify_isotype_outliers
from panpipes.funcs.plotting import adjust_x_axis

sc.settings.verbosity = 3

parser = argparse.ArgumentParser()

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
                    default="./figures/",
                    help="path to save the figures to")
parser.add_argument("--figure_suffix",
                    default="_qc-plot.png",
                    help="figures filename suffix to be appended to figures/umap")
parser.add_argument('--is_paired',
                    default=True,
                    help='If ATAC data comes from a paired multiome experiment or is a standalone ATAC.')
parser.add_argument('--paired_rna',
                    default="rna.h5ad",
                    help='If ATAC data comes from standalone ATAC, do you have an rna anndata with obs.interval that you can use for Tss')
parser.add_argument('--tss_coordinates',
                    default=None,
                    help='If ATAC data comes from standalone ATAC, do you have a TSS annotation file? mm10 or hg19 supported')
parser.add_argument("--atac_qc_metrics",
                    default="n_genes_by_counts,total_counts",
                    help="comma sep list of params to quantify and plot")
parser.add_argument("--use_muon",
                    default="",
                    help="")
                



args, opt = parser.parse_known_args()

L.info("running with args:")
L.debug(args)

figdir = args.figdir

if not os.path.exists(figdir):
    os.mkdir(figdir)

sc.settings.figdir = figdir
sc.set_figure_params(scanpy=True, fontsize=14, dpi=300, facecolor='white', figsize=(5,5))

if args.is_paired:
    args.use_muon= True
    L.info("I'm working on a multiome experiment")
else:
    L.info("qc'ing an atac standalone assay")

mdata = mu.read(args.input_anndata)
atac = mdata.mod['atac']


# CALCULATE SOME GENERAL QC metrics
sc.pp.calculate_qc_metrics(atac, percent_top=None, log1p=False, inplace=True)
ac.tl.nucleosome_signal(atac, n=1e6)
atac.obs['NS']=np.where(atac.obs['nucleosome_signal'] >3, 'trinucleosome', 'mono_trinucleosome')
L.debug("added nucleosome signal to atac")

if args.is_paired:
    #ac.tl.get_gene_annotation_from_rna(mdata['rna'])  # accepts MuData with 'rna' modality or mdata['rna'] AnnData directly if
    tss = ac.tl.tss_enrichment(mdata, n_tss=1000) # by default, features=ac.tl.get_gene_annotation_from_rna(mdata), but it could work with a custom annotation 
    tss.obs['tss_class'] = np.where(tss.obs["tss_score"]>2 , "High","Low") 
    ac.pl.tss_enrichment(tss, color="tss_class")
    plt.savefig(os.path.join(figdir, "tss_enrichment.png"))
    mdata.mod["tss"] = tss # add the tss anndata to the mudata, no point in saving it outside
    
else:
    #since we're here, change some names
    atac.obs["atac_peak_region_fragments"] = atac.obs["peak_region_fragments"]
    if atac.obs["mitochondrial"] is not None:
        atac.obs["atac_mitochondrial_reads"] = atac.obs["mitochondrial"]
    if atac.obs["passed_filters"] is not None:
        atac.obs['atac_fragments']= atac.obs['passed_filters']
    
    if args.paired_rna is not None:
        rna = mu.read(args.paired_rna)
        tss_features = ac.tl.get_gene_annotation_from_rna(rna)
    elif args.tss_coordinates is not None:
        tss_features = pd.read_csv(args.tss_coordinates, index_col=0, sep="\t")
    
    tss_data = ac.tl.tss_enrichment(atac, features= tss_features, n_tss=1000)
    tss_data.obs['tss_class'] = np.where(tss_data.obs["tss_score"]>2 , "High","Low") 
    ac.pl.tss_enrichment(tss_data, color="tss_class")
    plt.savefig(os.path.join(figdir, "tss_enrichment.png"))
    mdata.mod["tss"] = tss_data
    
    if args.paired_rna is None & args.tss_coordinates is None:
        L.warning("""you didn't provide an rna anndata or the features coordinate files to
        run the tss enrichment analysis, this will be skipped
        """)

# calc pct fragments in peaks 
if atac.obs["atac_fragments"] is not None:
    L.debug("adding percent fragments in peaks")
    atac.obs['pct_fragments_in_peaks']=atac.obs['atac_peak_region_fragments']/atac.obs['atac_fragments'] *100          

#PLOTS
mu.pl.histogram(atac, "nucleosome_signal", groupby='NS')
plt.savefig(os.path.join(figdir, "nucleosome.png"))

if args.atac_qc_metrics is not None:
    qc_vars=args.atac_qc_metrics.split(",")
    qc_vars= [a.strip() for a in qc_vars]
else:
    qc_vars = ["n_genes_by_counts","total_counts","pct_fragments_in_peaks"]

qc_vars_plot = [gg for gg in qc_vars if gg in atac.obs.columns]
sc.pl.violin(atac, qc_vars_plot, 
             jitter=0.4, multi_panel=True, save = "atac_metrics_violin.png")

mdata.update()

L.info("saving mudata and obs in a metadata tsv file")
write_obs(mdata, output_prefix=args.sampleprefix, 
        output_suffix="_cell_metadata.tsv")

mdata.write(args.outfile)





#Filter peaks which expression is not detected:
#mu.pp.filter_var(atac, 'n_cells_by_counts', lambda x: x >= 10)
# This is analogous to
#   sc.pp.filter_genes(rna, min_cells=10)
# but does in-place filtering and avoids copying the object

#Filter cells:
#mu.pp.filter_obs(atac, 'n_genes_by_counts', lambda x: (x >= 2000) & (x <= 15000))
# This is analogous to
#   sc.pp.filter_cells(atac, max_genes=15000)
#   sc.pp.filter_cells(atac, min_genes=2000)
# but does in-place filtering avoiding copying the object

#mu.pp.filter_obs(atac, 'total_counts', lambda x: (x >= 4000) & (x <= 40000))
