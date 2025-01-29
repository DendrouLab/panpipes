'''
scanpy QC spatial
'''
import sys
import logging
L = logging.getLogger()
L.setLevel(logging.INFO)
log_handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s: %(levelname)s - %(message)s')
log_handler.setFormatter(formatter)
L.addHandler(log_handler)

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import pandas as pd
import os
import argparse
import scanpy as sc
import muon as mu
import spatialdata as sd

from panpipes.funcs.io import write_obs


parser = argparse.ArgumentParser()
# required option
parser.add_argument("--input_anndata",
                    default="adata_raw.h5ad",
                    help="")
parser.add_argument("--spatial_filetype",
                    default="",
                    help="")
parser.add_argument("--outfile",
                    default="adata_unfilt.h5ad",
                    help="")
parser.add_argument("--figdir",
                    default="./figures/",
                    help="path to save the figures to")
parser.add_argument("--ccgenes",
                    default=None,
                    help="path to file containing cell cycle genes")
parser.add_argument("--customgenesfile",
                    default=None,
                    help="path to file containing list of genes to quantify")
parser.add_argument("--calc_proportions",
                    default=None,
                    help="which list of genes to use to calc proportion of mapped reads over total,per cell?")
parser.add_argument("--score_genes",
                    default=None,
                    help="which list of genes to use to scanpy.tl.score_genes per cell?")

args, opt = parser.parse_known_args()

sc.settings.verbosity = 3


L.info("Running with params: %s", args)
figdir = args.figdir

if not os.path.exists(figdir):
    os.mkdir(figdir)

sc.settings.figdir = figdir
sc.set_figure_params(scanpy=True, fontsize=14, dpi=300, facecolor='white', figsize=(5,5))


L.info("Reading in SpatialData from '%s'" % args.input_anndata)

#mdata = mu.read(args.input_anndata)
sdata = sd.read_zarr(args.input_anndata)
#spatial = mdata['spatial']

L.info("Spatial data is:")
print(sdata)
L.info("With sample id '%s'" % sdata["table"].obs["sample_id"].unique()[0])

qc_vars = []

if args.customgenesfile is not None:
    if os.path.exists(args.customgenesfile):
        L.info("Reading in custom genes file from '%s'" % args.customgenesfile)
        cat_dic = {}
        customgenes = pd.read_csv(args.customgenesfile)
        if not {'group', 'feature'}.issubset(customgenes.columns):
            L.error("The custom genes file needs to have both columns, 'group' and 'feature'.")
            sys.exit("The custom genes file needs to have both columns, 'group' and 'feature'.")
        custom_cat = list(set(customgenes['group'].tolist()))
        for cc in custom_cat:
            cat_dic[cc] = customgenes.loc[customgenes["group"] == cc,"feature"].tolist()
            
        # define qc_vars 
        if args.calc_proportions is not None:
            calc_proportions = args.calc_proportions.split(",")
            calc_proportions = [a.strip() for a in calc_proportions]
            for kk in calc_proportions:
                xname= kk
                gene_list = cat_dic[kk]
                sdata["table"].var[xname] = [x in gene_list for x in sdata["table"].var_names] 
                qc_vars.append(xname)
           
        # Score genes 
        if args.score_genes is not None:
            score_genes = args.score_genes.split(",")
            score_genes = [a.strip() for a in score_genes]
            for kk in score_genes:
                L.info("Computing gene scores for '%s'" % kk)
                xname= kk
                gene_list = cat_dic[kk]
                sc.tl.score_genes(sdata["table"], gene_list , 
                                    ctrl_size=min(len(gene_list), 50), 
                                    gene_pool=None, 
                                    n_bins=25, 
                                    score_name=kk + '_score', 
                                    random_state=0, 
                                    copy=False, 
                                    use_raw=None)
        
    else:
        L.error("The path of the custom genes file '%s' could not be found" % args.customgenesfile)
        sys.exit("The path of the custom genes file '%s' could not be found" % args.customgenesfile)

        
        
        
# Calculate QC metrics 
qc_info = ""
if qc_vars != []:
    qc_info = " and calculating proportions for '%s'" % qc_vars
L.info("Calculating QC metrics with scanpy.pp.calculate_qc_metrics()" + qc_info)
percent_top = [50, 100, 200, 500] #default
percent_top = [x for x in percent_top if x <= sdata["table"].n_vars]
sc.pp.calculate_qc_metrics(sdata["table"], qc_vars=qc_vars, percent_top=percent_top, inplace=True)

if (args.spatial_filetype == "vizgen") and ("blank_genes" in sdata["table"].obsm):    
    sdata["table"].obsm["blank_genes"].to_numpy().sum() / sdata["table"].var["total_counts"].sum() * 100

# Calculate cc scores 
if args.ccgenes is not None:
    if os.path.exists(args.ccgenes):
        L.info("Reading in cell cycle genes tsv file from '%s'" % args.ccgenes)
        ccgenes = pd.read_csv(args.ccgenes, sep='\t')
        if not {'cc_phase', 'gene_name'}.issubset(ccgenes.columns):
            L.error("The cell cycle genes file needs to have both columns, 'cc_phase' and 'gene_name'.")
            sys.exit("The cell cycle genes file needs to have both columns, 'cc_phase' and 'gene_name'.")
        sgenes = ccgenes[ccgenes["cc_phase"] == "s"]["gene_name"].tolist()
        g2mgenes = ccgenes[ccgenes["cc_phase"] == "g2m"]["gene_name"].tolist()
        L.info("Calculating cell cycle scores")
        sc.tl.score_genes_cell_cycle(sdata["table"], s_genes=sgenes, g2m_genes=g2mgenes)
    else: 
        L.error("The path of the  cell cycle genes tsv file '%s' could not be found" % args.ccgenes)
        sys.exit("The path of the  cell cycle genes tsv file '%s' could not be found" % args.ccgenes)


#TODO: we now need to update the mdata object to pick the calc proportion outputs made on 
# spatial = mdata['spatial']

#mdata.update()

single_id = os.path.basename(str(args.input_anndata))
single_id = single_id.replace("_raw.h5mu","")

L.info("Saving updated obs in a metadata tsv file to ./" + single_id + "_cell_metadata.tsv")
write_obs(sdata["table"], output_prefix=single_id, output_suffix="_cell_metadata.tsv")
L.info("Saving updated SpatialData to '%s'" % args.outfile)
sdata.write(args.outfile)

L.info("Done")

