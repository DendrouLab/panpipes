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

from panpipes.funcs.io import write_obs


parser = argparse.ArgumentParser()
# required option
parser.add_argument("--sampleprefix",
                    default="",
                    help="prefix to prepend when saving the metadata file")
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

L.info("running scanpy spatial qc pipeline")

sc.settings.verbosity = 3


L.info("running with args:")
L.info(args)
figdir = args.figdir

if not os.path.exists(figdir):
    os.mkdir(figdir)

sc.settings.figdir = figdir
sc.set_figure_params(scanpy=True, fontsize=14, dpi=300, facecolor='white', figsize=(5,5))


L.info("reading in data")

mdata = mu.read(args.input_anndata)
spatial = mdata['spatial']

L.info("spatial data is")
print(spatial)
L.info("sample id")
print(spatial.obs["sample_id"])

qc_vars = []

if args.customgenesfile is not None:
    if os.path.exists(args.customgenesfile):
        cat_dic = {}
        customgenes = pd.read_csv(args.customgenesfile)
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
                spatial.var[xname] = [x in gene_list for x in spatial.var_names] 
                qc_vars.append(xname)
        else:
            L.info("no genes to calc proportions for")
           
        # Score genes 
        if args.score_genes is not None:
            score_genes = args.score_genes.split(",")
            score_genes = [a.strip() for a in score_genes]
            for kk in score_genes:
                xname= kk
                gene_list = cat_dic[kk]
                sc.tl.score_genes(spatial, gene_list , 
                                    ctrl_size=min(len(gene_list), 50), 
                                    gene_pool=None, 
                                    n_bins=25, 
                                    score_name=kk + '_score', 
                                    random_state=0, 
                                    copy=False, 
                                    use_raw=None)
        else:
            L.info("no genes to calc scores for")    
        
    else:
        sys.exit("the path of the custom genes file does not exist")

        
        
        
# Calculate QC metrics 
L.info("calculating QC metrics")
percent_top = [50, 100, 200, 500] #default
percent_top = [x for x in percent_top if x <= spatial.n_vars]
sc.pp.calculate_qc_metrics(spatial, qc_vars=qc_vars, percent_top=percent_top, inplace=True)

if args.spatial_filetype == "vizgen":
    spatial.obsm["blank_genes"].to_numpy().sum() / spatial.var["total_counts"].sum() * 100

# Calculate cc scores 
if args.ccgenes is not None:
    ccgenes = pd.read_csv(args.ccgenes, sep='\t')
    sgenes = ccgenes[ccgenes["cc_phase"] == "s"]["gene_name"].tolist()
    g2mgenes = ccgenes[ccgenes["cc_phase"] == "g2m"]["gene_name"].tolist()
    L.info("calculating CC scores")
    sc.tl.score_genes_cell_cycle(spatial, s_genes=sgenes, g2m_genes=g2mgenes)

# Aug 2023: we now need to update the mdata object to pick the calc proportion outputs made on 
# spatial = mdata['spatial']
print(spatial.obs.columns)
print(mdata.obs.columns)

mdata.update()

single_id = os.path.basename(str(args.input_anndata))
single_id.replace("_raw.h5mu","")
L.info("updated metadata")
print(mdata.obs.columns)

L.info("saving obs in a metadata tsv file")
write_obs(mdata, output_prefix=(args.sampleprefix +"."+ single_id), output_suffix="_cell_metadata.tsv") #check this function writes one metadata for each sample input
L.info("saving mudata")
mdata.write(args.outfile)

L.info("done")

