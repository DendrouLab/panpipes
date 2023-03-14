'''
scanpy QC script GEX
order of QC:
- GEX
- ADT
- Repertoire
- ATAC
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
parser.add_argument("--scrubletdir",
                    default=None,
                    help="path to save the figures to")
parser.add_argument("--ccgenes",
                    default=None,
                    help="path to file containing cell cycle genes")
parser.add_argument("--customgenesfile",
                    default=None,
                    help="path to file containing list of genes to quantify")
parser.add_argument("--calc_proportions",
                    default="mitochondrial,ribosomal",
                    help="which list of genes to use to calc proportion of mapped reads over total,per cell?")
parser.add_argument("--score_genes",
                    default="MarkersNeutro",
                    help="which list of genes to use to scanpy.tl.score_genes per cell?")

args, opt = parser.parse_known_args()

L.info("Running scanpy gex qc pipeline")

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
rna = mdata['rna']


L.info("merge in the scrublet scores")
# load the scrublet scores into the anndata (if they have been run)
if args.scrubletdir is not None:
    scrub_dir = args.scrubletdir
    sample_ids = rna.obs[['sample_id']].drop_duplicates()
    [scrub_dir + '/' + ss + '_scrublet_scores.txt' for ss in sample_ids.sample_id]
    doubletscores = [pd.read_csv(scrub_dir + '/' + ss + '_scrublet_scores.txt', sep="\t", header=0) for ss in sample_ids.sample_id]
    doubletscores = pd.concat(doubletscores, keys=sample_ids.sample_id).reset_index(level="sample_id")
    # rename the barcodes to match up qwith what the rna.obs barcodes are
    if len(sample_ids) > 1:
        doubletscores['barcode'] = doubletscores['barcode'] + '-' + doubletscores['sample_id']
    doubletscores = doubletscores.set_index('barcode').drop('sample_id', axis=1)
    # merge with rna.obs
    rna.obs = rna.obs.merge(doubletscores, how="left", left_index=True, right_index=True)


# Cell-wise QC based on .var lists of genes 

qc_vars = []

if args.customgenesfile is not None:
    if os.path.exists(args.customgenesfile):
        cat_dic = {}
        customgenes = pd.read_csv(args.customgenesfile)
        custom_cat = list(set(customgenes['group'].tolist()))
        for cc in custom_cat:
            cat_dic[cc] = customgenes.loc[customgenes["group"] == cc,"feature"].tolist()
        if args.calc_proportions is not None:
            calc_proportions = args.calc_proportions.split(",")
            calc_proportions = [a.strip() for a in calc_proportions]
        else:
            L.warn("I don't have genes to calc proportions for")
            
        if args.score_genes is not None:
            score_genes = args.score_genes.split(",")
            score_genes = [a.strip() for a in score_genes]
        else:
            L.warn("I don't have genes to calc scores for")    
        
    else:

        sys.exit("You have not provided a list of custom genes to use for QC purposes")

for kk in calc_proportions:
    xname= kk
    gene_list = cat_dic[kk]
    rna.var[xname] = [x in gene_list for x in rna.var_names] # annotate the group of hb genes as 'hb'
    qc_vars.append(xname)

sc.pp.calculate_qc_metrics(rna, qc_vars=qc_vars, percent_top=None,log1p=True, inplace=True)

if args.score_genes is not None:
    for kk in score_genes:
        xname= kk
        gene_list = cat_dic[kk]
        sc.tl.score_genes(rna, gene_list , 
                            ctrl_size=min(len(gene_list), 50), 
                            gene_pool=None, 
                            n_bins=25, 
                            score_name=kk + '_score', 
                            random_state=0, 
                            copy=False, 
                            use_raw=None)


if args.ccgenes is not None:
    ccgenes = pd.read_csv(args.ccgenes, sep='\t')
    sgenes = ccgenes[ccgenes["cc_phase"] == "s"]["gene_name"].tolist()
    g2mgenes = ccgenes[ccgenes["cc_phase"] == "g2m"]["gene_name"].tolist()
    sc.tl.score_genes_cell_cycle(rna, s_genes=sgenes, g2m_genes=g2mgenes)



L.info("calculated scores and metrics")

L.info("saving anndata and obs in a metadata tsv file")
write_obs(mdata, output_prefix=args.sampleprefix, 
        output_suffix="_cell_metadata.tsv")
# CRITICAL to do WORK OUT WHICH QC SCRIPT TO USE  TO SAVE THE MDATA OR ANNDATA
mdata.write(args.outfile)

L.info("done")

