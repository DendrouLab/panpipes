'''
run scrublet on single channel
expects sample id and path to input data
'''
import scrublet as scr
import pandas as pd
import os
import argparse
import muon as mu

import sys
import logging





import matplotlib
matplotlib.use('Agg')
matplotlib.pyplot.ioff()

L = logging.getLogger()
L.setLevel(logging.INFO)
log_handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s: %(levelname)s - %(message)s')
log_handler.setFormatter(formatter)
L.addHandler(log_handler)

parser = argparse.ArgumentParser()
parser.add_argument("--sample_id",
                    default="sampleID",
                    help="name of the sample, usually just the name of the cellranger filtered folder")
parser.add_argument("--inputpath",
                    default="/path/to/anndata_object",
                    help="path to the single channel worth of cells in anndata format")
parser.add_argument("--outdir",
                    default="/gpfs3/well/combat/projects/preprocess/citeseq_final/FC_annotation/dotplot_minimal/",
                    help="string, b or t cells?")
parser.add_argument("--expected_doublet_rate",
                    default=0.06,
                    help="the expected fraction of transcriptomes that are doublets, typically 0.05-0.1. Results are not particularly sensitive to this parameter")
parser.add_argument("--sim_doublet_ratio",
                    default=2,
                    help="the number of doublets to simulate, relative to the number of observed transcriptomes. Setting too high is computationally expensive. Min tested 0.5")
parser.add_argument("--n_neighbors",
                    default=20,
                    help="Number of neighbors used to construct the KNN classifier of observed transcriptomes and simulated doublets. The default value of round(0.5*sqrt(n_cells)) generally works well")
parser.add_argument("--min_counts",
                    default=2,
                    help="Used for gene filtering prior to PCA. Genes expressed at fewer than `min_counts` in fewer than `min_cells` (see below) are excluded")
parser.add_argument("--min_cells",
                    default=3,
                    help="Used for gene filtering prior to PCA. Genes expressed at fewer than `min_counts` (see above) in fewer than `min_cells` are excluded.")
parser.add_argument("--min_gene_variability_pctl",
                    default=85,
                    help="Used for gene filtering prior to PCA. Keep the most highly variable genes (in the top min_gene_variability_pctl percentile), as measured by the v-statistic [Klein et al., Cell 2015]")
parser.add_argument("--n_prin_comps",
                    default=30,
                    help="Number of principal components used to embed the transcriptomes prior to k-nearest-neighbor graph construction")
parser.add_argument("--use_thr",
                    default=True,
                    help="use a user defined thr to define min doublet score to split true from false doublets? if false just use what the software produces")
parser.add_argument("--call_doublets_thr",
                    default=0.25,
                    help="if use_thr is True, this thr will be used to define doublets")
args, opt = parser.parse_known_args()

# loading data
adata = mu.read(args.inputpath + "/rna")


    
counts_matrix = adata.X
L.info('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
L.info('Number of genes in gene list: {}'.format(len(adata.var_names)))

L.info("now initializing the scrublet object with expected_doublet_rate:\n")
L.info(args.expected_doublet_rate)

scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=float(args.expected_doublet_rate))
L.info("predicting doublets with params: \nmincells: %s \nmingenes: %s \nmin_gene_variabilty_pctl: %s\nn_prin_comps: %s\n" % (args.min_counts, args.min_cells, args.min_gene_variability_pctl, args.n_prin_comps))
L.info("predicting using:\n")

doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=int(args.min_counts),
                                                          min_cells=int(args.min_cells),
                                                          min_gene_variability_pctl=float(args.min_gene_variability_pctl),
                                                          n_prin_comps=int(args.n_prin_comps))

# cgat pipelines will probably parse a string
if args.use_thr == "True":
    use_thr = True
else:
    use_thr = False



if use_thr:
    L.info("using provided or default threshold  %s to call doublets, instead of predicted" % args.call_doublets_thr)
    predicted_doublets = scrub.call_doublets(threshold=float(args.call_doublets_thr))
else:

    L.info("attempting to use predicted threshold for doublet prediction")
    try:
        use_threshold=round(scrub.threshold_, 2)
        L.info("using predicted threshold %s" % use_threshold)
    except AttributeError:
        use_threshold=float(args.call_doublets_thr)
        L.info("scrublet couldn't predict a threshold falling back on provided or default threshold %s" % use_threshold)

    predicted_doublets = scrub.call_doublets(threshold=float(use_threshold))
L.info("saving plots to outdir")

fig = scrub.plot_histogram()
fig[0].savefig(args.outdir + '/' + args.sample_id + '_' + "doubletScore_histogram.png",
               bbox_inches='tight', dpi=120)

L.info("cellnames and doublet scores and prediction")

data = pd.DataFrame({'doublet_scores': doublet_scores,
                     'predicted_doublets': predicted_doublets})

data['barcode'] = adata.obs_names
data.to_csv(os.path.join(args.outdir + "/" + args.sample_id + "_scrublet_scores.txt"),
            sep="\t", index=False)
L.info('done')



