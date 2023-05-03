import argparse
import yaml
# import scanpy as sc
import pandas as pd
# import numpy as np
# from scipy.sparse import csr_matrix
import muon as mu
import warnings
from muon._atac.tools import add_peak_annotation, locate_fragments


import sys
import logging

L = logging.getLogger()
L.setLevel(logging.INFO)
log_handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s: %(levelname)s - %(message)s')
log_handler.setFormatter(formatter)
L.addHandler(log_handler)

from panpipes.funcs.io import check_filetype, load_mdata_from_multiple_files, update_intersecting_feature_names, read_yaml
from panpipes.funcs.processing import check_for_bool, add_var_mtd, update_var_index, intersection, intersect_obs_by_mod


parser = argparse.ArgumentParser()
#parser.add_argument('--mode_dictionary',default=mode_dictionary,type=json.loads)
parser.add_argument('--mode_dictionary',default="",
                    help="accepts a yaml speficication of modalities dictionary", 
                    type=yaml.safe_load)
parser.add_argument('--sample_id',
                    default=None,
                    help='')
parser.add_argument('--rna_infile',
                    default=None,
                    help='')
parser.add_argument('--rna_filetype',
                    default=None,
                    help='')
parser.add_argument('--prot_infile', 
                    default=None,
                    help='')
parser.add_argument('--prot_filetype', 
                    default=None,
                    help='')
parser.add_argument('--subset_prot_barcodes_to_rna', 
                    default=True,
                    help='')
parser.add_argument('--atac_infile', 
                    default=None,
                    help='')
parser.add_argument('--atac_filetype', 
                    default=None,
                    help='')
parser.add_argument('--tcr_filtered_contigs', 
                    default=None,
                    help='')
parser.add_argument('--tcr_filetype', 
                    default=None,
                    help='')
parser.add_argument('--bcr_filtered_contigs', 
                    default=None,
                    help='')
parser.add_argument('--bcr_filetype',
 default=None,
                    help='')
parser.add_argument('--output_file',
                    default=None,
                    help='')
parser.add_argument('--per_barcode_metrics_file',
                    default=None,
                    help='ATAC/Multiome specific input file from csv')                    
parser.add_argument('--peak_annotation_file',
                    default=None,
                    help='ATAC/Multiome specific input file from csv')                    
parser.add_argument('--fragments_file',
                    default=None,
                    help='ATAC/Multiome specific input file from csv')                    


parser.set_defaults(verbose=True)
args, opt = parser.parse_known_args()
L.info(args)
# Make mdata object

# unimodal mu (check if all the modalities)
if isinstance(args.mode_dictionary, dict):
    mode_dictionary = args.mode_dictionary
else:
    mode_dictionary = read_yaml(args.mode_dictionary) 

permf = [key for key, value in mode_dictionary.items() if value == True]
all_files = {"rna": [args.rna_infile, args.rna_filetype], 
            "prot": [args.prot_infile, args.prot_filetype], 
            "atac": [args.atac_infile, args.atac_filetype], 
            "tcr":[args.tcr_filtered_contigs, args.tcr_filetype],
            "bcr":[args.bcr_filtered_contigs, args.bcr_filetype]}
#subset to the modalities we want from permf
all_files = {nm: x  for (nm, x) in all_files.items() if nm in permf}
# this will work for any combination of modalities or number of modalities in all_files.
[check_filetype(x[0], x[1]) for x in all_files.values()]
mdata = load_mdata_from_multiple_files(all_files)
                
# now lets do some extra processing on the different modalities
if 'rna' in mdata.mod.keys():
    mdata['rna'].var_names_make_unique()
    # keep only rna with at least 1 count (relevent when loading raw data)
    mdata['rna'].obs['total_counts'] = mdata['rna'].X.sum(axis=1) 
    mu.pp.filter_obs(mdata['rna'], 'total_counts', lambda x: (x > 1))
    # drop this columns as we don't need it anymore, and it will confuse things downstream
    mdata['rna'].obs = mdata['rna'].obs.drop(columns="total_counts")
    mdata.update()

if 'prot' in mdata.mod.keys():
    if 'rna' in mdata.mod.keys():
        if check_for_bool(args.subset_prot_barcodes_to_rna):
            intersect_obs_by_mod(mdata, ['rna', 'prot'])
            mdata.update()
    L.info(mdata['prot'].var.head())

if 'atac' in mdata.mod.keys():
    mdata['atac'].var_names_make_unique()
    if args.per_barcode_metrics_file is not None:
        L.info("Returning mudata with metadata from barcode metrics csv")
        per_barcode_metrics_file = pd.read_csv(args.per_barcode_metrics_file, index_col=0) #it's a csv by default and goes in atac obs
        mdata.mod['atac'].obs = mdata.mod['atac'].obs.merge(per_barcode_metrics_file, left_index=True, right_index=True)
    if args.peak_annotation_file is not None:
        L.info("Returning mudata with peak annotation in mdata['atac'].uns['atac']['peak_annotation']")
        add_peak_annotation(mdata, annotation=args.peak_annotation_file, sep="\t", return_annotation= False)
    if args.fragments_file is not None:
        L.info("Returning mudata with fragments file in mdata['atac'].uns['files']['fragments']")
        locate_fragments(mdata,fragments=args.fragments_file, return_fragments= False)


if 'rep' in mdata.mod.keys():
    col_update = mdata['rep'].obs.columns[mdata['rep'].obs.columns.str.contains("count")]
    mdata['rep'].obs[col_update] = mdata['rep'].obs[col_update].apply(pd.to_numeric)
    mdata.update()
    

#make var names unique
for mm in mdata.mod.keys():
    mdata[mm].var_names_make_unique()

mdata.obs['sample_id'] = str(args.sample_id)

# copy the  to each modality
for mm in mdata.mod.keys():
    print("saving sample_id to each modality")
    # mdata[mm].obs['sample_id'] = mdata.obs['sample_id']
    mdata[mm].obs['sample_id'] = mdata.obs.loc[mdata[mm].obs_names,:]['sample_id']

mdata.update()


L.info("saving data")
L.debug(mdata)
# this will write the file as anndata or muon depending on the class of mdata.
mdata.write(args.output_file)

L.info("done")

