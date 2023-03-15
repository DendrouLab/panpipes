import muon as mu
import argparse
import os

import sys
import logging
L = logging.getLogger()
L.setLevel(logging.INFO)
log_handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s: %(levelname)s - %(message)s')
log_handler.setFormatter(formatter)
L.addHandler(log_handler)


parser = argparse.ArgumentParser()
parser.add_argument('--preprocessed_mudata',
                    default='adata_scaled.h5mu',
                    help='')
parser.add_argument('--output_mudata',
                    default='mdata_corrected.h5mu',
                    help='')
parser.add_argument('--rna_correction_choice', default=None,
                    help='')
parser.add_argument('--prot_correction_choice', default=None,
                    help='')
parser.add_argument('--atac_correction_choice', default=None,
                    help='')
parser.add_argument('--multimodal_correction_choice', default=None,
                    help='')


args, opt = parser.parse_known_args()


L.info(args)

base_object = args.preprocessed_mudata

multi_mod = args.multimodal_correction_choice
if multi_mod is not None:
    base_object = "tmp/" + multi_mod + "_scaled_adata.h5mu"


    

uni_mod_paths = {}

if args.rna_correction_choice is not None:
    uni_mod_paths['rna']= "tmp/" + args.rna_correction_choice + "_scaled_adata_rna.h5ad"
if args.prot_correction_choice is not None:
    uni_mod_paths['prot']= "tmp/" + args.prot_correction_choice + "_scaled_adata_prot.h5ad"
if args.atac_correction_choice is not None:
    uni_mod_paths['atac'] = "tmp/" + args.atac_correction_choice + "_scaled_adata_atac.h5ad"
uni_mod_paths

base_obj = mu.read(base_object)


if multi_mod=="totalvi":
    totalvi_extra_obsm = {'prot': ['totalvi_protein_foreground_prob', 'totalvi_denoised_protein'], 'rna': [ 'totalvi_denoised_rna']}
    # we will want to keep these so extract them first
    totalvi_extra_obsm_keep = {}
    for k, v in totalvi_extra_obsm.items():
        print(k,v)
        totalvi_extra_obsm_keep[k] = {w:base_obj[k].obsm[w] for w in v}


# this overwrites the old modalities
for k,v in uni_mod_paths.items():
    mod_replace = mu.read(v)
    if mod_replace.shape[1] == base_obj[k].shape[1]:
        L.info('replacing the whole anndata object with batch corrected')
        base_obj.mod[k] = mod_replace
    else:
        # in this situation the vars in the mod has been subset but the 
        # original object contains all the genes.
        L.info("number of genes does not match base object, only integrating obsp and obsm into full mdata[%s] object" % k)
        base_obj.mod[k].obsm = mod_replace.obsm
        base_obj.mod[k].obsp = mod_replace.obsp
        base_obj.mod[k].uns = mod_replace.uns
        base_obj.mod[k].uns['batch_correct_gene_subset'] = mod_replace.var_names.tolist()
        # merge any extra obs columns
        new_obs = mod_replace.obs.iloc[:, ~ mod_replace.obs.columns.isin( base_obj.mod[k].obs.columns)]
        base_obj.mod[k].obs = base_obj.mod[k].obs.merge(new_obs, left_index=True, right_index=True)

if multi_mod=="totalvi":
    # restore the totalvi extras
    for k, v in totalvi_extra_obsm_keep.items():
        for w in v:
            base_obj[k].obsm[ w] = totalvi_extra_obsm_keep[k][w]

base_obj.update()
base_obj.write(args.output_mudata)

L.info("Done")

