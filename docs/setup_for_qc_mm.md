
# Inputs to QC_MM pipeline

For qc_mm the minimum required columns are

sample_id | gex_path | gex_filetype  
----------|----------|-------------


If you want to analyse other modalities, add columns to the input file

- adt_path/adt_filetype
- atac_path/atac_filetype
- tcr_path/tcr_filetype
- bcr_path/bcr_filetype

See [Supported input filetypes](##Supported-input-filetypes) to see the options for the {X}_filetype columns

example at `resources/sample_file_qc_mm.txt`

If giving a cellranger path, give the path folder containing all the cellranger outputs. Otherwise path should be the complete path to the file. 

If you have cellranger outputs which have gex and adt within the same files, specify the same path in gex_path and adt_path

To include sample level metadata, you can add additional columns to the submission file
e.g Tissue and Diagnoisis columns in `resources/sample_file_qc_mm.txt`
You will also need to list which additional metadata columns you want to include in your data object in the pipeline.yml for qc_mm.

## Additional file inputs for ATAC data
Include additional files from the cellranger outputs under the following three columns:
- fragments_file 
- peak_annotation_file
- per_barcode_metrics_file

## Supported input filetypes:

For each modality per sample, specify the value in the key column in the X_filetype columns

modality    |key       |description
------------|----------|----------
gex/adt/atac|cellranger| the "outs" folder produced by **cellranger count**
gex/adt/atac|cellranger_multi| the "outs" folder produced by **cellranger multi**
gex/adt/atac|10X_h5   | outs/filtered_feature_bc_matrix.h5 produced by cellranger
gex/adt/atac|hd5 | Read a generic .h5 (hdf5) file.
gex/adt/atac|h5ad  | Anndata h5ad objects (one per sample)
gex/adt/atac|txt_matrix  | tab-delimited file (one per sample)
gex/adt/atac|csv_matrix  | comma-delimited file (one per sample)
tcr/bcr     |cellranger_vdj| the "outs" folder produced by **cellranger vdj**
tcr/bcr     |tracer| data from [TraCeR](https://github.com/Teichlab/tracer)
tcr/bcr     |bracer| data from [BraCeR](https://github.com/Teichlab/bracer)
tcr/bcr     |airr  | the "outs" folder produced by **cellranger vdj**




# Combining data sets.
Note that if you are combining multiple datasets from different sources the final anndata object will only contain the intersection of the genes
from all the data sets. For example if the mitochondrial genes have been excluded from one of the inputs, they will be excluded from the final data set.
In this case it might be wise to run qc separately on each dataset, and them merge them together to create on h5ad file to use as input for
integration pipeline.



## Barcode level metadata 
If you have cell or barcode level metadata such as results from a demultiplexing algorithm, save it in a two column csv file containing 2 column; barcode_id, and sample_id. Specify the path to this file in the pipeline.yml file.

The sample_ids in this file must match the sample_id column in the submission file.
