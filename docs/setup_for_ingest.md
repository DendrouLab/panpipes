
# Inputs to QC_MM pipeline

To run the `ingest` workflow, the minimum required columns are:

sample_id | rna_path | rna_filetype  
----------|----------|-------------

If you want to analyse other modalities, add columns to the input file and specify the relevant file path or filetype

- prot_path/prot_filetype
- atac_path/atac_filetype
- tcr_path/tcr_filetype
- bcr_path/bcr_filetype

See [Supported input filetypes](##Supported-input-filetypes) to see the options for the {X}_filetype columns

For examples of ingestion files, see example at `resources/sample_file_ingest.txt` and the other sample files included in the tutorials and the resources folder.

**Please note:**  
If you have cellranger outputs which have rna and prot (a Citeseq experiment) or rna and atac (a multiome experiment) within the same files, you have to specify the same path in rna_path and prot_/atac_path.

If using as input a cellranger output folder and {X}_filetype "cellranger", specify the path to the folder containing all the cellranger outputs. This is normally referred to as the "outs" folder, which includes multiple output files including h5 files, bam files and filtered_bc or raw_bc matrices folders. 
The path should be specified until the outs folder `projectX_folder/data.dir/outs`. 
For all the other inputs, including 10X h5 objects, the path should be the complete path to the file. 

To include sample level metadata, you can add additional columns to the submission file
e.g Tissue and Diagnoisis columns in `resources/sample_file_ingest.txt`
You will also need to list which additional metadata columns you want to include in your data object in the pipeline.yml for ingest.

## Additional file inputs for ATAC data

Include additional files from the cellranger outputs under the following three columns:
- fragments_file
- peak_annotation_file
- per_barcode_metrics_file

## Supported input filetypes:

For each modality per sample, we require you to specify the value in the key column in the X_filetype columns

modality    |key       |description
------------|----------|----------
rna/prot/atac|cellranger| the "outs" folder produced by **cellranger count**
rna/prot/atac|cellranger_multi| the "outs" folder produced by **cellranger multi**
rna/prot/atac|10X_h5   | outs/filtered_feature_bc_matrix.h5 produced by cellranger
rna/prot/atac|hd5 | Read a generic .h5 (hdf5) file.
rna/prot/atac|h5ad  | Anndata h5ad objects (one per sample)
rna/prot/atac|txt_matrix  | tab-delimited file (one per sample)
rna/prot/atac|csv_matrix  | comma-delimited file (one per sample)
tcr/bcr     |cellranger_vdj| Path to filtered_contig_annotations.csv, all_contig_annotations.csv or all_contig_annotations.json. produced by **cellranger vdj** further [details](https://scverse.org/scirpy/latest/generated/scirpy.io.read_10x_vdj.html)
tcr/bcr     |tracer| data from [TraCeR](https://github.com/Teichlab/tracer) further [details](https://scverse.org/scirpy/latest/generated/scirpy.io.read_tracer.html)
tcr/bcr     |bracer| data from [BraCeR](https://github.com/Teichlab/bracer) further [details](https://scverse.org/scirpy/latest/generated/scirpy.io.read_bracer.html)
tcr/bcr     |airr  | airr formatted tsv further [details](https://scverse.org/scirpy/latest/generated/scirpy.io.read_airr.html#scirpy.io.read_airr)

For repertoire (tcr/bcr) inputs, panpipes uses the scirpy io functions https://scverse.org/scirpy/latest/api.html 
Review their documentation for more specific details about inputs.

# Combining data sets

Note that if you are combining multiple datasets from different sources the final anndata object will only contain the intersection of the genes
from all the data sets. For example if the mitochondrial genes have been excluded from one of the inputs, they will be excluded from the final data set.
In this case it might be wise to run qc separately on each dataset, and them merge them together to create on h5ad file to use as input for
integration pipeline.

## Barcode level metadata

If you have cell or barcode level metadata such as results from a demultiplexing algorithm, save it in a two column CSV file containing 2 columns; barcode_id, and sample_id. Specify the path to this file in the pipeline.yml file.

The sample_ids in this file must match the sample_id column in the submission file.
