
Panpipes sample submission file
===========================

*This section covers the ingestion of cell suspension datasets. For spatial trascriptomics ingestion please check [Sample submission file for the ingestion of spatial data](./setup_for_spatial_workflows.md)*

The multimodal QC pipeline (ingest) requires a sample submission file which it uses to ingest the data into the pipeline. This is a tab-separated file with a minimum of 3 required columns, and one sample per row.

The minimum required columns are

sample_id | rna_path | rna_filetype  
----------|----------|-------------

If you want to analyse other modalities, add additional columns to the input file

- `prot_path` and `prot_filetype` (when you have protein data)
- `atac_path` and `atac_filetype` (when you have atac data)
- `tcr_path` and `tcr_filetype`
- `bcr_path` and `bcr_filetype`

1. **sample id**: Each row must have a unique sample ID. It doesn't necessarily correspond to folder where the same is stored, will be used from the ingestion onward to identify your sample. 
For example, in your dataset with data from six individuals, you might have sample names like 'Sample_1' to 'Sample_6,' but you can choose to name your sample something more meaningful to you, like 'pbmc'

2. **{X}_paths**: If giving a cellranger path, give the path folder containing all the cellranger outputs, known as the `outs` folder. Otherwise path should be the complete path to the file. If you have cellranger outputs which have rna and prot within the same files, specify the same path in rna_path and prot_path. The same applies to multiome and cellranger multi outputs.


3. **{X}_filetype**: The "filetype" column tells panpipe how to read in the data. Panpipes supports a range of inputs. See the [supported input filetypes](#supported-input-filetypes) below to see the options for the {X}_filetype columns



## Sample metadata

To include sample level metadata, you can add additional columns to the submission file
e.g Tissue and Diagnosis columns in [sample_file_qc_mm.txt](sample_file_qc_mm)
To include the additional metadata columns, specify them in `metadatacols` in the pipeline.yml for `ingest`.

The order of the columns doesn't matter.

    `metadatacols`: Use this section if you wish to include additional columns specified in your submission file, such as 'sex,' 'batch,' 'diseases,' etc. Leave it empty if you don't want to include metadata.

## Example sample submission file


| sample_id | rna_path                           | rna_filetype | prot_path                            | prot_filetype | tissue | diagnosis |
|-----------|------------------------------------|--------------|-------------------------------------|--------------|--------|-----------|
| Sample1   | Sample1_gex.csv                    | csv_matrix   | Sample1_adt.csv                     | csv_matrix   | pbmc   | healthy   |
| Sample2   | cellranger_count/Sample2_GEX/outs/ | cellranger   | cellranger_count/Sample2_CITE/outs/ | cellranger   | pbmc   | diseased  |


Download this file: [sample_file_qc_mm.txt](sample_file_qc_mm.txt)

View other examples on sample submission files on our [github page](https://github.com/DendrouLab/panpipes/tree/main/panpipes/resources)


## Additional file inputs for ATAC data
Include additional files from the cellranger outputs under the following three columns:
- fragments_file 
- peak_annotation_file
- per_barcode_metrics_file

## Supported input filetypes

For each modality per sample, specify the value in the key column in the X_filetype columns

modality    |key       |description
------------|----------|----------
rna/prot/atac|cellranger| the "outs" folder produced by **cellranger count**
rna/prot/atac|cellranger_multi| the "outs" folder produced by **cellranger multi**
rna/prot/atac|10X_h5   | outs/filtered_feature_bc_matrix.h5 produced by cellranger
rna/prot/atac|hd5 | Read a generic .h5 (hdf5) file.
rna/prot/atac|h5ad  | Anndata h5ad objects (one per sample)
rna/prot/atac|txt_matrix  | tab-delimited file (one per sample)
rna/prot/atac|csv_matrix  | comma-delimited file (one per sample)
tcr/bcr     |cellranger_vdj| Path to filtered_contig_annotations.csv, all_contig_annotations.csv or all_contig_annotations.json.  produced by **cellranger vdj** further [details](https://scverse.org/scirpy/latest/generated/scirpy.io.read_10x_vdj.html)
tcr/bcr     |tracer| data from [TraCeR](https://github.com/Teichlab/tracer) further [details](https://scverse.org/scirpy/latest/generated/scirpy.io.read_tracer.html)
tcr/bcr     |bracer| data from [BraCeR](https://github.com/Teichlab/bracer) further [details](https://scverse.org/scirpy/latest/generated/scirpy.io.read_bracer.html)
tcr/bcr     |airr  | airr formatted tsv further [details](https://scverse.org/scirpy/latest/generated/scirpy.io.read_airr.html#scirpy.io.read_airr)

For repertoire (tcr/bcr) inputs, panpipes uses the scirpy (v.0.12.0) functions [scirpy]https://scverse.org/scirpy/latest/api.html. Review their documentation for more specific details about inputs.


## Combining data sets.
Note that if you are combining multiple datasets from different sources the final anndata object will only contain the intersection of the genes
from all the data sets. For example if the mitochondrial genes have been excluded from one of the inputs, they will be excluded from the final data set.
In this case it might be wise to run qc separately on each dataset, and them merge them together to create on h5ad file to use as input for
integration pipeline.



## Barcode level metadata 
If you have cell or barcode level metadata such as results from a demultiplexing algorithm, save it in a two column csv file containing 2 column; barcode_id, and sample_id. Specify the path to this file in the pipeline.yml file.

The sample_ids in this file must match the sample_id column in the submission file.
