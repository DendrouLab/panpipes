
Panpipes sample submission file
===========================

*This section covers the ingestion of cell suspension datasets. For spatial trascriptomics ingestion please check [Sample submission file for the ingestion of spatial data](./setup_for_spatial_workflows.md)*

The multimodal quality control pipeline (ingest) requires a sample submission file which is used to ingest the data from various formats into a `MuData` multimodal container that the pipeline can interact with. 
The sample submission file is a tab-separated file with a minimum of 3 required columns, and one sample per row

The minimum required columns are

sample_id | rna_path | rna_filetype  
----------|----------|-------------

If you want to ingest multiple modalities, you need add additional columns to the input file. The order of the columns is irrelevant, but the names are fixed as follow:


- `rna_path` and `rna_filetype` for rna modality
- `prot_path` and `prot_filetype` when you have protein data
- `atac_path` and `atac_filetype` when you have atac data

for the repertoire data
- `tcr_path` and `tcr_filetype`
- `bcr_path` and `bcr_filetype`

## Submission file columns

1. **sample id**: Each row must have a unique sample ID. It doesn't necessarily correspond to folder where the input sample is stored, but will be used from the ingestion onward to identify this sample.
For example, in your dataset with data from six individuals, you might have sample names like 'Sample_1' to 'Sample_6,' but you can choose to name your sample something more meaningful to you, like 'pbmc_control','pbmc_treated'.

1. **{X}_paths**: If giving a cellranger path, give the path folder containing all the cellranger outputs, known as the `outs` folder. Otherwise path should be the complete path to the file. If you have cellranger outputs which have rna and prot within the same files, specify the same path in rna_path and prot_path. The same applies to multiome and cellranger multi outputs.

2. **{X}_filetype**: The "filetype" column tells panpipe how to read in the data. Panpipes supports a range of inputs. See the [supported input filetypes](#supported-input-filetypes) below to see the options for the {X}_filetype columns

## Additional input files when processing atac data

Panpipes supports reading ATAC/multiome data from a single sample. If you have multiple samples you have to merge them before ingesting with panpipes to ensure that the peak coordinates are harmonized across samples, since simple concatenation would induce undesired batch effects. We recommend using [cellranger](https://kb.10xgenomics.com/hc/en-us/articles/6057890578829-Does-cellranger-atac-aggr-redo-peak-calling-).  or [cellatac](https://github.com/cellgeni/cellatac) to aggregate the peaks across samples. Check the [signac documentation](https://stuartlab.org/signac/articles/merging) for a detailed explanation of the steps (you can use this approach to merge samples if you don't have the cellranger outputs)


When available, you can include additional files from the cellranger outputs under the following three columns:

- **peak_annotation_file**: file containing peaks mapped to gene based on the genomic location of the nearby gene. This is a tabular format tab separated (usually `.tsv`). See [general 10X info](https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/pipelines/latest/output/peak-annotations?src=social&lss=youtube&cnm=&cid=NULL)
- **per_barcode_metrics_file** a tabular format `.csv` file containing metrics for every observed barcode. See [general 10X info](https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/pipelines/latest/output/per_barcode_metrics?src=social&lss=youtube&cnm=&cid=NULL)

- **fragments_file** : cellranger outputs this large [BED-like tabular file](https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/pipelines/latest/output/fragments?src=social&lss=youtube&cnm=&cid=NULL), where each line represents a unique ATAC-seq fragment captured by the assay. Normally, the pipeline outs/ folder contains fragments.tsv.gz and a fragments.tsv.gz.tbi index generated with tabix. panpipes expects both there, and it will fail to read the atac modality if you specified a fragment file but the index is missing. You can generate your own index file using [tabix](https://www.htslib.org/doc/tabix.html).

**Note** panpipes supports both tab or comma-separated tables. 

## Custom columns 

### Sample-level metadata

To include sample level metadata, you can add additional columns to the submission file
e.g tissue and tiagnosis columns in the [example](#example-sample-submission-file) below. There is no requirement for the name or the order of these columns.

To make sure these columns are ingested in the `MuData` object, you must include them by specifying them in `metadatacols` in the `ingest` workflow's `pipeline.yml` configuration file.  

```yaml 

    `metadatacols`: Use this section if you wish to include additional columns specified in your submission file, such as 'sex,' 'batch,' 'diseases,' etc. Leave it empty if you don't want to include metadata.
```
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

For repertoire (tcr/bcr) inputs, panpipes uses the scirpy (v.0.12.0) functions [scirpy]<https://scverse.org/scirpy/latest/api.html>. Review their documentation for more specific details about inputs.


## Example sample submission file

| sample_id | rna_path                           | rna_filetype | prot_path                            | prot_filetype | tissue | diagnosis |
|-----------|------------------------------------|--------------|-------------------------------------|--------------|--------|-----------|
| Sample1   | Sample1_gex.csv                    | csv_matrix   | Sample1_adt.csv                     | csv_matrix   | pbmc   | healthy   |
| Sample2   | cellranger_count/Sample2_GEX/outs/ | cellranger   | cellranger_count/Sample2_CITE/outs/ | cellranger   | pbmc   | diseased  |

Download this file: [sample_file_qc_mm.txt](sample_file_qc_mm.txt)
View other examples on sample submission files on our [github page](https://github.com/DendrouLab/panpipes/tree/main/panpipes/resources)


## Barcode-level metadata

If you have cell or barcode level metadata such as the results from a [cell-demultiplexing pipeline](https://hadge.readthedocs.io/en/latest/), you can provide it to panpipes to be ingested with the rest of the data. Save this file in a `.csv` file containing at least 2 columns: 
- **barcode_id**: the column with the cell barcodes 
- **sample_id**: that indicates which sample the cells belong to.
The file can have additional custom metadata columns to be specified in the configuration file.
**Note** The sample_ids in this file **must** match the sample_id column in the submission file. If you are ingesting multiple files (i.e. you have multiple sample_ids) the barcode-level metadata file **must be the concatenated metadata for all the samples** in your sample submission file.

To read in the barcode-level metadata, specify the path to this file in the ingestion `pipeline.yml` configuration file.

```yaml
barcode_mtd:
  include: True
  path: path_to_file
  metadatacols: percent_mito,tissue,lab
```




## A note on genes when combining data sets

Note that if you are combining multiple datasets from different sources the final anndata object will only contain the intersection of the genes
from all the data sets. For example if the mitochondrial genes have been excluded from one of the inputs, they will be excluded from the final data set.
In this case it might be wise to run qc separately on each dataset, and them merge them together to create on h5ad file to use as input for
integration pipeline.

