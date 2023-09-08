Panpipes sample submission file for spatial trascriptomics analysis
===========================

The spatial trascriptomics ingestion pipeline requires a sample submission file which tells the pipeline where to look for input files. 
This is a tab-separated file with a minimum of 4 required columns, and one sample per row.

The minimum required columns are

sample_id | spatial_path | spatial_filetype | spatial_counts
----------|--------------|------------------|---------------



**sample id**: Each row must have a unique sample ID. 

**spatial_path**: If giving a cellranger path, give the path folder containing all the cellranger outputs. For vizgen, this is the the root directory containing Vizgen files.

**spatial_filetype**: The "filetype" column tells panpipe how to read in the data. Panpipes for spatial trascriptomics currently supports `visium` and `vizgen` formats.

**spatial_counts**: The file containing the matrix of features per spot-barcode of counts. Usually **filtered_feature_bc_matrix.h5** or **raw_feature_bc_matrix.h5** for a `visium` dataset. For `vizgen` inputs, this file typically ends with `_cell_by_gene.csv.`



| sample_id | spatial_path             | spatial_filetype   | spatial_counts                |
| --------- |--------------------------|--------------------|-------------------------------|
| Human_LN  | data/V1_Human_Lymph_Node | visium             | filtered_feature_bc_matrix.h5 |

For `vizgen` inputs, an extra column **spatial_metadata** is required.

Optionally, a transformation_file can be supplied as an extra column **spatial_transformation**

| sample_id | spatial_path | spatial_filetype | spatial_counts                          | spatial_metadata                         | spatial_trasformation |
| --------- |--------------|------------------|-----------------------------------------|------------------------------------------|--------------------|
| Sample3   | data.dir     | vizgen           | Slice1_Replicate1_cell_by_gene_S1R1.csv | Slice1_Replicate1_cell_metadata_S1R1.csv | Slice1_Replicate1_images_micron_to_mosaic_pixel_transform.csv |


