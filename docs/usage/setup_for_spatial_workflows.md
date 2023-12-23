# Sample submission file for the ingestion of spatial data

The spatial transcriptomics ingestion workflow requires a sample submission file that specifies the location of the input files. The sample submission file is a tab-separated file with one row per sample. Panpipes currently supports the ingestion of `Visium` and `Vizgen` data.

The 6 columns of the sample submission file are:

**sample id**: Unique sample ID.

**spatial_path**: The root directory containing the data files. Please note, that the folder structure of the root directory needs to be structured as expected by the [squidpy.read.visium](https://squidpy.readthedocs.io/en/stable/api/squidpy.read.visium.html) (for `Visium` data) or [squidpy.read.vizgen](https://squidpy.readthedocs.io/en/stable/api/squidpy.read.vizgen.html) (for `Vizgen` data) functions.

**spatial_filetype**: Either "vizgen" or "visium".

**spatial_counts**: The count matrix file. Usually `filtered_feature_bc_matrix.h5` or `raw_feature_bc_matrix.h5` for a `Visium` dataset. For `Vizgen` inputs, this file typically ends with `_cell_by_gene.csv.`

**spatial_metadata**: The metadata csv-file for `Vizgen` data. Leave empty for `Visium` data.

**spatial_transformation**: The transformation csv-file for `Vizgen` data. This column is **optional** for `Vizgen` data. Leave empty for `Visium` data.

**Note, that the columns, `sample_id`, `spatial_path`, `spatial_filetype`, and `spatial_counts` are required for both, `Visium` and `Vizgen` data. The `spatial_metadata`(required) and `spatial_transformation`(optional) columns are `Vizgen`-specific and should be left empty for `Visium` data.**

### <u>Example submission file</u>

| sample_id | spatial_path | spatial_filetype | spatial_counts                          | spatial_metadata                         | spatial_transformation |
| --------- |--------------|------------------|-----------------------------------------|------------------------------------------|--------------------|
| V1_Human_Heart |./data_visium/V1_Human_Heart |visium |V1_Human_Heart_filtered_feature_bc_matrix.h5 |
| V1_Human_Lymph_Node |./data_visium/V1_Human_Lymph_Node| visium | V1_Human_Lymph_Node_filtered_feature_bc_matrix.h5 |
Mouse_Brain  | ./data_vizgen | vizgen | cell_by_gene_S1R1.csv | cell_metadata_S1R1.csv | images_micron_to_mosaic_pixel_transform.csv
