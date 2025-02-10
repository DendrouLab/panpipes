Sample submission file for the ingestion of spatial data
===========================

The spatial transcriptomics ingestion workflow requires a sample submission file that specifies the location of the input files. The sample submission file is a tab-separated file with one row per sample. Panpipes currently supports the ingestion of `Visium`, `Vizgen`, and `Xenium` data. The data of different technologies needs to be ingested separately with different sample submission files. 


The minimum required (non-optional) columns for each submission file are

**sample id**: Unique sample ID.

**spatial_path**: The root directory containing the data files. Please note, that the folder structure of the root directory needs to be structured as expected by the [spatialdata_io.visium](https://spatialdata.scverse.org/projects/io/en/latest/generated/spatialdata_io.visium.html) (for `Visium` data), [spatialdata_io.merscope](https://spatialdata.scverse.org/projects/io/en/latest/generated/spatialdata_io.merscope.html) (for `Vizgen` data), or [spatialdata_io.xenium](https://spatialdata.scverse.org/projects/io/en/latest/generated/spatialdata_io.xenium.html) (for `Xenium` data) functions.

**spatial_filetype**: Either "vizgen", "visium", or "xenium".


## Visium

The 7 columns of the Visium sample submission file are:

sample_id |	spatial_path |	spatial_filetype |	visium_feature_bc_matrix |	visium_fullres_image_file |	visium_tissue_positions_file |	visium_scalefactors_file	
----------|----------|------------|-----------|----------|-------------|-------------

The following 4 columns are **optional**:

**visium_feature_bc_matrix**: Name of the counts file. Corresponds to the `counts_file` parameter of [spatialdata_io.visium](https://spatialdata.scverse.org/projects/io/en/latest/generated/spatialdata_io.visium.html)

**visium_fullres_image_file**: Path to the full-resolution image. Corresponds to the `fullres_image_file` parameter of [spatialdata_io.visium](https://spatialdata.scverse.org/projects/io/en/latest/generated/spatialdata_io.visium.html)

**visium_tissue_positions_file**: Path to the tissue positions file. Corresponds to the `tissue_positions_file` parameter of [spatialdata_io.visium](https://spatialdata.scverse.org/projects/io/en/latest/generated/spatialdata_io.visium.html)

**visium_scalefactors_file**:	Path to the scalefactors file. Corresponds to the `scalefactors_file` parameter of [spatialdata_io.visium](https://spatialdata.scverse.org/projects/io/en/latest/generated/spatialdata_io.visium.html)

#### [Example submission file](https://github.com/DendrouLab/panpipes-tutorials/blob/sarah_spatialData/docs/ingesting_visium_data/sample_file_qc_visium.txt)


## Vizgen

The 6 columns of the Vizgen sample submission file are:  

sample_id |	spatial_path |	spatial_filetype |	vpt_cell_by_gene    |	vpt_cell_metadata	|	vpt_cell_boundaries
----------|----------|------------|----------|-------------|-------------

The following 3 columns are **optional**:

**vpt_cell_by_gene**: The file name of the output of the vizgen-postprocessing-tool. See [spatialdata_io.merscope](https://spatialdata.scverse.org/projects/io/en/latest/generated/spatialdata_io.merscope.html)

**vpt_cell_metadata**: The file name of the output of the vizgen-postprocessing-tool. See [spatialdata_io.merscope](https://spatialdata.scverse.org/projects/io/en/latest/generated/spatialdata_io.merscope.html)

**vpt_cell_boundaries**: The file name of the output of the vizgen-postprocessing-tool. See [spatialdata_io.merscope](https://spatialdata.scverse.org/projects/io/en/latest/generated/spatialdata_io.merscope.html)


#### Example submission files [MERFISH](https://github.com/DendrouLab/panpipes-tutorials/blob/sarah_spatialData/docs/ingesting_merfish_data/sample_file_qc_merfish.txt) [MERSCOPE](https://github.com/DendrouLab/panpipes-tutorials/blob/sarah_spatialData/docs/ingesting_merscope_data/sample_file_qc_merscope.txt)

## Xenium

The 3 columns of the Xenium sample submission file are:

sample_id |	spatial_path |	spatial_filetype |
----------|----------|------------

#### [Example submission file](https://github.com/DendrouLab/panpipes-tutorials/blob/sarah_spatialData/docs/ingesting_xenium_data/sample_file_qc_xenium.txt)






