# Preprocessing spatial data

The `preprocess_spatial` workflow filters the data and preprocesses the data by normalization, HVG selection, and PCA computation. Multiple `SpatialData` objects of the same assay (`Visium`, `Vizgen`, or `Xenium`) can be filtered and preprocessed in one run.

## Steps

If multiple `SpatialData` objects are provided, the following steps are run for each **with the same parameter setting.**

- `SpatialData` object is filtered by the specified thresholds in the pipeline.yml.  Note, that the filtering step is **optional**. You can avoid filtering by setting the `run` parameter in the pipeline.yml under `filtering` to `False`.
- Post-filter plotting is performed (only when data was filtered, i.e. `run: True`). Specified metrics in the pipeline.yml are plotted in violin and spatial embedding plots. Plots are saved into the `./figures/spatial` directory.
- Data is normalized  and HVGs are selected.
  Before normalization, raw counts are saved into `.layers["raw_counts"]`, if not present already. Normalized counts are saved into `.X` and `.layers["lognorm"]` or `.layers["norm_pearson_resid"]`, depending on the chosen normalization. HVGs are saved into `.var["highly_variable"]`.
- PCA is computed and plotted. PCA plots are also saved into the `./figures/spatial` directory.
- Final `SpatialData` object is saved into the `./filtered.data` directory as a `zarr` file

If multiple `SpatialData` objects have been preprocessed, you have the option to concatenate all of them in the last step of the preprocessing: 

- (Optional) All `SpatialData` objects in `./filtered.data` are concatenated and saved to `./concatenated.data/concatenated.zarr`

## Steps to run

1. Activate conda environment `conda activate pipeline_env`
2. Generate yaml and log file `panpipes preprocess_spatial config`
3. Specify the parameter setting in the pipeline.yml file
4. Run complete preprocess workflow with `panpipes preprocess_spatial make full --local`

The [Preprocessing spatial data](https://panpipes-tutorials.readthedocs.io/en/latest/preprocess_spatial_data/preprocess_spatial_data_with_panpipes.html) tutorial guides you through the preprocessing step by step.
