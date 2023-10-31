Preprocessing spatial data
==========================

The `preprocess_spatial` workflow filters the data and preprocesses the data by normalization, HVG selection, and PCA computation. Multiple `MuData` objects can be filtered and preprocessed in one run. 

## Steps


- `MuData` object is filtered by the specified thresholds in the pipeline.yml.
- Post-filter plotting is performed. Specified metrics in the pipeline.yml are plotted in violin and spatial embedding plots. Plots are saved into the `./figures/spatial` directory.
- Data is normalized  and HVGs are selected. 
  Before normalization, raw counts are saved into `.layers["raw_counts"]`, if not present already. Normalized counts are saved into `.X` and `.layers["lognorm"]` or `.layers["norm_pearson_resid"]`, depending on the chosen normalization. HVGs are saved into `.var["highly_variable"]`.
  **Note: whether normalization is performed before HVG selection or vice versa depends on the chosen normalization. More information about the normalization and HVG selection can be found in the pipeline.yml**
- PCA is computed and plotted. PCA plots are also saved into the `./figures/spatial` directory.



## Steps to run:

1.  Activate conda environment `conda activate pipeline_env`
2.  Generate yaml and log file `panpipes preprocess_spatial config`
3.  Edit the pipeline.yml file for your dataset
4.  Run complete preprocess workflow with `panpipes preprocess_spatial make full --local`

The [Preprocessing spatial data]() tutorial guides you through the preprocessing step by step. 