
# Spatial preprocessing yaml

In this documentation, the parameters of the `preprocess_spatial` yaml file are explained. In general, the user can leave parameters empty to use defaults. <br> To understand the individual steps run by the pipeline, please refer to the [spatial preprocess workflow](). 



## 0. Compute resource options

| `resources` |  |
| --- | --- |
| `threads_high` | __`int`__ (default: 1) <br> Number of threads used for high intensity computing tasks. |
| `threads_medium` | __`int`__ (default: 1) <br> Number of threads used for medium intensity computing tasks. For each thread, there must be enough memory to load your mudata and do computationally light tasks. |
| `threads_low` | __`int`__ (default: 1) <br> Number of threads used for low intensity computing tasks. For each thread, there must be enough memory to load text files and do plotting, requires much less memory than the other two.|

|  |  |
| ---- | --- |
| `condaenv` | __`str`__ (default: None) <br> Path to conda environment that should be used to run `Panpipes`. Leave blank if running native or your cluster automatically inherits the login node environment. |


## 1. Input options

With the preprocess_spatial workflow, one or multiple `MuData` objects can be preprocessed in one run. The workflow reads in all `.h5mu` objects of a directory. The `MuData` objects in the directory need to be of the same assay (vizgen or visium). The workflow then runs the preprocessing of each `MuData` object separately with the same parameters that are specified in the yaml file. 

| |  |
| ---- | --- |
| `input_dir` | __`str`__ (not optional) <br> Path to the folder containing all input `h5mu` files. |
| `assay` | __`str` [`"visium"`, `"vizgen"`]__ (default: "visium") <br> Spatial transcriptomics assay of the `h5mu` files in `input_dir`.|


## 2. Filtering options


| `filtering` |  |
| --- | --- |
| `run` | __`bool`__ (default: False) <br> Whether to run filtering. If `False` will not filter the data and will not produce post-filtering plots. |
| `keep_barcodes` | __`str`__ (default: None) <br> Path to a csv-file that has no header containing barcodes you want to keep. Barcodes that are not in the file, will be removed from the dataset. |


____________

------TODO-------




