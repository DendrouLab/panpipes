# Panpipes - multimodal single cell pipelines 

## Table of contents
- [Introduction](#introduction)
- [Installation](#installation)
- [General principles for running pipelines](#general-principles-for-running-pipelines)
    - [Conda environments](#conda-environments)
- [Running the complete pipeline](#running-the-complete-pipeline)
    - [QC](#qc)
    - [Integration](#integration)
    - [Clustering](#clustering)
- [Inputs to QC pipelines](#inputs-to-qc-pipelines)
  - [Demultiplexing data](#demultiplexing-data)
- [Running pipeline modules separately](#running-pipeline-modules-separately)
  - [Integration](#integration-1)
  - [Clustering](#clustering-1)
- [Future Plans](#future-plans)



# Introduction
These pipelines use cgat-core pipeline software

Available pipelines:
- qc_mm (recommended qc pipeline)
- integration
- clustering


Maintained by Charlotte Rich-Griffin and Fabiola Curion
Contributors: Charlotte Rich-Griffin, Tom Thomas and Fabiola Curion

# Installation

Oxford BMRC Rescomp users find additional advice in [docs/installation_rescomp](https://github.com/DendrouLab/sc_pipelines/blob/master/docs/installation_rescomp.md)

### Python dependencies
It is advisable to run everything in a virtual environment either pip or conda.

Navigate to where you want to create your virtual environment  and follow the steps below to create a virutal environment
```
mkdir panpipes
cd panpipes
python3 -m venv --prompt=panpipes python3-venv-panpipes/
# This will create a panpipes/venv folder
```

activate the environment

```
cd my-project
source python3-venv-panpipes/bin/activate
```

**OR** with conda

```
conda create --name pipeline_env python=3.8.6
conda activate pipeline_env
```

we include an environment.yml for a conda environment tested on all the pipelines packaged in this version of Panpipes.


Download and install this repo

```
git clone https://github.com/DendrouLab/sc_pipelines
cd sc_pipelines
pip install --upgrade pip
pip install -r requirements.txt # installs required python packages
python setup.py develop
```
The pipelines are now installed as a local python package.


### installing R requirements
The pipelines use R (mostly for ggplot visualisations). 

If you are using a venv virtual environment,  the pipeline will call a local R installation, so make sure R is installed and install the required pacakges by running the run the code in  `R/install_R_libs.R`

If you are using a conda virtual environment, R *and the required packages (check this)* will be installed along with the python packages. 

### pipeline configuration
Create a yml file for the cgat core pipeline software to read

```
vim ~/.cgat.yml
```
containing at least the following information
```
cluster:
  queue_manager: sge
condaenv:
```
See [cgat-core documentation](https://cgat-core.readthedocs.io/en/latest/getting_started/Cluster_config.html) for cluster specific additional configuration instructions.

Note that there is extra information on the .cgat.yml for Oxford BMRC rescomp users in [docs/installation_rescomp](https://github.com/DendrouLab/sc_pipelines/blob/master/docs/installation_rescomp.md)


To check the installation was successful run the following line
```
sc_pipelines --help
```
A list of available pipelines should appear!

# General principles for running pipelines

Run the pipeline from the login node on your server, it will use in built the job submission system to submit jobs.

Navigate to the directory where you want to run your analysis (this should not be within the dendrou_pipelines folder, or your virutal environment folder)

```
mkdir data_dir/
cd data_dir/
sc_pipelines qc_mm config
```
This will produce two files, `pipeline.log` and `pipeline.yml`

Edit `pipeline.yml` as appropriate for your data, following the instructions within the yml file.

Then check which jobs will run with the command
```
sc_pipelines qc_mm show full
```
The output of this will show a list of tasks that will be run as part of the pipeline.

To run use the command
```
sc_pipelines qc_mm make full
```


Occasionally you might want to run tasks individually (e.g. to assess outputs before deciding the parameters for the next step)
In order to do this you can run any task in the `show full` list such as:

```
sc_pipelines qc_mm make plot_tenx_metrics
```

### Conda environments

If one or more conda environments are needed to run each of the pipelines, (i.e. one pipeline = one environment) the environment (s) should be specified in the .cgat.yml file or in the pipeline.yml configuration file and it will be picked up by the pipeline as the default environment.
If no environment is specified, the default behaviour of the pipeline is to inherit environment variables from the node where the pipeline is run. However there have been reported issues on SLURM clusters where this was not the default behaviour. In those instances we recommend to add the conda environment param in the .cgat.yml file or in each of the pipeline.yml independently.

i.e. :

```
condaenv: pipeline_env
cluster:
    queue_manager: slurm
    queue: cpu_p
    options: --qos=xxx --exclude=compute-node-0[02-05,08-19],compute-node-010-0[05,07,35,37,51,64,68-71]
```


#  Running the complete pipeline

Run each of pipeline qc, integration and clustering in separate folders.
### QC 
1. Generate config file (`sc_pipelines qc_mm config`) and inputs 
2. Run complete qc pipeline with `sc_pipelines qc_mm make full `
3. Use outputs to decide filtering thresholds. Note that the actual filtering occurs in the first step of integration pipeline

### Integration
1. In a new folder, generate config file for integration, `sc_pipelines integration config` and edit the pipeline.yml file.
2. Run `sc_pipelines integration make plot_pcas` and assess the post filtering qc plots, and pca outputs
3. Run batch correction with `sc_pipelines integration make batch_correction` (or run steps 2 and 3 in one go with `sc_pipelines integration make full`)
4. Use outputs to decide on the best batch correction method
5. Edit the integration pipeline yml with your preferred batch correction 
6. Run `sc_pipelines integration make merge_batch_correction`

### Clustering
TODO: update for multimodal

1. In a new folder, generate config file for integration, `sc_pipelines clustering config` and edit the pipeline.yml file. 
2. Run the clustering pipeline  `sc_pipelines clustering make cluster_analysis`. This will do the initial nearest neighbours and clustering for the parameters you specify. 
3. Decide on the best values for k nearest neighbours based on UMAPs and clustree results. Once decided delete the folders for the parameters you don't need and delete those from the pipeline.yml.
4. Find markers for each of your cluster resolutions with `sc_pipelines clustering make marker_analysis` 
(Again you could run all the clustering pipeline at once with `sc_pipelines clustering make full` but by making decisions along the way you'll reduce the computation and file size burden of the pipeline) 



# Inputs to QC pipelines

For qc_mm the minimum required columns are

sample_id   |    gex_path    |    gex_filetype  
---|---|---

#TODO: Example at `resources/sample_file_qc_mm_gex_only.txt`

where type is one of `cellranger`, `10X_h5`, `h5ad`, `csv_matrix` or `txt_matrix`. If giving a cellranger path, give the path folder containing all the cellranger outputs. Otherwise path should be the complete path to the file. 

If you want to analyse other modalities, add columns to the input file
e.g. adt
- adt_path
- adt_filetype
example at `resources/sample_file_qc_mm.txt`

If you have cellranger outputs which have gex and adt within the same files, specify the same path in gex_path and adt_path

adding tcr:
- tcr_path

To inclue sample level metadata, add additional columns to the submission file e.g. batch.
You will also need to list which additional metadata columns you want to include in your data object in the pipeline.yml for qc_mm.
e.g Tissue and Diagnoisis in `resources/sample_file_qc_mm.txt`

## Barcode level metadata TODO: fix this.
If you have barcode level metadata (e.g. demultiplexing data) you can include two extra columns in your samples 
file (e.g. resources/sample_file_inc_demultiplexing.txt); demultiplex_map_file  and   
demultiplex_mtd_file examples at resource/demult_map.csv and resources/demult_mtd.csv, 
espectively. demultiplexing_map_file should contain all barcodes, "antibody" 
demultimplexing_mtd_file should contain "antibody" which corresponds to the map file and any other metadata associated to your samples
Other metadata that you want to add to the amnndata object can be specified in 
the pipeline.yml.

# Combining data sets.
Note that if you are combining multiple datasets from different sources the final anndata object will only contain the intersection of the genes
from all the data sets. For example if the mitochondrial genes have been excluded from one of the inputs, they will be excluded from the final data set.
In this case it might be wise to run qc separately on each dataset, and them merge them together to create on h5ad file to use as input for
integration pipeline.

# Running pipeline modules separately

For circumstances where you arelady have a qc'd anndata object

## Integration
It is possible to run the integration pipeline starting from one combined anndata object containing all your samples containing raw counts in the X slot, 
either with or without running filtering first.
If your data is set to run simply call your anndata object [PROJECT_PREFIX]_filt.h5ad and set filtering_run: False.
You must have a column called sample_id which groups the data in some way (otherwise the plotting script will break) TODO: Fix this

If you have not filtered your data then you can set run filtering_run: True, and set the remaining parameters, 
BUT you must ensure that your obs columns names which you want to use for filtering match the column names in [resources/qc_column_names.txt](https://github.com/DendrouLab/sc_pipelines/blob/matrix_start/resources/qc_column_names.txt)

## Clustering
To run clustering_scanpy without the prior steps, you will need to produce 2 anndata objects
[PROJECT_PREFIX]_log1p.h5ad and [PROJECT_PREFIX]_scaled.h5ad 

[PROJECT_PREFIX]_log1p.h5ad  
- log normalised data in the adata.X slot
- highly variabel genes calculated

[PROJECT_PREFIX]_scaled.h5ad
- log normalised data saved in adata.raw.X
- scaled data (optionally regress) in adata.X
- pca

Minimal code:
```
sc.pp.normalize_total(adata, target_sum=1e4);
sc.pp.log1p(adata))
sc.pp.highly_variable_genes(adata)
adata.write("anndata_log1p.h5ad")

# optional steps
# adata = adata[:, adata.var.highly_variable]
# sc.pp.regress_out()

sc.pp.scale(adata)
sc.tl.pca(adata)
adata.write("anndata_scaled.h5ad")
``` 

If you want to use a specific batch correction then fit it into the above minimal code as appropriate

(It's probably easier to just run integration again with `filtering_run: False`)




CRG 2022-03-14

