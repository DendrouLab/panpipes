
# Installation of panpipes

## Step 1: create virtual environment

We recommend running panpipes within a virtual environment to maintain reproducibility


### Option 1: create conda environment (Recommended)

We create a conda environment with R and python
Panpipes has a lot of dependencies, so you may want to consider [`mamba`](https://mamba.readthedocs.io/en/latest/index.html) instead of `conda for installation.

```
conda config --add channels conda-forge
conda config --set channel_priority strict
# you should remove the strict priority afterwards!
conda search r-base
conda create --name pipeline_env python=3.9 r-base=4.3.0
```
now we activate the environment

```
conda activate pipeline_env
```

This follows the suggestions made here: [https://www.biostars.org/p/498049/](https://www.biostars.org/p/498049/) 

Install specific dependencies

```
conda install -c conda-forge pynndescent
```

Install R packages
```
conda install -c conda-forge r-tidyverse r-optparse r-ggforce r-ggraph r-xtable r-hdf5r r-clustree
```

Panpipes requires the unix package `time`, in conda you can install it with:

You can check if it installed with 

```
dpkg-query -W time
```
if this is not already installed on your conda env with: 

```
conda install time
```
or

```
apt-get install time
```

You can install `panpipes` directly from `PyPi` with:

```
pip install panpipes
```

If you intend to use panpies for spatial analysis, instead install:
```
pip install 'panpipes[spatial]'
```
The extra `[spatial]` includes squidpy and cell2location packages.



#### Nightly versions of panpipes.

If you would prefer to use the most recent dev version, install from github

```
git clone https://github.com/DendrouLab/panpipes
cd panpipes
pip install -e .
```

### Option 2: python venv environment:

Navigate to where you want to create your virtual environment  and follow the steps below to create a pip virtual environment

```
python3 -m venv --prompt=panpipes python3-venv-panpipes/
# This will create a panpipes/venv folder
```

activate the environment

```
source python3-venv-panpipes/bin/activate
```

As explained in the conda installation, you can install `panpipes` with:
```
pip install panpipes
```

If you would prefer to use the most recent dev version, install from github

```
git clone https://github.com/DendrouLab/panpipes
cd panpipes
pip install -e .
```



If you are using a venv virtual environment,  the pipeline will call a local R installation, so make sure R is installed and install the required packages with the command we provide below.
(This executable requires that you specify  a CRAN mirror in your `.Rprofile`)

 ```
panpipes install_r_dependencies
 ```

If you want more control over your installation use the [script on github](https://github.com/DendrouLab/panpipes/blob/main/panpipes/R_scripts/install_R_libs.R).
Running with the option `--vanilla` or `--no-site-file` prevents R from reading your `.Renvironment` or `.Rprofile` in case you want to use different settings from you local R installation.
You can expect the installation of R libraries to take quite some time, this is not something related to `panpipes` but how R manages their libraries and dependencies!


#### Check installation

To check the installation was successful run the following line
```
panpipes --help
```
A list of available pipelines should appear!


You're all set to run `panpipes` on your local machine.
If you want to configure it on a HPC server, jump to [step 2](#step-2-pipeline-configuration)


## Step 2 pipeline configuration 

(For SGE or SLURM clusters)
*Note: You won't need this for a local installation of panpipes.*

Create a yml file for the cgat core pipeline software to read

```
vim ~/.cgat.yml
```

For SGE servers:
```
cluster:
  queue_manager: sge
  queue: short
condaenv:
```


For Slurm servers:
```
cluster:
    queue_manager: slurm
    queue: short
```

There are likely other options that relate to **your specific server**, e.g. 
These are added as 
```
cluster:
    queue_manager: slurm
    queue: short
    options: --qos=xxx --exclude=compute-node-0[02-05,08-19],compute-node-010-0[05,07,35,37,51,64,68-71]

```



See [cgat-core documentation](https://cgat-core.readthedocs.io/en/latest/getting_started/Cluster_config.html) for cluster specific additional configuration instructions.

Note that there is extra information on the .cgat.yml for Oxford BMRC rescomp users in [docs/installation_rescomp](https://github.com/DendrouLab/panpipes/blob/main/docs/installation_rescomp.md)

### DRMAA library path


if `echo $DRMAA_LIBRARY_PATH` does not return anything, add DRMAA_LIBRARY_PATH to your bash environment (the line below is specific to the rescomp server)
You need to find the path to this libdrmaa.so.1.0 file on your server, ask your sys admin if you cannot find it!

```
PATH_TO_DRMAA = ""
echo "export DRMAA_LIBRARY_PATH=$PATH_TO/libdrmaa.so.1.0" >> ~/.bashrc
```

### Specifying Conda environments to run panpipes
If using conda environments, you can use one single big environment (the instructions provided do that) or create one for each of the workflows in panpipes, (i.e. one workflow = one environment) 
The environment (s) should be specified in the .cgat.yml global configuration file or in each of the single workflows pipeline.yml configuration files and it will be picked up by the pipeline as the default environment. 
Please note that if you specify the conda environment in the workflows configuration file this will be the first choice to run the pipeline. 



If no environment is specified, the default behaviour of the pipeline is to inherit environment variables from the node where the pipeline is run. However there have been reported issues on SLURM clusters where this was not the default behaviour. In those instances we recommend to add the conda environment param in the .cgat.yml file or in each of the pipeline.yml independently.

i.e. :

```

cluster:
    queue_manager: slurm
    queue: cpu_p
    options: --qos=xxx --exclude=compute-node-0[02-05,08-19],compute-node-010-0[05,07,35,37,51,64,68-71]
condaenv: /path/to/pipeline_env
```
or 

```
# ----------------------- #
# Visualisation pipeline DendrouLab
# ----------------------- #
# written by Charlotte Rich-Griffin and Fabiola Curion

# WARNING: Do not edit any line with the format `continuous_vars: &continuous_vars` or `continuous_vars: *continuous_vars`

# ------------------------
# compute resource options
# ------------------------
resources:
  # Number of threads used for parallel jobs
  # this must be enough memory to load your mudata and do computationally intensive tasks
  threads_high: 1
  # this must be enough memory to load your mudata and do computationally light tasks
  threads_medium: 1
  # this must be enough memory to load text files and do plotting, requires much less memory than the other two
  threads_low: 1

# path to conda env, leave blank if running native or your cluster automatically inherits the login node environment
condaenv: /path/to/env
```

