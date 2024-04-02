
# Installation of panpipes

>Note: Oxford BMRC Rescomp users find additional advice on the installation [here](https://github.com/DendrouLab/panpipes/blob/main/docs/installation_rescomp.md).

## Create virtual environment and install panpipes

We recommend running panpipes within a virtual environment to prevent conflicts.

### Option 1: Installation in conda environment (Recommended)
>Note: For installation instructions on Apple machines with M chips, scroll down.

To run panpipes, we install it in a conda environment with R and python.
Panpipes has a lot of dependencies, so you may want to consider the faster [`mamba`](https://mamba.readthedocs.io/en/latest/index.html) instead of `conda` for installation.
Panpipes can be installed via different methods, either from PyPi or from the Github repository.
We recommend using Option 1.1 to install the nightly version of panpipes.

#### Option 1.1: Nightly panpipes version with preconfigured conda config file (Recommended)
We recommend installing a nightly version of panpipes.
For that, we provide a minimal conda config file in `pipeline_env.yaml`.
First, clone this repository and navigate to the root directory of the repository:

```
git clone https://github.com/DendrouLab/panpipes.git
cd panpipes
```

Then, create the conda environment and install the nightly version of panpipes using the following command:

```
conda env create --file=pipeline_env.yaml 
conda activate pipeline_env
pip install -e .
```

#### Option 1.2: Manual conda environment creation
As an alternative to the preconfigured conda environment, you can create a conda environment manually.

```bash
#This follows the suggestions made here: [https://www.biostars.org/p/498049/](https://www.biostars.org/p/498049/) 
conda config --add channels conda-forge
conda config --set channel_priority strict
# you should remove the strict priority afterwards!
conda search r-base
conda create --name pipeline_env python=3.10 r-base=4.3.0
```

Next, we activate the environment:

```bash
conda activate pipeline_env
```

Let's first install the R packages:

```bash
conda install -c conda-forge r-tidyverse r-optparse r-ggforce r-ggraph r-xtable r-hdf5r r-clustree r-cowplot
```

Finally, we install panpipes, which you can install either from PyPi or from the Github repository.

##### Installing panpipes from PyPi

You can install `panpipes` directly from `PyPi` with:

```bash
pip install panpipes
```

If you intend to use panpipes for spatial analysis, instead install:

```bash
pip install 'panpipes[spatial]'
```
The extra `[spatial]` includes squidpy, cell2location, and tangram-sc packages.

##### Nightly version of panpipes

If you prefer to use the most recent dev version, install panpipes from Github:

```bash
git clone https://github.com/DendrouLab/panpipes
cd panpipes
pip install -e .
```

Panpipes requires the unix package `time`. 
You can check if it installed with `dpkg-query -W time`. If time not already installed, you can 

```bash
conda install time
```

or

```bash
apt-get install time
```

#### Installation on Apple Silicon M chips
If you intend to install panpipes via conda on a macOS machine with M-Chip, you might face issues when installing or using certain workflows of panpipes.
This is because panpipes relies on [scvi-tools], which currently only supports execution on Apple Silicon machines when installed using a native Python version (owing to a dependency on JAX).

Follow these steps to install pertpy on an Apple Silicon machine:

1. Install [Homebrew](https://brew.sh/)

2. Install Apple Silicon version of Mambaforge (If you already have Anaconda/Miniconda installed, make sure
   having both mamba and conda won't cause conflicts). Additionally, we need clang which is included in llvm, so we install that as well.

```bash
brew install --cask mambaforge
brew install llvm
```

3. Create a new environment using mamba (here with python 3.10) and activate it

```bash
conda config --add channels conda-forge
conda config --set channel_priority strict
# you should remove the strict priority afterwards!
mamba search r-base
mamba create --name pipeline_env
mamba activate pipeline_env
```

4. Add the osx-64 channel to the environment, then install Python and R
Because not all R packages are available via the ARM64 channel, we need to specify the osx-64 channel to install all required R packages.

```bash
conda config --env --set subdir osx-64
mamba install python=3.10 r-base=4.3.0
```

5. Install dependencies

```bash
conda install -c conda-forge r-tidyverse r-optparse r-ggforce r-ggraph r-xtable r-hdf5r r-clustree r-cowplot
pip install panpipes
```


### Option 2: python venv environment

Navigate to where you want to create your virtual environment and follow the steps below to create a pip virtual environment.

```bash
python3 -m venv --prompt=panpipes python3-venv-panpipes/
# This will create a panpipes/venv folder
```

Activate the environment

```bash
source python3-venv-panpipes/bin/activate
```

As explained in the conda installation, you can install `panpipes` with:

```bash
pip install panpipes
```

or install a nightly version of panpipes by cloning the Github repository.

#### R packages installation in python venv

If you are using a venv virtual environment, the pipeline will call a local R installation, so make sure R is installed and install the required packages with the command we provide below.
(This executable requires that you specify a CRAN mirror in your `.Rprofile`).
for example, add this line to your `.Rprofile` to automatically fetch the preferred mirror:

*remember to customise with your preferred [R mirror](https://cran.r-project.org/mirrors.html).*

```R
  options(repos = c(CRAN="https://cran.uni-muenster.de/"))
```

Now, to automatically install the R dependecies, run:

 ```bash
panpipes install_r_dependencies
 ```

If you want more control over your installation use the [script on github](https://github.com/DendrouLab/panpipes/blob/main/panpipes/R_scripts/install_R_libs.R).
Running with the option `--vanilla` or `--no-site-file` prevents R from reading your `.Renvironment` or `.Rprofile` in case you want to use different settings from you local R installation.
You can expect the installation of R libraries to take quite some time, this is not something related to `panpipes` but how R manages their libraries and dependencies!

#### Check installation

To check the installation was successful run the following line

```bash
panpipes --help
```

A list of available pipelines should appear!

You're all set to run `panpipes` on your local machine.
If you want to configure it on a HPC server, follow the next instructions.

## Pipeline configuration for HPC clusters

(For SGE or SLURM clusters)
*Note: You only need this configuration step if you want to use an HPC to dispatch individual task as separate parallel jobs. You won't need this for a local installation of panpipes.*

Create a yml file for the cgat core pipeline software to read

```bash
vim ~/.cgat.yml
```

For SGE servers:

```bash
cluster:
  queue_manager: sge
  queue: short
condaenv:
```

For Slurm servers:

```bash
cluster:
    queue_manager: slurm
    queue: short
```

There are likely other options that relate to **your specific server**.
These are added as:

```yaml
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

```bash
PATH_TO_DRMAA = ""
echo "export DRMAA_LIBRARY_PATH=$PATH_TO/libdrmaa.so.1.0" >> ~/.bashrc
```

### Specifying Conda environments to run panpipes

If using Conda environments, you can use one single big environment (the instructions provided do just that) or create one for each of the workflows in panpipes, (i.e. one workflow = one environment).
The environment (s) should be specified in the .cgat.yml global configuration file or in each of the single workflows pipeline.yml configuration files, and it will be picked up by the pipeline as the default environment.
Please note that if you specify the Conda environment in the workflows configuration file this will be the first choice to run the pipeline.

If no environment is specified, the default behaviour of the pipeline is to inherit environment variables from the node where the pipeline is run. However there have been reported issues on SLURM clusters where this was not the default behaviour.
In such instances we recommend adding the Conda environment param in the .cgat.yml file or in each of the pipeline.yml independently.

```yaml
cluster:
    queue_manager: slurm
    queue: cpu_p
    options: --qos=xxx --exclude=compute-node-0[02-05,08-19],compute-node-010-0[05,07,35,37,51,64,68-71]
condaenv: /path/to/pipeline_env
```

or:

```yaml
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
