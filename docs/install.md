
# Installation of panpipes

We recommend running panpipes within a virtual environment to prevent conflicts.
In the following, we provide instructions on how to do this using conda, mamba, or python venv.

> **Note**: For installation instructions on **Apple machines with M chips**, scroll down.
> If you are working on a **Windows** machine, please use Windows Subsystem for Linux (WSL) to run panpipes in order to have access to the latest r-base version.
> **Oxford BMRC Rescomp** users find additional advice on the installation [here](https://github.com/DendrouLab/panpipes/blob/main/docs/installation_rescomp.md).

## Option 1: Installation in manually configured conda environment
To run panpipes, we install it in a conda environment with R and python.
Panpipes has a lot of dependencies, so you may want to consider the faster [`mamba`](https://mamba.readthedocs.io/en/latest/index.html) instead of `conda` for installation.
Panpipes can be installed via different methods, either from `PyPi` or from the Github repository.
Option 1 describes the installation via `PyPi` in a manually configured conda environment (no cloning of the repository necessary).

### Create conda environment
First, create a conda environment to run panpipes in. For that, we replicate the suggestions made [here](https://www.biostars.org/p/498049/):

```bash 
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

### Install panpipes
Finally, we install panpipes.
You can install it either from `PyPi` (shown here) or a nightly version from the Github repository (shown in Option 2):

```bash
# Install panpipes from PyPi
pip install panpipes
```

If you intend to use panpipes for spatial analysis, instead install:

```bash
pip install 'panpipes[spatial]'
```
The extra `[spatial]` includes the `squidpy`, `cell2location`, and `tangram-sc` packages.

The packages we support are in active development. We provide installation instructions to deal with the issues encountered, which we test in our [github actions](https://github.com/DendrouLab/panpipes)
- To install a version of panpipes that has a working MultiVI please install panpipes with `pip install 'panpipes[multivipatch]'`
- To install a version of panpipes to run old scvi-tools models in the refmap workflow, please install panpipes with `pip install 'panpipes[refmap_old]'`


## Option 2: Install nightly panpipes version with preconfigured conda config file

If you prefer to use the most recent development version, install the nightly panpipes version from the Github repository.
To make the installation easier, we provide a minimal conda config file in `pipeline_env.yaml`.

### Clone the repository
First, clone the [panpipes repository](https://github.com/DendrouLab/panpipes) and navigate to the root directory of the repository:

```bash
git clone https://github.com/DendrouLab/panpipes.git
cd panpipes
```

### Create conda environment and install nightly panpipes version
Then, create the conda environment and install the nightly version of panpipes using the following command:

```bash
conda env create --file=pipeline_env.yaml 
conda activate pipeline_env
pip install -e .
```
To install the spatial dependencies from the repo, use

```bash
pip install .[spatial]
```

The packages we support are in active development. We provide installation instructions to deal with the issues encountered, which we test in our [github actions](https://github.com/DendrouLab/panpipes)
- To install a version of panpipes that has a working MultiVI please install panpipes with `pip install '.[multivipatch]'`
- To install a version of panpipes to run old scvi-tools models in the refmap workflow, please install panpipes with `pip install '.[refmap_old]'`


Panpipes requires the unix package `time`. 
You can check if it installed with `dpkg-query -W time`.
If `time` is not already installed, you can install it using:

```bash
conda install time
```

or

```bash
apt-get install time
```


## Option 3: python venv environment
As an alternative to a conda environment, you can also install panpipes in a python virtual environment.
Navigate to where you want to create your virtual environment and follow the steps below to create a `pip` virtual environment.

```bash
# Create a panpipes/venv folder
python3 -m venv --prompt=panpipes python3-venv-panpipes/
```

Activate the environment:

```bash
source python3-venv-panpipes/bin/activate
```

As explained above, you can install panpipes from `PyPi` with:

```bash
pip install panpipes
```

Alternatively, you can install a nightly version of panpipes by cloning the Github repository (see instructions above).

### R packages installation in python venv

If you are using a venv virtual environment, the pipeline will call a local R installation, so make sure R is installed and install the required packages with the command we provide below.
(This executable requires that you specify a CRAN mirror in your `.Rprofile`).
for example, add this line to your `.Rprofile` to automatically fetch the preferred mirror:

> **Note:** Remember to customise with your preferred [R mirror](https://cran.r-project.org/mirrors.html).

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

### Check installation

To check the installation was successful, run the following line:

```bash
panpipes --help
```

A list of available pipelines should appear!


## Installation on Apple Silicon M chips
If you intend to install panpipes via conda on a macOS machine with M-Chip, you might face issues when installing or using certain workflows of panpipes.
This is because panpipes relies on `scvi-tools`, which currently [only supports execution on Apple Silicon machines when installed using a native Python version](https://docs.scvi-tools.org/en/stable/installation.html#apple-silicon) (due to a dependency on JAX).

Follow these steps to install panpipes on an Apple Silicon machine:

1. Install [Homebrew](https://brew.sh/)

2. Install Apple Silicon version of Mambaforge (If you already have Anaconda/Miniconda installed, make sure
   having both mamba and conda won't cause conflicts). Additionally, we need clang which is included in llvm, so we install that as well:

```bash
brew install --cask mambaforge
brew install llvm
```

3. Create a new environment using mamba (here with python 3.10) and activate it:

```bash
conda config --add channels conda-forge
conda config --set channel_priority strict
# you should remove the strict priority afterwards!
mamba search r-base
mamba create --name pipeline_env
mamba activate pipeline_env
```

4. Add the osx-64 channel to the environment, then install Python and R
Because not all R packages are available via the ARM64 channel, we need to specify the osx-64 channel to install all required R packages:

```bash
conda config --env --set subdir osx-64
mamba install python=3.10 r-base=4.3.0
```

5. Install R dependencies and panpipes itself:

```bash
mamba install -c conda-forge r-tidyverse r-optparse r-ggforce r-ggraph r-xtable r-hdf5r r-clustree r-cowplot
pip install panpipes # or install the nightly version by cloning the repository as described above
```

6. Ensure the `time` package is installed
If the `time` package is not already installed, you can install it using:

```bash
mamba install time
```
You're all set to run `panpipes` on your local machine.

>Note: If you encounter an issue with the `jax` dependency when running panpipes, try reinstalling it from source:
>```bash
>pip uninstall jax jaxlib
>mamba install -c conda-forge jaxlib
>mamba install -c conda-forge jax
>```

If you want to configure it on an HPC server, follow the instructions in the following section.

## Pipeline configuration for HPC clusters

This section is for users who want to run panpipes on a High-Performance Computing (HPC) cluster with a job scheduler like SGE or SLURM.

>Note: You only need this configuration step if you want to use an HPC to dispatch individual task as separate parallel jobs. You won't need this for a local installation of panpipes.

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
