# BMRC cluster specific advice

Created: 2020-09-21 CRG
Last edited: 2023-08-18 (updated for slurm)

## modules to use

These will need to be loaded every time you want to run the pipeline software


```bash
module purge

module load Python/3.9.5-GCCcore-10.3.0
module load R-bundle-Bioconductor/3.13-foss-2021a-R-4.1.0
```

You can choose different versions of Python and R but the GCCcore that they use must be the same.

### Step 1 - create virutal environment

It is advisable to run everything in a virtual environment either pip or Conda.

Using pip venv:
Navigate to where you want to create your virtual environment and follow the steps below to create a pip virtual environment. 
On rescomp, it's advsiable to make a folder where you store these environments that is separate from your data analysis.

e.g. /well/{project}/users/{user}/devel

where {project} is your workspace name
and {user} is your bmrc ID.

```bash
project=""
user=""
mkdir /well/${project}/users/${user}/devel
cd /well/${project}/users/${user}/devel

# create the evironment.
python3 -m venv --prompt=panpipes python3-venv-panpipes/
# This will create a panpipes/venv folder
```

activate the environment

```bash
source python3-venv-panpipes/bin/activate
```

### Step 2 Download and install this repo

If you have not already set up SSH keys for github first follow these [instructions](https://github.com/DendrouLab/panpipes/docs/set_up_ssh_keys_for_github.md):


```bash
cd /well/${project}/users/${user}/devel
git clone https://github.com/DendrouLab/panpipes
cd panpipes
pip install .
```

OR if you anticipate doing development on panpipes:

```bash
pip install --editable .
```

The pipelines are now installed as a local python package.

### Step 3 Modifying PYTHONPATH 

On the BMRC, the path to packages for loaded python SciPy-bundle supercedes the python venv python. This affects the numpy package version being used when the venv is activated, with the module Numpy version being much older than the one panpipes depends on ("numpy>=1.22.4"). Not doing this might lead to numpy relatad errors while running the various panpipes workflows downstream.
So we need to make sure that after activating the python venv, the `PYTHONPATH` has the path to the venv python packages as the first path. to do this do the following:

```bash
path_to_virtual_env=/well/${project}/users/${user}/python3-venv-panpipes
source ${path_to_virtual_env}/bin/activate
export PYTHONPATH=${path_to_virtual_env}/lib/python3.9/site-packages:$PYTHONPATH
```

One can check if we have the correct Numpy version by doing `pip list` after setting the PYTHONPATH this way to be sure.

### Step 4 installing R requirements

The pipelines use R (mostly for ggplot visualisations).

If you are using a venv virtual environment,  the pipeline will call a local R installation, so make sure R is installed and install the required packages with the following command

From within the panpipes repo folder run:

```bash
cd /well/${project}/users/${user}/devel
cd panpipes
Rscript panpipes/R_scripts/install_R_libs.R
```

<!-- 
```
pip install git+https://github.com/DendrouLab/panpipes
``` -->

If you intend to use both architectures on rescomp then R packages need to be installed twice, once on rescomp1 and once on rescomp3. If you restrict jobs submission to just a or e nodes as described above, you only need to install on rescomp1.  

Details on how to use R on the server:
https://www.medsci.ox.ac.uk/divisional-services/support-services-1/bmrc/r-and-rstudio-on-the-bmrc-cluster

### Step 5 pipeline configuration (for SGE or SLURM)

Create a yml file for the cgat core pipeline software to read using:

```bash
vim ~/.cgat.yml
```

add the following information (pressing `i` to start inserting text)

```yaml
cluster:
    queue_manager: slurm 
    queue: short
    options: --constraint=skl 

```

Currently, I restrict jobs to either a nodes or e nodes (short.qc@@short.qe) because of the multiple architectures on rescomp.
If you want the code to be able to run on any of the nodes on rescomp then, install everything on rescomp3, specify `queue: short` in .cgat.yml and run the pipeline from rescomp3.
This works because ivybridge is the older hardware, so any software compiled on ivybridge should work in skylake, but not the other way round.

(to exit vim press ESC then type `:x`)


#### DRMAA_LIBRARY

if `echo $DRMAA_LIBRARY_PATH` does not return anything, add DRMAA_LIBRARY_PATH to your bash environment (the line below is specific to the rescomp server)

```bash
echo "export SLURM_CONF=/run/slurm/conf/slurm.conf"  >> ~/.bashrc
echo "export DRMAA_LIBRARY_PATH=/usr/lib64/libdrmaa.so"  >> ~/.bashrc
echo "export SBATCH_ACCOUNT=${project}.prj" >> ~/.bashrc

tail ~/.bashrc
```

Congratulations you have now installed panpipes

### Running jobs in the background

To run jobs in the background it is recommended to use a tmux session https://tmuxcheatsheet.com
to prevent the jobs from hanging up when you log off the server. It's a bit of a learning curve to work out tmux, but trust me it is worth it!
Open a tmux session, activate your modules and virtual environment

Alternatively use:

```bash
nohup sc_pipelines clustering make full &
```

### A note on GPUS

GPU usage is supported in `panpipes integration` and in `panpipes refmap`.
GPUs are utilised by scvi-tools and mofa.

GPU resources on bmrc are included here: https://www.medsci.ox.ac.uk/for-staff/resources/bmrc/gpu-resources-slurm

There is a section in the integration `pipeline.yml` for the specific configuration of the gpu queue:

```yaml
queues:
  long: long
  gpu:  "gpu_short --gres=gres:gpu:1"
```
