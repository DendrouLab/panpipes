
## Rescomp specific advice
Created: 2020-09-21 CRG
Last edited: 2023-03-29 (updated for slurm)


### modules to use
These will need to be loaded every time you want to run the pipeline software
```
ARCH=`/apps/misc/utils/bin/get-cpu-software-architecture.py`
echo $ARCH

module purge
module use -a /apps/eb/dev/$ARCH/modules/all

module load Python/3.9.5-GCCcore-10.3.0
module load R-bundle-Bioconductor/3.13-foss-2021a-R-4.1.0
```




#### editing the .cgat.yml

Create a yml file for the cgat core pipeline software to read using:

```
vim ~/.cgat.yml
```
add the following information (pressing `i` to start inserting text)
```
cluster:
    queue_manager: slurm 
    queue: short
    options: --constraint=skl 

```
Currently I restrict jobs to either a nodes or e nodes (short.qc@@short.qe) because of the multiple architectures on rescomp.
If you want the code to be able to run on any of the nodes on rescomp then, install everything on rescomp3, specify `queue: short` in .cgat.yml and run the pipeline from rescomp3. (I think this works).
This works because ivybridge is the older hardware, so any softare compiled on ivybridge should work in skylake, but not the other way round. 

(to exit vim press ESC then type `:x`)

#### DRMAA_LIBRARY

if `echo $DRMAA_LIBRARY_PATH` does not return anything, add DRMAA_LIBRARY_PATH to your bash environment (the line below is specific to the rescomp server)

```
vim ~/.bashrc
```

```
export SLURM_CONF=/run/slurm/conf/slurm.conf  
export DRMAA_LIBRARY_PATH=/usr/lib64/libdrmaa.so
export SBATCH_ACCOUNT=*project*.prj
```

#### R packages
If you intend to use both architectures on rescomp then R packages need to be installed twice, once on rescomp1 and once on rescomp3. If you restrict jobs submission to just a or e nodes as described above, you only need to install on rescomp1.  

Details on how to use R on the server:
https://www.medsci.ox.ac.uk/divisional-services/support-services-1/bmrc/r-and-rstudio-on-the-bmrc-cluster


### Running jobs in the background
To run jobs in the background is is recommended to use a tmux session https://tmuxcheatsheet.com
to prevent the jobs from hanging up when you log off the server. It's a bit of a learning curve to work out tmux, but trust me it is worth it!
Open a tmux session, activate your modules and virtual environment


Alternatively use:
```
nohup sc_pipelines clustering make full &
```

### A note on GPUS

GPU usage is supported in `panpipes integration` and in `panpipes refmap`.
GPUs are utilised by scvi-tools and mofa.

GPU resources on bmrc are included here: https://www.medsci.ox.ac.uk/for-staff/resources/bmrc/gpu-resources-slurm

There is a section in the integration `pipeline.yml` for the specific configuration of the gpu queue:
```
queues:
  long: long
  gpu:  "gpu_short --gres=gres:gpu:1"
```


