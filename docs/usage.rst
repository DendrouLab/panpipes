Usage guidelines
==============================

See `installation instrcutions
here <https://github.com/DendrouLab/panpipes/blob/main/docs/install.md>`__

Oxford BMRC Rescomp users find additional advice in
`docs/installation_rescomp <https://github.com/DendrouLab/panpipes/blob/main/docs/installation_rescomp.md>`__

General principles for running pipelines
''''''''''''''''''''''''''''''''''

Run the pipeline from the login node on your server, it will use in
built the job submission system to submit jobs.

Navigate to the directory where you want to run your analysis (this
should not be within the panpipes folder, or your virutal environment
folder)

::

   mkdir data_dir/
   cd data_dir/
   panpipes qc_mm config

This will produce two files, ``pipeline.log`` and ``pipeline.yml``

Edit ``pipeline.yml`` as appropriate for your data, following the
instructions within the yml file.

Then check which jobs will run with the command

::

   panpipes qc_mm show full

The output of this will show a list of tasks that will be run as part of
the pipeline.

To run use the command

::

   panpipes qc_mm make full

Occasionally you might want to run tasks individually (e.g. to assess
outputs before deciding the parameters for the next step) In order to do
this you can run any task in the ``show full`` list such as:

::

   panpipes qc_mm make plot_tenx_metrics

Pipeline components
''''''''''''''''''''''''''''''''''

Run each of pipeline qc, integration and clustering in separate folders.
### QC

1. Generate sample submission file

   -  `more details on creating the submission
      file <https://github.com/DendrouLab/panpipes/blob/main/docs/setup_for_qc_mm.md>`__

2. Generate qc genelists

   -  `more
      details <https://github.com/DendrouLab/panpipes/blob/main/docs/gene_list_format.md>`__

3. For adt assay - generate the protein metadata file
   `example <(https://github.com/DendrouLab/panpipes/blob/main/resources/protein_metadata_w_iso.md)>`__.
   This file is integrated into the mdata[‘prot’].var slot.
4. Generate config file (``panpipes qc_mm config``)
5. Edit the pipeline.yml file for your dataset

   -  this is explained step by step within the pipeline.yml file

6. Run complete qc pipeline with ``panpipes qc_mm make full``
7. Use outputs to decide filtering thresholds.

   -  **Note that the actual filtering occurs in the first step of
      Preprocess pipeline**
   -  TODO: create doc to explain the pipeline outputs

The h5mu file outputted from ``qc_mm`` contains concatenated raw counts
from all samples in the submission file, plus qc metrics are computed,
and these qc metrics are visualised in a variety of plots to aid the
user to determine data quality and filtering thresholds.

Preprocess
~~~~~~~~~~

1. In a new folder, generate config file for integration,
   ``panpipes preprocess config``
2. edit the pipeline.yml file

   -  The filtering options are dynamic depending on your qc_mm inputs
      `more details
      here <https://github.com/DendrouLab/panpipes/blob/main/docs/filter_dict_instructions.md>`__
   -  There are lots of options for normalisation explained in the
      pipeline.yml

3. Run complete preprocess pipeline with
   ``panpipes preprocess make full``

The h5mu outputted from ``preprocess`` is filtered and normalised, and
for rna highly variable genes are computed.

Integration
~~~~~~~~~~~

1. In a new folder, generate config file for integration,
   ``panpipes integration config`` and edit the pipeline.yml file.
2. Run ``panpipes integration make plot_pcas`` and assess the post
   filtering qc plots, and pca outputs
3. Run batch correction with
   ``panpipes integration make batch_correction`` (or run steps 2 and 3
   in one go with ``panpipes integration make full``)
4. Use pipeline outputs to decide on the best batch correction method
5. Edit the integration pipeline yml with your preferred batch
   correction
6. Run ``panpipes integration make merge_batch_correction``

Refmap
~~~~~~

1. In a new folder, generate config file for integration,
   ``panpipes refmap config`` and edit the pipeline.yml file.
2. Run complete refmap pipeline with ``panpipes refmap make full``

Clustering
~~~~~~~~~~

1. In a new folder, generate config file for integration,
   ``panpipes clustering config`` and edit the pipeline.yml file.
2. Run the clustering pipeline
   ``panpipes clustering make cluster_analysis``. This will do the
   initial nearest neighbours and clustering for the parameters you
   specify.
3. Decide on the best values for k nearest neighbours based on UMAPs and
   clustree results. Once decided delete the folders for the parameters
   you don’t need and delete those from the pipeline.yml.
4. Find markers for each of your cluster resolutions with
   ``panpipes clustering make marker_analysis`` (Again you could run all
   the clustering pipeline at once with
   ``panpipes clustering make full`` but by making decisions along the
   way you’ll reduce the computation and file size burden of the
   pipeline)

Vis
~~~

1. In a new folder, generate config file for integration,
   ``panpipes vis config`` and edit the pipeline.yml file.
2. Prepare plotting gene list files

   -  `more
      details <https://github.com/DendrouLab/panpipes/blob/main/docs/gene_list_format.md>`__

3. Run complete refmap pipeline with ``panpipes vis make full``

To repeat the pipeline after editing the pipeline.yml, delete the files
in log and repeat step 3.

Running pipeline modules from different entry points.
''''''''''''''''''''''''''''''''''

`see
details <https://github.com/DendrouLab/panpipes/blob/main/docs/different_entry_points.md>`__
