Clustering 
===========


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