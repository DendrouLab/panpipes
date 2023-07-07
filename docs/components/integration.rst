Integration
============


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