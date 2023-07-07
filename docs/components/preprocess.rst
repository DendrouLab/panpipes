Preprocessing
=============



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