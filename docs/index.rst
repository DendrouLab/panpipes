 Panpipes - multimodal single cell pipelines 
==================================================


**Panpipes** Created and Maintained by Charlotte Rich-Griffin and Fabiola Curion  
Additional contributors: Devika Agarwal and Tom Thomas 

See our [preprint](https://www.biorxiv.org/content/10.1101/2023.03.11.532085v1):  
Panpipes: a pipeline for multiomic single-cell data analysis  
Charlotte Rich-Griffin, Fabiola Curion, Tom Thomas, Devika Agarwal, Fabian J. Theis, Calliope A. Dendrou.  
bioRxiv 2023.03.11.532085;  
doi: https://doi.org/10.1101/2023.03.11.532085


# Introduction
These pipelines use cgat-core pipeline software

Available pipelines:
1. "qc_mm" : for the ingestion of data and computation of QC metrics' 
2. "preprocess" : for filtering and normalising of each modality
3. "integration" : integrate and batch correction using  single and multimodal methods
4. "clustering" : cell clustering on single modalities
5. "refmap" : transfer scvi-tools models from published data to your data
6. "vis" : visualise metrics from other pipelines in context of experiment metadata

Check out the :doc:`usage` section for further information, including
how to :ref:`installation` the project.
