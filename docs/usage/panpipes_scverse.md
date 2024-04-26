# Panpipes and the scverse


Panpipes is implemented in python to rely on the `scverse`, an ecosystem of tools for single cell 'omics analysis in python.
Checkout the [scverse](https://scverse.org/) documentation and how to contribute your package! 



Panpipes has at its core [`AnnData`](https://anndata.readthedocs.io/en/latest/) and [`MuData`](https://mudata.readthedocs.io/en/latest/), for handling annotated data matrices in memory and on disk, with the best of the pandas and xarray functionalities.

<figure>
    <img src="https://github.com/DendrouLab/panpipes/blob/main/docs/img/anndata_schema.svg?raw=true" alt="img1" width="40%">
    <figcaption>`AnnData` is anndata is a container for handling annotated data matrices objects.</figcaption>
</figure>




<figure>
    <img src="https://github.com/DendrouLab/panpipes/blob/main/docs/img/mudata_paper.svg?raw=true" alt="img2" width="40%">
    <figcaption>`MuData` is a dictionary of `AnnData` objects.</figcaption>
</figure>


Please use these links to familiarize with these data structures


