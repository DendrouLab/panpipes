<style>
  .parameter {
    border-top: 4px solid lightblue;
    background-color: rgba(173, 216, 230, 0.2);
    padding: 4px;
    display: inline-block;
    font-weight: bold;
  }
</style>
# Integration methods implemented in panpipes

The panpipes integration workflow implements a variety of tools to batch correct individual modalities and/or integrate across modalities to produce a reduced dimension representation of the dataset.<br>
There are different tools available for each modality:  RNA (also referred to as GEX), PROT (can be referred to as ADT) and ATAC which can be integrated into any preferred combination customising the integration workflow configuration file and running `panpipes integration make full`. After the results of the integration are inspected, the final object is created with `panpipes integration make merge_integration`.<br> 

The ideal way to run `panpipes integration` is to use the output `MuData`file from `panpipes preprocess` since it will already be in the required format. 
However, if using independent MuData the object should contain normalised data in the X slot of each modality, a ‘raw_counts’ layer in each modality, and a sample_id column in each slot of the obs and the outer obs. 

Users can choose which integration method they want to apply based on their experiment, their experience with the tools or available benchmarks: we link all the relevant resources below. 
We don't believe in "one method fits all", we instead offer the user a framework to run multiple tools efficiently, keeping a record of the hyperparameters in the configuration files, so you can safely re-run your analysis and share it with collaborators. We will continue to update the integration methods offered in `panpipes` and we invite you to contribute yours!
The following table describes the different methods currently supported and their specificities: 

| Method    | type of integration         | modalities      | code                                                                              | references                                                                                           | benchmarks paper                                                                                           |
|-----------|-----------------------------|-----------------|-----------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------------|
| harmony   | unimodal (batch correction) | rna, atac, prot | [harmony](https://github.com/immunogenomics/harmony)                              | [Korsunsky  et al 2019](https://www.nature.com/articles/s41592-019-0619-0)                           | [Luecken et al 2022](https://www.nature.com/articles/s41592-021-01336-8)                                   |
| BBKNN     | unimodal (batch correction) | rna, atac, prot | [BBKNN](https://github.com/Teichlab/bbknn)                                        | [Polański et al 2020](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9883685/)                         | [Luecken et al 2022](https://www.nature.com/articles/s41592-021-01336-8)                                   |
| Scanorama | unimodal (batch correction) | rna             | [Scanorama](https://github.com/brianhie/scanorama)                                | [Hie, Bryson, and  Berger 2019](https://pubmed.ncbi.nlm.nih.gov/31061482/)                           | [Luecken et al 2022](https://www.nature.com/articles/s41592-021-01336-8)                                   |
| scVI      | unimodal (batch correction) | rna             | [scVI](https://github.com/scverse/scvi-tools)                                     | [Gayoso et al 2022](https://www.nature.com/articles/s41587-021-01206-w)                              | [Luecken et al 2022](https://www.nature.com/articles/s41592-021-01336-8)                                   |
| MultiVI   | multimodal                  | atac, rna       | [MultiVI](https://github.com/scverse/scvi-tools)                                  | [Ashuach et al 2023](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10406609/)                         | [Lee, Kaestner,  and Li 2023](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-023-03073-x) |
| totalVI   | multimodal                  | prot, rna       | [totalVI](https://github.com/scverse/scvi-tools)                                  | [Gayoso  et al. 2021](https://pubmed.ncbi.nlm.nih.gov/33589839/)                                     | [Makrodimitris et al 2024](https://academic.oup.com/bib/article/25/1/bbad416/7450271)                      |
| MOFA      | multimodal                  | rna, atac, prot | [MOFA](https://github.com/bioFAM/mofapy2)                                        | [Argelaguet et al 2020](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02015-1) | [Lee, Kaestner,  and Li 2023](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-023-03073-x) |
| WNN       | multimodal                  | rna, atac, prot | [WNN](https://muon.readthedocs.io/en/latest/api/generated/muon.pp.neighbors.html) | [Hao et al 2021](https://pubmed.ncbi.nlm.nih.gov/34062119/)                                          | [Lee, Kaestner,  and Li 2023](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-023-03073-x) |


