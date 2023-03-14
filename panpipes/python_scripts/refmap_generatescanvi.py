
if args.generate_scanvi:
    if args.adata_reference is not None:
        reference_data = os.path.basename(args.adata_reference)
        mdata=read_anndata(args.adata_reference, use_muon=True, modality="all")
        if reference_data.endswith(".h5mu"):
            if "rna" not in mdata.mod.keys():
                sys.exit("we only support querying using RNA but your mdata doesn't contain rna")
            else:
                adata_ref = mdata["rna"].copy()
        else:
            adata_ref = mdata.copy()

        adata_ref.obs.loc[:, 'is_reference'] = 'Reference'
        adata_query.obs.loc[:, 'is_reference'] = 'Query'
        adata_query = adata_query[:, adata_ref.var_names].copy()
     
        #just to be clear, what i don't understand is why i can't split a scvi training from a scanvi update by saving the model and realoading afterwards. 
        #if i want to create a scanvi model from a scvi saved model there is no way to do it unless i retrain from scratch
        L.info(""" Generating scanvi reference from initialized SCVI, 
                    we'll train a scvi model first and
                    we expect the cell type labels into a cell_type column
                    """)   
        scvi.model.SCVI.setup_anndata(adata_ref, batch_key="tech", layer="counts")
        
        arches_params = dict(
            use_layer_norm="both",
            use_batch_norm="none",
            encode_covariates=True,
            dropout_rate=0.2,
            n_layers=2,
        )

        vae_ref = scvi.model.SCVI(
            adata_ref,
            **arches_params
        )
        vae_ref.train()    
        os.makedirs("models")
        dir_path_scan=os.path.join("models", "scvi_model")
        vae_ref.save(dir_path_scan, anndata=False)

        vae_ref_scan = scvi.model.SCANVI.from_scvi_model(
            vae_ref, #this is scvi reference model
            unlabeled_category="Unknown", #how to call the cells with no label in the "cell_type" column
            labels_key="cell_type")
        vae_ref_scan.train(max_epochs=20, n_samples_per_label=100) #to do add training args options
        
        dir_path_scan=os.path.join("models", "scanvi_model")
        vae_ref_scan.save(dir_path_scan, anndata=False)
        
        vae_q = scvi.model.SCANVI.load_query_data(
        adata_query,
        dir_path_scan
        )

        vae_q._unlabeled_indices = np.arange(adata_query.n_obs)
        vae_q._labeled_indices = []
        latent_choice= "X_scANVI"
        vae_q.train(**train_kwargs)
        adata_query.obsm["X_scANVI"] = vae_q.get_latent_representation()
        adata_query.obs["predictions"] = vae_q.predict()
        adata_full = adata_query.concatenate(adata_ref)
    else:
        sys.exit("To generate a scanvi reference i need a reference dataset to start from")
        

L.info("Done")

