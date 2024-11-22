<table>
  <tr>
    <th colspan="3">Task ingest</th>
  </tr>
  <tr>
    <th>threads_high</th>
    <th>threads_medium</th>
    <th>threads_low</th>
  </tr>
  <tr>
    <td>Creating h5mu from filtered data files</td>
    <td>load_mudatas</td>
    <td>run_repertoire_qc</td>
  </tr>
  <tr>
    <td>Creating h5mu from bg data files</td>
    <td>load_bg_mudatas</td>
    <td>run_atac_qc</td>
  </tr>
  <tr>
    <td>rna QC</td>
    <td>downsample_bg_mudatas</td>
    <td>plot_qc</td>
  </tr>
  <tr>
    <td>prot QC</td>
    <td>run_scrublet</td>
    <td>10X metrics plotting</td>
  </tr>
  <tr>
    <td>prot QC</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <th colspan="3">Task preprocessed</th>
  </tr>
  <tr>
    <th>threads_high</th>
    <th></th>
    <th></th>
  </tr>
  <tr>
    <td>assess background</td>
    <td></td>
    <td>threads_low</td>
  </tr>
  <tr>
    <td>rna_preprocess</td>
    <td></td>
    <td>filter_mudata</td>
  </tr>
  <tr>
    <td>prot_preprocess</td>
    <td></td>
    <td>postfilterplot</td>
  </tr>
  <tr>
    <td>atac_preprocess</td>
    <td></td>
    <td>downsample</td>
  </tr>
  <tr>
    <th colspan="3">Task integration</th>
  </tr>
  <tr>
    <th>threads_high</th>
    <th>threads_medium</th>
    <th>threads_low</th>
  </tr>
  <tr>
    <td>run_no_batch_correct_rna</td>
    <td>Evaluation</td>
    <td>run_lisi</td>
  </tr>
  <tr>
    <td>run_bbknn_rna</td>
    <td>plot_umaps</td>
    <td></td>
  </tr>
  <tr>
    <td>run_harmony_rna</td>
    <td>run_scib_metrics</td>
    <td></td>
  </tr>
  <tr>
    <td>run_combat_rna</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td>run_scanorama_rna</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td>run_scvi_rna</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td>run_no_batch_correct_prot</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td>run_harmony_prot</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td>run_bbknn_prot</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td>run_combat_prot</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td>run_no_batch_correct_atac</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td>run_harmony_atac</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td>run_bbknn_atac</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td>run_totalvi</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td>run_multivi</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td>run_mofa</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td>run_wnn</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td>merge_integration</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <th colspan="3">Task clustering</th>
  </tr>
  <tr>
    <th>threads_high</th>
    <th>threads_medium</th>
    <th>threads_low</th>
  </tr>
  <tr>
    <td>run_neighbors</td>
    <td>run_clustering</td>
    <td>clustering</td>
  </tr>
  <tr>
    <td>run_umap</td>
    <td>collate_mdata</td>
    <td>clustering</td>
  </tr>
  <tr>
    <td>find_markers</td>
    <td>plot_cluster_umaps</td>
    <td></td>
  </tr>
  <tr>
    <td></td>
    <td>plot_markers</td>
    <td></td>
  </tr>
  <tr>
    <th colspan="3">Task vis</th>
  </tr>
  <tr>
    <th>threads_high</th>
    <th></th>
    <th>threads_low</th>
  </tr>
  <tr>
    <td>plot_custom_markers_per_group</td>
    <td></td>
    <td>plot_metrics</td>
  </tr>
  <tr>
    <td>plot_custom_markers_umap</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td>plot_categorical_umaps</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td>write_obs</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td>plot_scatters</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td>Task refmap</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td>threads_high</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td>run_refmap_scvi</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td>run_scib_refmap</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <th colspan="3">Task preprocessed spatial</th>
  </tr>
  <tr>
    <th>threads_high</th>
    <th></th>
    <th>threads_low</th>
  </tr>
  <tr>
    <td>spatial_preprocess</td>
    <td></td>
    <td>filter_mudata</td>
  </tr>
  <tr>
    <th colspan="3">Task Spatial</th>
  </tr>
  <tr>
    <th>threads_high</th>
    <th></th>
    <th>threads_low</th>
  </tr>
  <tr>
    <td>load_mudata</td>
    <td></td>
    <td>plotQC_spatial</td>
  </tr>
</table>