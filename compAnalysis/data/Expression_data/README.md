This folder contains files with the expression data for the angiosperm (DevSeq data set) and mammalian [(Brawand data set)](https://pubmed.ncbi.nlm.nih.gov/22012392/) 1-1 ortholog genes. The columns are semicolon-separated and contain headers.

### Detailed description

| File name        | Data set           | Normalization  |   Metric  |
| ------------- |:-------------:| -----:| -----:|
| comp_Brass_AT_inter_tpm_mat_deseq_sample_names.csv      | DevSeq angiosperm | DESeq2 inter-organ | TPM  | 
| Brawand_inter_tpm_mat_deseq_sample_names_0_5_threshold.csv      | Brawand mammalian (re-analyzed) | DESeq2 inter-organ | TPM  | 
| Brawand_intra_tpm_mat_deseq_sample_names_0_5_threshold.csv      | Brawand mammalian (re-analyzed) | DESeq2 intra-organ | TPM  | 
| Brawand_inter_count_mat_vsd_sample_names_0_5_threshold.csv      | Brawand mammalian (re-analyzed) | DESeq2 inter-organ | VST counts  | 
| Brawand_intra_count_mat_vsd_sample_names_0_5_threshold.csv      | Brawand mammalian (re-analyzed) | DESeq2 intra-organ | VST counts  | 



`inter-organ` indicates that samples were normalized between organs and species, whereas `intra-organ` indicates that samples were not normalized between organs.

`Brawand_Supplementary_Data1` - this folder contains the processed RNA-Seq data published by Brawand et al.(2011). Expression values are RPKM and the files are tab-separated. Check the README.txt file within the folder for more details.
