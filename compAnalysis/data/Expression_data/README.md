This folder contains files with the expression data for the angiosperm (DevSeq data set) and mammalian [(Brawand data set)](https://pubmed.ncbi.nlm.nih.gov/22012392/) 1-1 ortholog genes. The columns are semicolon-separated and contain headers. The first column contains the identifiers of the orthologous genes. 

### Detailed description

`Brawand_inter_tpm_mat_deseq_sample_names_0_5_threshold.csv` contains the TPM expression values normalized by DESeq2. Expression data of the 1-1 ortholog genes was normalized between organs and species. Data was re-analyzed from the original raw sequencing data published in Brawand et al.(2011).

`Brawand_intra_tpm_mat_deseq_sample_names_0_5_threshold.csv` Same as above, but the data was normalized within organs between species (intra-organ normalization).

`Brawand_inter_count_mat_vsd_sample_names_0_5_threshold.csv` contains the count data normalized by DESeq2 using the variance stabilizing transformation (VST). Expression data of the 1-1 ortholog genes was normalized between organs and species. Data was re-analyzed from the original raw sequencing data published in Brawand et al.(2011).

`Brawand_intra_count_mat_vsd_sample_names_0_5_threshold.csv` Same as above, but the data was normalized within organs between species (intra-organ normalization).

`Brawand_Supplementary_Data1` - this folder contains the processed RNA-Seq data published by Brawand et al.(2011). Expression values are RPKM and the files are tab-separated. Check the README.txt file within the folder for more details.
