This folder contains files with the expression data for the angiosperm (DevSeq) and mammalian [(Brawand 2011)](https://pubmed.ncbi.nlm.nih.gov/22012392/) 1-1 ortholog genes. The columns are semicolon-separated and contain headers.

### Sample description

| <sub> File name  </sub>                                             | <sub> Data set </sub>  | <sub>Query species</sub>|<sub> Normalization </sub> | <sub> Metric </sub> |
| :------------------------------------------------------------------ | :----------------------------------------| :---- | :------------------------- | :-------------- |
| <sub> AT_brass_inter_count_mat_vsd_sample_names.csv </sub>          |<sub> DevSeq Brassicaceae </sub>          |<sub>AT<sub>|<sub>DESeq inter-organ</sub>|<sub>VST counts</sub>| 
| <sub> AT_brass_inter_tpm_mat_deseq_sample_names.csv </sub>          |<sub> DevSeq Brassicaceae </sub>          |<sub>AT<sub>|<sub>DESeq inter-organ</sub>|<sub>TPM</sub>| 
| <sub> AT_brass_intra_count_mat_vsd_sample_names.csv </sub>          |<sub> DevSeq Brassicaceae </sub>          |<sub>AT<sub>|<sub>DESeq intra-organ</sub>|<sub>VST counts</sub>| 
| <sub> AT_brass_intra_tpm_mat_deseq_sample_names.csv </sub>          |<sub> DevSeq Brassicaceae </sub>          |<sub>AT<sub>|<sub>DESeq intra-organ</sub>|<sub>TPM</sub>| 
| <sub> AT_core_inter_count_mat_vsd_sample_names.csv  </sub>          |<sub> DevSeq Angiosperm </sub>            |<sub>AT<sub>|<sub>DESeq inter-organ</sub>|<sub>VST counts</sub>| 
| <sub> AT_core_inter_tpm_mat_deseq_sample_names.csv </sub>           |<sub> DevSeq Angiosperm </sub>            |<sub>AT<sub>|<sub>DESeq inter-organ</sub>|<sub>TPM</sub>| 
| <sub> AT_core_intra_count_mat_vsd_sample_names.csv </sub>           |<sub> DevSeq Angiosperm </sub>            |<sub>AT<sub>|<sub>DESeq intra-organ</sub>|<sub>VST counts</sub>| 
| <sub> AT_core_intra_tpm_mat_deseq_sample_names.csv </sub>           |<sub> DevSeq Angiosperm </sub>            |<sub>AT<sub>|<sub>DESeq intra-organ</sub>|<sub>TPM</sub>| 
| <sub> AL_brass_inter_count_mat_vsd_sample_names.csv </sub>          |<sub> DevSeq Brassicaceae </sub>          |<sub>AL<sub>|<sub>DESeq inter-organ</sub>|<sub>VST counts</sub>| 
| <sub> AL_brass_inter_tpm_mat_deseq_sample_names.csv </sub>          |<sub> DevSeq Brassicaceae </sub>          |<sub>AL<sub>|<sub>DESeq inter-organ</sub>|<sub>TPM</sub>| 
| <sub> AL_brass_intra_count_mat_vsd_sample_names.csv </sub>          |<sub> DevSeq Brassicaceae </sub>          |<sub>AL<sub>|<sub>DESeq intra-organ</sub>|<sub>VST counts</sub>| 
| <sub> AL_brass_intra_tpm_mat_deseq_sample_names.csv </sub>          |<sub> DevSeq Brassicaceae </sub>          |<sub>AL<sub>|<sub>DESeq intra-organ</sub>|<sub>TPM</sub>| 
| <sub> AL_core_inter_count_mat_vsd_sample_names.csv  </sub>          |<sub> DevSeq Angiosperm </sub>            |<sub>AL<sub>|<sub>DESeq inter-organ</sub>|<sub>VST counts</sub>| 
| <sub> AL_core_inter_tpm_mat_deseq_sample_names.csv </sub>           |<sub> DevSeq Angiosperm </sub>            |<sub>AL<sub>|<sub>DESeq inter-organ</sub>|<sub>TPM</sub>| 
| <sub> AL_core_intra_count_mat_vsd_sample_names.csv </sub>           |<sub> DevSeq Angiosperm </sub>            |<sub>AL<sub>|<sub>DESeq intra-organ</sub>|<sub>VST counts</sub>| 
| <sub> AT_core_intra_tpm_mat_deseq_sample_names.csv </sub>           |<sub> DevSeq Angiosperm </sub>            |<sub>AL<sub>|<sub>DESeq intra-organ</sub>|<sub>TPM</sub>| 
|<sub>Brawand_inter_tpm_mat_deseq_sample_names_0_5_threshold.csv</sub>|<sub>Brawand mammalian (re-analyzed)</sub>|<sub>hsa<sub>|<sub>DESeq inter-organ</sub>|<sub>TPM</sub>| 
|<sub>Brawand_intra_tpm_mat_deseq_sample_names_0_5_threshold.csv</sub>|<sub>Brawand mammalian (re-analyzed)</sub>|<sub>hsa<sub>|<sub>DESeq intra-organ</sub>|<sub>TPM</sub>| 
|<sub>Brawand_inter_count_mat_vsd_sample_names_0_5_threshold.csv</sub>|<sub>Brawand mammalian (re-analyzed)</sub>|<sub>hsa<sub>|<sub>DESeq inter-organ</sub>|<sub>VST counts</sub>| 
|<sub>Brawand_intra_count_mat_vsd_sample_names_0_5_threshold.csv</sub>|<sub>Brawand mammalian (re-analyzed)</sub>|<sub>hsa<sub>|<sub>DESeq intra-organ</sub>|<sub>VST counts</sub>| 

<br/>

`inter-organ` indicates that samples were normalized between organs and species, whereas `intra-organ` indicates that samples were not normalized between organs.

Query species: `AT` = *Arabidopsis thaliana*; `AL` = *Arabidopsis lyrata*; `hsa` = *Homo sapiens* 

`Brawand_Supplementary_Data1` - this folder contains the processed RNA-Seq data published by Brawand et al.(2011). Expression values are RPKM and the files are tab-separated. Check the README.txt file within the folder for more details.
