These files provide the expression values (RPKM), normalized on the 1000 most rank conserved genes, for the amniotes and primates 1-1 ortholog datasets. 

NormalizedRPKM_ConstitutiveExons_Amniote1to1Orthologues.txt - this file provides the RPKM values for the 5,636 genes present in the amniote 1-1 orthologues set. The expression levels were computed on the constitutive exons, defined independently for each species. 

NormalizedRPKM_ConstitutiveExons_Primate1to1Orthologues.txt - this file provides the RPKM values for the 13,277 genes present in the primate 1-1 orthologues set. The expression levels were computed on the constitutive exons, defined independently for each species. 

NormalizedRPKM_ConstitutiveAlignedExons_Amniote1to1Orthologues.txt - this file provides the RPKM values for the 5,636 genes present in the amniote 1-1 orthologues set. The expression levels were computed on the "constitutive aligned" exons (i.e., those parts of the constitutive exons that were perfectly aligned - without gaps - for all the species in the dataset).

NormalizedRPKM_ConstitutiveAlignedExons_Primate1to1Orthologues.txt - this file provides the RPKM values for the 13,277 genes present in the primate 1-1 orthologues set. The expression levels were computed on the "constitutive aligned" exons (i.e., those parts of the constitutive exons that were perfectly aligned - without gaps - for all the primate species in the dataset).

Note that we excluded from these datasets the genes for which no constitutive exons could be defined, or for which the constitutive exons could not be aligned. 

Note also that for the analysis of lineage-specific expression changes we did not consider those genes whose exons overlapped with those of neighboring genes. This list of genes is provided in the second supplementary dataset, which provides also expression levels computed with TopHat read alignments, on Ensembl annotations. 

The columns are tab-separated, and the files contain headers. The first columns  (10 for the amniote dataset, 5 for the primate dataset) provide the Ensembl gene identifiers for the species included in that dataset. The other columns correspond to different samples.
