Integration of gene expression data into the genome-scale metabolic model of L. kluyveri

1. Download the series matrix format of gene expression data from GEO Omnibus-GSE48135

2. The following R code was used to check for the normality of the gene expression data:

Code: NormalCheckExpressionDataset.R 
Inputs: GSE48135.txt, UniversalIDmap.csv
Outputs: GeneExpressionNitrogenSources.csv

4. The expression dataset with IDs corresponding to iPN730 were imported into MATLAB 2017b having COBRA Toolbox v3.0. The following code was run to map the expression
to the model using GIMME algorithm, conduct flux variability analysis and enrichment of altered reactions. 

Code: GIMMEcode.m
Inputs: iPN730_etac.xml, GeneExpressionNitrogenSources.csv
Outputs: EnrichedReactionsFluxRanges.csv, 

5. The following code in R was used to make the heatmap for the FVA data obtained in 4:

Code: Heatmap.R
Inputs: EnrichedReactionsFluxRanges.csv
