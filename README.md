## Description

To assess our filtering strategy, we analyzed a subset of samples from TCGA and compared fusions retained with annoFuse filtering and prioritization to those deemed the final call set in a previously published analysis by The Fusion Analysis Working Group (13). A group of 160 samples were randomly selected BLCA(10), BRCA(11), CESC(5), COAD(11),  ESCA(5) , GBM(7), HNSC(10), KIRP(9), LGG(9), LIHC(9), LUAD(5), LUSC(11), OV(9), PAAD(8), PCPG(2), PRAD(14), SARC(6), SKCM(9), TGCT(6), THCA(4). We first standardized STAR-Fusion and Arriba fusion calls using fusion_standardization, then performed artifact and QC filtering using fusion_filtering_QC. Using a default spanningDelta (spanningFragCount - JunctionReadCount) of 10, annoFuse retained only 40% of the fusions in the final call set (Table 3). Therefore, we visualized the distribution of spanningDelta (spanningFragCount - JunctionReadCount) across fusions called from the TCGA and PBTA cohorts to assess a cutoff for spanningDelta (Additional file 1 : Figure S2). We found that sensitivity reaches 96% at cutoff of 100 (Additional file 1 : Figure S34 and Table 3) for TCGA 50 to 76 bp read length RNAseq data. Therefore, we have implemented a default spanningDelta of 100 and made this a customizable input parameter. 


## Steps 
get_tcga_sample_set.R : get TCGA sample to run rnaseq and annoFuse
rnaseq_tcga_runs.R : rnaseq runs for TCGA samples
annofuse_tcga_runs.R : running arriba formatting and annotataion and annoFuse
TCGA_validation_plots.Rmd : TCGA validation plots and tables
