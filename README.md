## Overview 
Code for the cell cycle chromatin paper. Relevant data sets are under `Data` directory. Codes are under `Scripts` directory.  
  
## Scripts usage   
`rna_tpm.R`: combine RNA-seq datasets of corresponding time points between replicates and calculate the TPM of each gene. 
  
`sample_bam_by_fragment.R`: sample bam files to obtain equal fragment number among samples for each fragment size. The sampled data are stored as GenomicRange objects and used for downstream analysis.  
  
`nucleosome_score.R`: calculate the whole-genome nucleosome scores for a sample using its GenomicRange object produced by `sample_bam_by_fragment.R`.  
  
`tf_score.R`: calculate the whole-genome TF scores for a sample using its GenomicRange object produced by `sample_bam_by_fragment.R`. 
  
`copy_number_correction.R`: normalize nucleosome and TF scores by local DNA copy number.  
  
`chromatin_featuer_profiling.R`: calculate the three types of chromatin features - nucleosome score of the first 500 bp of gene body; TF score for the promoter; and nucleosome entropy of the first 500 bp of gene body. 
  
`fourier_analysis.R`: perform Fourier analysis to determine periodicity of RNA levels, genebody nucleosome scores, and promoter tf scores. 
  
`RNA_and_MNase_region_plot.R`: plot RNA-seq data with histogram and MNase-seq data with typhoon plot for visualization of specific genome regions.  
  
`RNA_and_MNase_gene_heatmap.R`: plot heatmaps of RNA-seq data and MNase-seq data for specific genes.  


`gausian_process.R`: construction of Gaussian process models for predicting RNA levels with chromatin features.  
  
