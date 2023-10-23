### Calculate copy number for each base
library(Rsamtools)
library(GenomicRanges)
chr.lengths=scanBamHeader('example.bam')[[1]]$targets
for (i in 1:length(sample_filenames)){
  filename=sample_filenames[i]
  reads.gr=readRDS(filename)
  mid.gr=reads.gr
  start(mid.gr)=start(reads.gr)+width(reads.gr)/2  
  width(mid.gr)=1
  for (c in 1:1){
    regions.gr=GRanges(seqnames = paste0('chr',as.roman(c)),
                       ranges = IRanges(start=seq(1,chr.lengths[c],1)-500,width = 1001))
    
    regions.df=data.frame(chr=seqnames(regions.gr),start=start(regions.gr),end=end(regions.gr),mid=start(regions.gr)+500,RPKM=countOverlaps(regions.gr,mid.gr,ignore.strand=T)*10^9/width(regions.gr)/length(mid.gr))
    write.table(regions.df,file=paste0(filename,'_RPKM_1bpby1001bp.txt'),row.names = F,append = F)
  }
  for (c in 2:16){
    regions.gr=GRanges(seqnames = paste0('chr',as.roman(c)),
                       ranges = IRanges(start=seq(1,chr.lengths[c],1)-500,width = 1001))
    
    regions.df=data.frame(chr=seqnames(regions.gr),start=start(regions.gr),end=end(regions.gr),mid=start(regions.gr)+500,RPKM=countOverlaps(regions.gr,mid.gr,ignore.strand=T)*10^9/width(regions.gr)/length(mid.gr))
    write.table(regions.df,file=paste0(filename,'_RPKM_1bpby1001bp.txt'),row.names = F,append = T,col.names = F)
  }
  
}

### Normalize nucleosome and TF scores by local copy number
library(data.table)
G1_RPKM.df=fread('alpha_factor_RPKM_1bpby1001bp.txt',header = T)
for (i in 1:length(sample_filenames)){
  RPKM.df=fread(paste0(sample_filenames[i],'_RPKM_1bpby1001bp.txt'),header = T)
  scale.v=(RPKM.df$RPKM+1)/(G1_RPKM.df$RPKM+1)
  nuc.df=fread(paste0(sample_filenames[i],'_nuc_cross_cor.txt'),header = T)
  nuc_norm.df=nuc.df
  nuc_norm.df$cross_cor=nuc.df$cross_cor/scale.v
  write.table(nuc_norm.df,file=paste0(sample_filenames[i],'_nuc_cross_cor_norm.txt'),row.names = F)
  tf.df=fread(paste0(sample_filenames[i],'_tf_cross_cor.txt'),header = T)
  tf_norm.df=tf.df
  tf_norm.df$cross_cor=tf.df$cross_cor/scale.v
  write.table(tf_norm.df,file=paste0(sample_filenames[i],'_tf_cross_cor_norm.txt'),row.names = F)
}