### Average nucleosome score for the first 500bp of gene body
library(data.table)
gene.df=read.csv('../Data/gene_coordinates_annotation.csv',header = T)
nuc.df=data.frame(name=gene.df$name,
                  chr=gene.df$chrom,
                  start=gene.df$txStart,end=gene.df$txEnd,strand=gene.df$strand,
                  alias=gene.df$alias,tss=gene.df$tss,pas=gene.df$pas)
nuc_score.m=matrix(0,nrow(nuc.df),length(sample_filenames))
for (filename in sample_filenames){
  cross_cor.df=fread(paste0(filename,'_nuc_cross_cor_norm.txt'),header = T)
  chr_cross_cor.l=list()
  for (i in 1:16){
    chr_cross_cor.l[[i]]=cross_cor.df[chrom==paste0('chr',as.roman(i))]$cross_cor
  }
  names(chr_cross_cor.l)=sprintf('chr%s',as.roman(1:16))
  for (i in 1:nrow(nuc.df)){
    if (nuc.df[i,]$strand=='+'){
      gene_cross_cor.v=chr_cross_cor.l[[nuc.df$chr[i]]][nuc.df$tss[i]:min(nuc.df$tss[i]+500,nuc.df$end[i])]
      nuc_score.m[i,which(DM_IDs==DM_ID)]=mean(gene_cross_cor.v)
    }
    
    if (nuc.df[i,]$strand=='-'){
      gene_cross_cor.v=chr_cross_cor.l[[nuc.df$chr[i]]][max(nuc.df$tss[i]-500,nuc.df$start[i]):nuc.df$tss[i]]
      nuc_score.m[i,which(DM_IDs==DM_ID)]=mean(gene_cross_cor.v)
    }
  }
  
}
nuc.df=cbind(nuc.df,nuc_score.m)
colnames(nuc.df)[9:ncol(nuc.df)]=sample_filenames
write.csv(nuc.df,file=paste0('genebody_nuc_score.csv'),row.names = F)

### TF score for promoter
library(data.table)
gene.df=read.csv('../Data/gene_coordinates_annotation.csv',header = T)
tf.df=data.frame(name=gene.df$name,
                  chr=gene.df$chrom,
                  start=gene.df$txStart,end=gene.df$txEnd,strand=gene.df$strand,
                  alias=gene.df$alias,tss=gene.df$tss,pas=gene.df$pas)

tf_score.m=matrix(0,nrow(tf.df),length(sample_filenames))
for (filename in sample_filenames){
  cross_cor.df=fread(paste0(filename,'_tf_cross_cor_norm.txt'),header = T)
  chr_cross_cor.l=list()
  for (i in 1:16){
    chr_cross_cor.l[[i]]=cross_cor.df[chrom==paste0('chr',as.roman(i))]$cross_cor
  }
  names(chr_cross_cor.l)=sprintf('chr%s',as.roman(1:16))
  for (i in 1:nrow(tf.df)){
    if (tf.df[i,]$strand=='+'){
      gene_cross_cor.v=chr_cross_cor.l[[tf.df$chr[i]]][max(tf.df$tss[i]-250,1):(tf.df$tss[i])]
      tf_score.m[i,which(DM_IDs==DM_ID)]=mean(gene_cross_cor.v)
    }
    
    if (tf.df[i,]$strand=='-'){
      gene_cross_cor.v=chr_cross_cor.l[[tf.df$chr[i]]][(tf.df$tss[i]):(tf.df$tss[i]+250)]
      tf_score.m[i,which(DM_IDs==DM_ID)]=mean(gene_cross_cor.v)
    }
  }
  
}
tf.df=cbind(tf.df,tf_score.m)
colnames(tf.df)[9:ncol(tf.df)]=sample_filenames
write.csv(tf.df,file=paste0('promoter_tf_score.csv'),row.names = F)

### Entropy of nuc score for the first 500bp of gene body
library(data.table)
library(entropy)
gene.df=read.csv('../Data/gene_coordinates_annotation.csv',header = T)
nuc.df=data.frame(name=gene.df$name,
                  chr=gene.df$chrom,
                  start=gene.df$txStart,end=gene.df$txEnd,strand=gene.df$strand,
                  alias=gene.df$alias,tss=gene.df$tss,pas=gene.df$pas)
nuc_score.m=matrix(0,nrow(nuc.df),length(sample_filenames))

for (filename in sample_filenames){
  cross_cor.df=fread(paste0(filename,'_nuc_cross_cor_norm.txt'),header = T)
  chr_cross_cor.l=list()
  for (i in 1:16){
    chr_cross_cor.l[[i]]=cross_cor.df[chrom==paste0('chr',as.roman(i))]$cross_cor
  }
  names(chr_cross_cor.l)=sprintf('chr%s',as.roman(1:16))
  for (i in 1:nrow(nuc.df)){
    if (nuc.df[i,]$strand=='+'){
      gene_cross_cor.v=chr_cross_cor.l[[nuc.df$chr[i]]][nuc.df$tss[i]:min(nuc.df$tss[i]+500,nuc.df$end[i])]
      nuc_score.m[i,which(DM_IDs==DM_ID)]=entropy.empirical(gene_cross_cor.v)
    }
    
    if (nuc.df[i,]$strand=='-'){
      gene_cross_cor.v=chr_cross_cor.l[[nuc.df$chr[i]]][max(nuc.df$tss[i]-500,nuc.df$start[i]):nuc.df$tss[i]]
      nuc_score.m[i,which(DM_IDs==DM_ID)]=entropy.empirical(gene_cross_cor.v)
    }
  }
  
}
nuc.df=cbind(nuc.df,nuc_score.m)
colnames(nuc.df)[9:ncol(nuc.df)]=sample_filenames
write.csv(nuc.df,file=paste0('genebody_nuc_entropy.csv'),row.names = F)
