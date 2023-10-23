### Plot heatmap of chromatin scores for specific genes (e.g., Figure 1D), using nucleosome and TF cross correlation scores
library(pheatmap)
library(data.table)
library(grid)
library(gridExtra)
difference_heatmap_MNase_by_region_from_gr<-function(MNase_filenames,chr,pos,left_window_size,right_window_size,time_labels=NULL,gene_label=NULL,gene_strand='+',legend_low=-0.1,legend_high=0.1){
  
  
  
  source('~/functions/create_typhoon_plot_functions.R')
  
  
  cross_cor.m=matrix(0,length(MNase_filenames),length((pos-left_window_size):(pos+right_window_size)))
  
  for (i in 1:(length(MNase_filenames))){
    
    nuc_cross_cor.df=fread(MNase_filenames[i],'_nuc_cross_cor_norm.txt'),header = T)
colnames(nuc_cross_cor.df)[2]='position'
region_nuc_cross_cor.df=subset(nuc_cross_cor.df,chrom==chr&position>=(pos-left_window_size)&position<=(pos+right_window_size))
sub_cross_cor.df=fread(paste0(MNase_filenames[i],'_sub_cross_cor_norm.txt'),header = T)
colnames(sub_cross_cor.df)[2]='position'

region_sub_cross_cor.df=subset(sub_cross_cor.df,chrom==chr &position>=(pos-left_window_size)&position<=(pos+right_window_size))
cross_cor.m[i,]=region_nuc_cross_cor.df$cross_cor-region_sub_cross_cor.df$cross_cor
  }
  t= c(paste0('-',left_window_size),rep('',left_window_size-1),'TSS',rep('',right_window_size-1),paste0('+',right_window_size))
  
  assignInNamespace(
    x = "draw_colnames",
    value = "draw_colnames_45",
    ns = asNamespace("pheatmap")
  )
  
  colnames(cross_cor.m)=t
  saveRDS(cross_cor.m,file=paste0(gene_label,'cross_cor.m.rds'))
  MNase_ph=pheatmap(-cross_cor.m,cluster_rows = F,cluster_cols = F,show_colnames = T,labels_row = time_labels,breaks = seq(legend_low,legend_high,length.out = 101),legend_breaks = c(-0.05,0,0.05,legend_high),legend_labels = c(-0.05,0,0.05,'chr score'),main=paste0(gene_label,' chromatin score'),fontsize = 10,cellwidth = 0.25,cellheight = 12,silent = F)
  saveRDS(MNase_ph,file=paste0(gene_label,'MNase_ph.rds')) 
  total_TPM.df=read.csv('../Data/RNA_TPM.csv',header = T,check.names = F)
  total_TPM.m=as.matrix(total_TPM.df[,9:23])
  rownames(total_TPM.m)=total_TPM.df$alias
  
  for (i in 1:nrow(total_TPM.m)){
    total_TPM.m[i,]=(total_TPM.m[i,]+1)/(mean(total_TPM.m[i,])+1)
    
  }
  
  total_TPM.m=log2(total_TPM.m)
  gene_TPM.m=as.matrix(total_TPM.m[gene_label,])
  colnames(gene_TPM.m)=gene_label
  gene_TPM_ph=pheatmap(gene_TPM.m,cluster_rows = F,cluster_cols = F,show_rownames=F,show_colnames = T,breaks = seq(-1.2,1.2,length.out = 101),legend_breaks = c(-1,-0.5,0,0.5,1,1.2),legend_labels =c(-1,-0.5,0,0.5,1,'log2 fc') ,main='RNA',fontsize = 10,cellwidth=20,cellheight = 12,silent = F)
  
  
  grid.arrange(gene_TPM_ph[[4]],MNase_ph[[4]],ncol=2,widths=c(1,5),respect=T,padding=unit(0.0,'npc'))   
  downViewport('arrange.1-2-1-2')
  downViewport('matrix.4-3-4-3')
  
  grid.lines((left_window_size+1)/(left_window_size+right_window_size+1),c(0,1))
  
  gene.df=return_gene_df(chr,pos-left_window_size,pos+right_window_size)
  for (i in 1:nrow(gene.df)){
    if (gene.df[i,]$strand=='+'){
      grid.rect(x=(gene.df[i,]$left_coord-(pos-left_window_size)+1)/(right_window_size+left_window_size+1),y=unit(-20,'points'),width = (min(gene.df[i,]$right_coord,pos+right_window_size)-gene.df[i,]$left_coord+1)/(right_window_size+left_window_size+1),height = unit(15,'points'),just=0,gp=gpar(fill=rgb(1,0,0,0.5)))
      grid.text(gene.df[i,]$alias,x=(gene.df[i,]$left_coord-(pos-left_window_size)+1)/(right_window_size+left_window_size+1),y=unit(-20,'points'),hjust = -0.2)
    }
    
    if (gene.df[i,]$strand=='-'){
      grid.rect(x=(max(gene.df[i,]$left_coord,pos-left_window_size)-(pos-left_window_size)+1)/(right_window_size+left_window_size+1),y=unit(-40,'points'),width = (gene.df[i,]$right_coord-max(gene.df[i,]$left_coord,pos-left_window_size)+1)/(right_window_size+left_window_size+1),height = unit(15,'points'),just=0,gp = gpar(fill=rgb(0,0,1,0.5)))
      grid.text(gene.df[i,]$alias,x=(gene.df[i,]$right_coord-(pos-left_window_size)+1)/(right_window_size+left_window_size+1),y=unit(-40,'points'),hjust =1.2)
    }
  }
}