### Gausian Process modeling for cell cycle genes with log2 fold change from timepoint 10 min
nuc.df=read.csv('../Data/genebody_nuc_score.csv',header = T,check.names = F)
sub.df=read.csv('../Data/promoter_tf_score.csv',header = T,check.names = F)
RNA.df=read.csv('../Data/RNA_TPM.csv',header = T,check.names = F)
nuc_entropy.df=read.csv('../Data/genebody_nuc_entropy.csv',header = T,check.names = F)


nuc.m=as.matrix(nuc.df[,-c(1:8)])
rownames(nuc.m)=nuc.df$Gene_name
sub.m=as.matrix(sub.df[,-c(1:8)])
rownames(sub.m)=sub.df$Gene_name
RNA.m=as.matrix(RNA.df[,-c(1:8)])
rownames(RNA.m)=RNA.df$Gene_name
nuc_entropy.m=as.matrix(nuc_entropy.df[,-c(1:8)])
rownames(nuc_entropy.m)=RNA.df$Gene_name

RNA.m=log2(RNA.m+1)

tp=2
pad=1e-6
RNA.m=RNA.m-RNA.m[,tp]
nuc.m=log2((nuc.m+pad)/(nuc.m[,tp]+pad))
sub.m=log2((sub.m+pad)/(sub.m[,tp]+pad))
nuc_entropy.m=log2((nuc_entropy.m+pad)/(nuc_entropy.m[,tp]+pad))



cycle_gene.df=read.csv('../Data/RNA_f_score.csv',header = T)
cycle_gene.df=cycle_gene.df[order(cycle_gene.df$FDR,decreasing = F),]
cyc_gene.v=cycle_gene.df[which(cycle_gene.df$FDR<0.01),]$Gene_name

cycle_idx=match(cyc_gene.v,rownames(nuc.m))
library(kernlab)

cyc_R2.df=data.frame(RNA=rep(NA,13),nuc=rep(NA,13),sub=rep(NA,13),entropy=rep(NA,13),chr=rep(NA,13),full=rep(NA,13))
gene.idx=cycle_idx
for (f in c(1:6)){
  feature_idx=list(1,2,3,4,c(2:4),c(1:4))[[f]]
  for (t in 3:15){
    tp_idx=t
    x_mat=cbind(RNA.m[,tp_idx],nuc.m[,tp_idx],sub.m[,tp_idx],nuc_entropy.m[,tp_idx])
    
    fit<-gausspr(x=x_mat[gene.idx,feature_idx],y=RNA.m[gene.idx,tp_idx],scaled. =T,type='regression',kernel='rbfdot',cross=10,var=0.001)
    
    cyc_R2.df[t-2,f]=round(cor(RNA.m[gene.idx,tp_idx],predict(fit,x_mat[gene.idx,feature_idx]),use = 'complete.obs')^2,2)
  }
}
write.csv(cyc_R2.df,file='794_cyc_genes_Rsquare_10cv_l2fc_vs_tp2.csv')
# Randomly select 794 non-cycling genes
cyc_R2.df=data.frame(RNA=rep(NA,13),nuc=rep(NA,13),sub=rep(NA,13),entropy=rep(NA,13),chr=rep(NA,13),full=rep(NA,13))
set.seed(1)
x=1:nrow(RNA.m)
gene.idx=sample(x[!x%in%cycle_idx],length(cycle_idx),replace = F)

for (f in c(1:6)){
  feature_idx=list(1,2,3,4,c(2:4),c(1:4))[[f]]
  for (t in 3:15){
    tp_idx=t
    x_mat=cbind(RNA.m[,tp_idx],nuc.m[,tp_idx],sub.m[,tp_idx],nuc_entropy.m[,tp_idx])
    
    fit<-gausspr(x=x_mat[gene.idx,feature_idx],y=RNA.m[gene.idx,tp_idx],scaled. =T,type='regression',kernel='rbfdot',cross=10,var=0.001)
    
    cyc_R2.df[t-2,f]=round(cor(RNA.m[gene.idx,tp_idx],predict(fit,x_mat[gene.idx,feature_idx]),use = 'complete.obs')^2,2)
  }
}
write.csv(cyc_R2.df,file='794_noncyc_genes_Rsquare_10cv_l2fc_vs_tp2.csv')