###  Calculating Fourier score
fourier_score<-function(T,expr,time_series){
  f_score=rep(0,nrow(expr))
  p=2*pi/T*time_series
  for (i in 1:nrow(expr)){
    f_score[i]=sqrt((sum(sin(p)*expr[i,]))^2+(sum(cos(p)*expr[i,]))^2)
  }
  return(f_score)
}

### Fourier analysis of RNA level
RNA.df=read.csv('../Data/RNA_TPM.csv',header = T,check.names = F)
RNA.m=as.matrix(RNA.df[,-c(1:9)])
rownames(RNA.m)=RNA.df$Gene_name
for (i in 1:nrow(RNA.m)){
  RNA.m[i,]=(RNA.m[i,]-mean(RNA.m[i,]))/sd(RNA.m[i,])
}
RNA.m[is.na(RNA.m)]<-0
time_series=seq(10,140,10)
f_score.v=fourier_score(T=63,RNA.m,time_series)
# Permutation to calculate p-values
r=1e5
f_score.m=matrix(NA,nrow(RNA.m),r)
for (i in 1:r){
  set.seed(i)
  idx=sample(ncol(RNA.m),ncol(RNA.m),replace = F)
  RNA_temp.m=RNA.m[,idx]
  f_score.m[,i]=fourier_score(63,RNA_temp.m,time_series)
}
p.v=rep(NA,nrow(RNA.m))
for (i in 1:nrow(RNA.m)){
  p.v[i]=length(which(f_score.m[i,]>=f_score.v[i]))/r
}
q.v=p.adjust(p.v,method = 'fdr')
RNA_f_score.df=data.frame(Gene_ID=RNA.df$Gene_ID,Gene_name=rownames(RNA.m),Fourier_score=f_score.v,P_value=p.v,FDR=q.v)
write.csv(RNA_f_score.df,file='RNA_f_score.csv',row.names = F)

### Fourier analysis of genebody nuc score
nuc.df=read.csv('../Data/genebody_nuc_score.csv',header = T,check.names = F)
nuc.m=as.matrix(nuc.df[,-c(1:9)])
rownames(nuc.m)=nuc.df$Gene_name
for (i in 1:nrow(nuc.m)){
  nuc.m[i,]=(nuc.m[i,]-mean(nuc.m[i,]))/sd(nuc.m[i,])
}
nuc.m[is.na(nuc.m)]<-0
f_score.v=fourier_score(T=63,nuc.m,time_series)
# Permutation to calculate p-values
r=1e5
f_score.m=matrix(NA,nrow(nuc.m),r)
for (i in 1:r){
  set.seed(i)
  idx=sample(ncol(nuc.m),ncol(nuc.m),replace = F)
  nuc_temp.m=nuc.m[,idx]
  f_score.m[,i]=fourier_score(63,nuc_temp.m,time_series)
}
p.v=rep(NA,nrow(nuc.m))
for (i in 1:nrow(nuc.m)){
  p.v[i]=length(which(f_score.m[i,]>=f_score.v[i]))/r
}
q.v=p.adjust(p.v,method = 'fdr')
nuc_f_score.df=data.frame(Gene_ID=nuc.df$Gene_ID,Gene_name=rownames(nuc.m),Fourier_score=f_score.v,P_value=p.v,FDR=q.v)
write.csv(nuc_f_score.df,file='nuc_f_score.csv',row.names = F)

### Fourier analysis of promoter tf score
tf.df=read.csv('../Data/promoter_tf_score.csv',header = T,check.names = F)
tf.m=as.matrix(tf.df[,-c(1:9)])
rownames(tf.m)=tf.df$Gene_name
for (i in 1:nrow(tf.m)){
  tf.m[i,]=(tf.m[i,]-mean(tf.m[i,]))/sd(tf.m[i,])
}
tf.m[is.na(tf.m)]<-0
f_score.v=fourier_score(T=63,tf.m,time_series)
# Permutation to calculate p-values
r=1e5
f_score.m=matrix(NA,nrow(tf.m),r)
for (i in 1:r){
  set.seed(i)
  idx=sample(ncol(tf.m),ncol(tf.m),replace = F)
  tf_temp.m=tf.m[,idx]
  f_score.m[,i]=fourier_score(63,tf_temp.m,time_series)
}
p.v=rep(NA,nrow(tf.m))
for (i in 1:nrow(tf.m)){
  p.v[i]=length(which(f_score.m[i,]>=f_score.v[i]))/r
}
q.v=p.adjust(p.v,method = 'fdr')
tf_f_score.df=data.frame(Gene_ID=tf.df$Gene_ID,Gene_name=rownames(tf.m),Fourier_score=f_score.v,P_value=p.v,FDR=q.v)
write.csv(tf_f_score.df,file='tf_f_score.csv',row.names = F)