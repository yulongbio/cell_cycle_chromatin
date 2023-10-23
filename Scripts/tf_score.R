### Derive the aggregate MNase count matrix for ABF1 binding sites as the reference TF matrix
library(GenomicRanges)
tf_sites.df=read.table('../Data/MacIsaac_sacCer3_liftOver.bed')
colnames(tf_sites.df)=c('chr','start','end','TF','strand')
tf_site.df=subset(tf_sites.df,TF=='ABF1')
# Specify sample file names ("sample_filenames") produced by sample_bam_by_fragment.R
counts_matrix.l=list()
for (l in 1:length(sample_filenames)){
  filename=sample_filenames[l]
  read.gr=readRDS(filename)
  mid.df=data.frame(chr=seqnames(read.gr),mid=start(read.gr)+floor(width(read.gr)/2),length=width(read.gr))
  chr4_counts.m=matrix(0,250,601)
  
  for (i in 1:nrow(tf_site.df)){
    
    counts.m=matrix(0,250,601)
    tf_site_mid=round(tf_site.df$start[i]/2+tf_site.df$end[i]/2,0)
    mid_nuc.df=subset(mid.df,chr==tf_site.df[i,]$chr&mid>=(tf_site_mid-300)&mid<=(tf_site_mid+300))
    for (size in 1:250){
      mid_nuc_size.df=subset(mid_nuc.df,length==size)
      counts.m[size,]=hist(mid_nuc_size.df$mid,seq(tf_site_mid-300,tf_site_mid+300+1,1),plot = F,right = F)$counts
    }
    chr4_counts.m=chr4_counts.m+counts.m
  }
  
  counts_matrix.l[[l]]=chr4_counts.m
}
names(counts_matrix.l)=sample_filenames
saveRDS(counts_matrix.l,file='ABF1_aggregate_matrix.l.RDS')

### Calculate TF 2D correlation score with reference TF kernel
tf_kernel<-function(mu,Sigma,x_width=151,ybot=20,ytop=250){
  library(mvtnorm)
  library(RColorBrewer)
  library(reshape)
  library(lattice)
  
  lengths <- seq(ybot, ytop)
  pos <- seq(1, x_width)
  dens.df <- as.data.frame(expand.grid(lengths, pos))
  colnames(dens.df)=c('length','pos')
  get.density <- function(row)
  {
    return(dmvnorm(c(row[1], row[2]), mean=mu, sigma=Sigma))
  }
  dens.df$density <- apply(dens.df, 1, get.density)
  kernel.df <- cast(dens.df, length ~ pos, value='density')
  kernel.mat <- as.matrix(kernel.df)
  x.scale=list(at=seq(1,x_width,10))
  y.scale=list(at=seq(ybot,ytop,40))
  levelplot(t(kernel.mat),xlab='Position',ylab='Fragment length',scales=list(x=x.scale,y=y.scale))
  
  return(kernel.mat)
}



calculate_tf_cross_cor_with_gr<-function(filename,chr,start_pos,end_pos,ybot=20,ytop=250,x_width=151,whole_genome='no'){
  library(GenomicRanges)
  library(Rsamtools)
  library(reticulate)
  use_python('/usr/lib/python2.7/')
  signal<-import('scipy.signal')
  reads.gr=readRDS(filename)
  broggard_matrix.l=readRDS('ABF1_aggregate_matrix.l.RDS')
  ID.m=broggard_matrix.l[[filename]]
  matrix_to_points<-function(v){
    point.v=vector()
    for (i in 1:length(v)){
      point.v=c(point.v,rep(i,v[i]))
    }
    return(point.v)
  }
  ID.m=ID.m[,(301-(x_width-1)/2):(301+(x_width-1)/2)]
  x.v=colSums(ID.m)
  x.t=matrix_to_points(x.v)
  y.v=rowSums(ID.m)
  y.t=matrix_to_points(y.v)
  bivn.m<-nuc_kernel(mu=c(quantile(y.t,0.3),mean(x.t)),x_width = 51,Sigma = matrix(c(var(y.t)/4,0,0,var(x.t)/4),2))
  
  
  if (whole_genome=='no') {
    range.gr=GRanges(seqnames = chr,
                     ranges = IRanges(start=start_pos-100,end=end_pos+100))
    
    mnase.df=data.frame(mid=start(reads.gr)+floor((width(reads.gr)-1)/2),length=width(reads.gr))
    mnase_count.m=matrix(0,ytop-ybot+1,end_pos-start_pos+1+200)
    colnames(mnase_count.m)=c((start_pos-100):(end_pos+100))
    for (i in 1:(ytop-ybot+1)){
      mnase_size_i.df=subset(mnase.df,length==(i+ybot-1)&mid<=(end_pos+100)&mid>=(start_pos-100))
      br=seq(start_pos-100,end_pos+100+1,1)
      mnase_count.m[i,]=hist(mnase_size_i.df$mid,br,plot = F,right = F)$counts
    }
    cross_cor.m=signal$correlate2d(mnase_count.m[,as.character((start_pos-0.5*(x_width-1)):(end_pos+0.5*(x_width-1)))],bivn.m,mode='valid')
    
    cross_cor.df=data.frame(chrom=chr,pos=c(start_pos:end_pos),cross_cor=cross_cor.m[1,])
    write.table(cross_cor.df,file=paste0(filename,'_',chr,'_',start_pos,'_',end_pos,'_tf_cross_cor.txt'),row.names = F)
    
  }
  else {
    chr_lengths=scanBamHeader('example.bam')[[1]]$targets
    cross_cor.df=data.frame(matrix(ncol=3,nrow=0))
    colnames(cross_cor.df)=c('chrom','pos','cross_cor')
    write.table(cross_cor.df,file=paste0(filename,'_tf_cross_cor.txt'),row.names = F)
    
    for (ch in 1:16){
      start_seq=seq(1,chr_lengths[ch],100000)
      for (s in 1:(length(start_seq)-1)){
        start_pos=start_seq[s];end_pos=start_pos+100000-1
        range.gr=GRanges(seqnames = paste0('chr',as.roman(ch)),
                         ranges = IRanges(start=start_pos-100,end=end_pos+100))
        reads_range.gr=subsetByOverlaps(reads.gr,range.gr)
        mnase.df=data.frame(mid=start(reads_range.gr)+floor((width(reads_range.gr)-1)/2),length=width(reads_range.gr))
        mnase_count.m=matrix(0,ytop-ybot+1,end_pos-start_pos+1+200)
        colnames(mnase_count.m)=c((start_pos-100):(end_pos+100))
        for (i in 1:(ytop-ybot+1)){
          mnase_size_i.df=subset(mnase.df,length==(i+ybot-1)&mid<=(end_pos+100)&mid>=(start_pos-100))
          br=seq(start_pos-100,end_pos+100+1,1)
          mnase_count.m[i,]=hist(mnase_size_i.df$mid,br,plot = F,right = F)$counts
        }
        cross_cor.m=signal$correlate2d(mnase_count.m[,as.character((start_pos-0.5*(x_width-1)):(end_pos+0.5*(x_width-1)))],bivn.m,mode='valid')
        
        cross_cor.df=data.frame(chrom=paste0('chr',as.roman(ch)),pos=c(start_pos:end_pos),cross_cor=cross_cor.m[1,])
        write.table(cross_cor.df,file=paste0(filename,'_sub_cross_cor.txt'),row.names = F,append = T,col.names = F)
      }
      start_pos=tail(start_seq,1);end_pos=chr_lengths[ch]
      range.gr=GRanges(seqnames = paste0('chr',as.roman(ch)),
                       ranges = IRanges(start=start_pos-100,end=end_pos+100))
      reads_range.gr=subsetByOverlaps(reads.gr,range.gr)
      mnase.df=data.frame(mid=start(reads_range.gr)+floor((width(reads_range.gr)-1)/2),length=width(reads_range.gr))
      mnase_count.m=matrix(0,ytop-ybot+1,end_pos-start_pos+1+200)
      colnames(mnase_count.m)=c((start_pos-100):(end_pos+100))
      for (i in 1:(ytop-ybot+1)){
        mnase_size_i.df=subset(mnase.df,length==(i+ybot-1)&mid<=(end_pos+100)&mid>=(start_pos-100))
        br=seq(start_pos-100,end_pos+100+1,1)
        mnase_count.m[i,]=hist(mnase_size_i.df$mid,br,plot = F,right = F)$counts
      }
      cross_cor.m=signal$correlate2d(mnase_count.m[,as.character((start_pos-0.5*(x_width-1)):(end_pos+0.5*(x_width-1)))],bivn.m,mode='valid')
      
      cross_cor.df=data.frame(chrom=paste0('chr',as.roman(ch)),pos=c(start_pos:end_pos),cross_cor=cross_cor.m[1,])
      write.table(cross_cor.df,file=paste0(filename,'_sub_cross_cor.txt'),row.names = F,append = T,col.names = F)
    }
  }
}
for (filename in sample_filenames){
  calculate_tf_cross_cor_with_gr(filename,whole_genome = 'yes',x_width = 51)
}