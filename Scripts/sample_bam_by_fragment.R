sample_bam_by_fragment<-function(min_size,max_size,sample_filenames){
  library(Rsamtools)
  library(GenomicRanges)
  count.df=data.frame(frag_size=c(min_size:max_size))
  p=ScanBamParam(what = c("rname","pos", "strand", "isize"))
  for (i in 1:length(sample_filenames)){
    filename=sample_filenames[i]
    reads.l=scanBam(filename,param = p)
    count.v=vector()
    for (l in min_size:max_size){
      count.v=c(count.v,length(which(reads.l[[1]][['isize']]==l)))
    }
    count.df=cbind(count.df,count.v)
    colnames(count.df)[ncol(count.df)]=basename(filename)
  }
  count_min.v=apply(count.df[,2:ncol(count.df)],1,min)
  count_min.df=data.frame(count.df$frag_size,count_min.v)
  for (i in 1:length(DM_id)){
    filename=sample_filenames[i]
    reads.l=scanBam(filename,param = p)
    reads.gr<-GRanges(seqnames=reads.l[[1]][['rname']],
                      ranges=IRanges(start=reads.l[[1]][['pos']],width=reads.l[[1]][['isize']]),
                      strand=reads.l[[1]][['strand']])
    set.seed(1)
    sample_idx=vector()
    for (l in min_size:max_size){
      l_idx=which(width(reads.gr)==l)
      l_idx_sample=sample(l_idx,count_min.df[which(count_min.df[,1]==l),2])
      sample_idx=c(sample_idx,l_idx_sample)
    }
    sample.gr=reads.gr[sample_idx]
    saveRDS(sample.gr,file=paste0('../Data/',filename,'_sample_gr.RDS'))
  }
}
sample_bam_by_fragment(20,250,sample_filenames)  
