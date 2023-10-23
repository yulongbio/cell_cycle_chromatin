### combine corresponding timepoints and calculate RNA levels
rep1=sample_filenames_rep1
rep2=sample_filenames_rep2
library(Rsamtools)
library(GenomicRanges)
gene.df=read.csv('../Data/gene_coordinates_annotation.csv',header = T)
pos.idx=which(gene.df$strand=='+')
neg.idx=which(gene.df$strand=='-')
gene.gr=GRanges(seqnames = gene.df$chrom,
                ranges = IRanges(start=gene.df$txStart,end=gene.df$txEnd),
                strand = gene.df$strand)
ranges(gene.gr[pos.idx])=IRanges(start=gene.df$tss[pos.idx],end = gene.df$pas[pos.idx])
ranges(gene.gr[neg.idx])=IRanges(start=gene.df$pas[neg.idx],end=gene.df$tss[neg.idx])

reads.df=data.frame(name=gene.df$name,
                    chr=gene.df$chrom,
                    start=gene.df$txStart,end=gene.df$txEnd,strand=gene.df$strand,
                    alias=gene.df$alias,tss=gene.df$tss,pas=gene.df$pas)


for (i in 1:length(rep1)){
  filename1=sample_filenames_rep1[i]
  filename2=sample_filenames_rep2[i]
  p = ScanBamParam(what = c("rname","pos", "strand", "qwidth"))
  reads1.1= scanBam(filename1, param = p)
  IP1.gr = GRanges(seqnames = reads1.1[[1]][["rname"]],
                   ranges = IRanges(start = reads1.1[[1]][["pos"]], width = 1),
                   strand = reads1.1[[1]][["strand"]])
  idx = which(strand(IP1.gr) == "-")
  ranges(IP1.gr[idx]) = 
    IRanges(start = start(IP1.gr[idx]) + reads1.1[[1]][["qwidth"]][idx] - 1,
            width = 1)
  IP1.gr = IP1.gr[order(start(IP1.gr))]
  
  reads2.1= scanBam(filename2, param = p)
  IP2.gr = GRanges(seqnames = reads2.1[[1]][["rname"]],
                   ranges = IRanges(start = reads2.1[[1]][["pos"]], width = 1),
                   strand = reads2.1[[1]][["strand"]])
  idx = which(strand(IP2.gr) == "-")
  ranges(IP2.gr[idx]) = 
    IRanges(start = start(IP2.gr[idx]) + reads2.1[[1]][["qwidth"]][idx] - 1,
            width = 1)
  IP2.gr = IP2.gr[order(start(IP2.gr))]
  
  IP.gr=c(IP1.gr,IP2.gr)
  
  pos.gr=IP.gr[which(strand(IP.gr)=='+')]
  neg.gr=IP.gr[which(strand(IP.gr)=='-')]
  counts.v=rep(0,nrow(reads.df))
  
  
  reads.df=data.frame(reads.df,counts.v)
  reads.df[pos.idx,ncol(reads.df)]=countOverlaps(gene.gr[pos.idx],neg.gr,ignore.strand=T)
  reads.df[neg.idx,ncol(reads.df)]=countOverlaps(gene.gr[neg.idx],pos.gr,ignore.strand=T)
  
}
for (c in 9:ncol(reads.df)){
  reads.df[,c]=reads.df[,c]*10^9/sum(reads.df[,c])/(reads.df$end-reads.df$start)
}
colnames(reads.df)[9:ncol(reads.df)]=sample_filenames
write.csv(reads.df,file = 'RNA_TPM.csv',row.names = F)