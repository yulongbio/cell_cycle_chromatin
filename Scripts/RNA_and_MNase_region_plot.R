### Plot histogram for RNA-seq data and typhoon plot for MNase-seq data to visualize specific genome region (e.g., Figure 1C)
library('Rsamtools')
library('GenomicRanges')
library('pheatmap')
library('data.table')
library('RColorBrewer')
# Function to plot gene schematic
make_gene_schematic = function(feature_chr,
                               feature_start,
                               feature_end,
                               y_low = 0,
                               y_high = 1,
                               cex_title = 1,
                               bg_type = "white",
                               proteinCoding = T,
                               geneName = T,
                               omit_genes = NA,
                               x_pos_title = 50)
{
  # Set up the plot
  plot(
    0,
    0,
    type = "n",
    bty = "n",
    bg = bg_type,
    xlim = c(feature_start, feature_end),
    xaxs = "i",
    xaxt = "n",
    ylim = c(0, 1),
    yaxs = "i",
    yaxt = "n",
    ann = F
  )
  
  # Load the gene dataframe
  gene.df = read.csv("../Data/sacCer3_genes_for_making_schematic.csv",
                     header = T)
  
 
  # Convert to a GenomicRanges object
  gene.gr = GRanges(
    seqnames = gene.df$chr,
    ranges = IRanges(start = gene.df$left_coord, end = gene.df$right_coord),
    strand = gene.df$strand
  )
  names(gene.gr) = gene.df$alias
  

  # Create the feature gr
  feature.gr = GRanges(seqnames = feature_chr,
                       ranges = IRanges(start = feature_start, end = feature_end))
  
  # Find the overlaps
  gene_overlaps.df = as.data.frame(as.matrix(findOverlaps(feature.gr, gene.gr)))
  
  if (any(nrow(gene_overlaps.df))) {
    # Enter in the genes
    for (i in 1:nrow(gene_overlaps.df)) {
      plot_gene(
        gene.df[gene_overlaps.df$subjectHits[i], ],
        y_low,
        y_high,
        feature_start,
        feature_end,
        cex_title,
        geneName,
        x_pos_title
      )
    }
    
  }
}
  
# Get dataframe of fragment midpoint
get_midpoint_dataframe_from_gr<-function(filename,chr,start_pos,end_pos){
  range.gr=GRanges(seqnames = chr,
                   ranges = IRanges(start=start_pos-500,end=end_pos+500))
  reads.gr=readRDS(file = filename)
  reads.gr=subsetByOverlaps(reads.gr,range.gr)
  reads.gr=reads.gr[width(reads.gr)<=250]
  mnase.df=data.frame(mid=start(reads.gr)+floor((width(reads.gr)-1)/2),length=width(reads.gr))
  mnase_mid.df=subset(mnase.df,mid>=start_pos&mid<=end_pos)
  return(mnase_mid.df)
}

# Create density color
densColors<-function(x, y = NULL, nbin = 128, bandwidth, transformation = function(x) x^1, colramp = colorRampPalette(blues9), z_factor = 1)
{
  xy <- xy.coords(x, y)
  select <- is.finite(xy$x) & is.finite(xy$y)
  x <- cbind(xy$x, xy$y)[select, ]
  map <- grDevices:::.smoothScatterCalcDensity(x, nbin, bandwidth)
  mkBreaks <- function(u) u - diff(range(u))/(length(u) - 1)/2
  xbin <- cut(x[, 1], mkBreaks(map$x1), labels = FALSE)
  ybin <- cut(x[, 2], mkBreaks(map$x2), labels = FALSE)
  dens <- map$fhat[cbind(xbin, ybin)]
  dens[is.na(dens)] <- 0
  dens[] <- transformation(dens)
  colpal <- cut(dens, length(dens), labels = FALSE)
  colpal<-ceiling(as.integer(colpal/z_factor)+1)
  cols <- rep(NA_character_, length(select))
  cols[select] <- colramp(length(dens))[colpal]
  cols
}

# Plot RNA and MNase profiles
difference_plot_RNA_MNase_by_region_from_gr <-
  function(RNA_filenames,
           MNase_filenames,
           chr,
           pos,
           left_window_size,
           right_window_size,
           time_labels = NULL,
           gene_label = NULL,
           strand) {
    range.gr = GRanges(
      seqnames = chr,
      ranges = IRanges(
        start = pos - left_window_size - 2000,
        width = left_window_size + 2000
      )
    )
    
    par(
      mfcol = c(length(RNA_filenames) + 1, 2),
      mar = c(1, 3, 1, 2),
      oma = c(1, 1, 0, 0),
      cex = 1.2
    )
    make_gene_schematic(chr, pos - left_window_size, pos + right_window_size)
    for (i in 1:(length(RNA_filenames) - 1)) {
      RNA_filename = RNA_filenames[i]
      p = ScanBamParam(
        what = c("rname", "pos", "strand", "qwidth"),
        which = GRanges(
          seqnames = chr,
          ranges = IRanges(
            start = pos - left_window_size - 2000,
            end = pos + right_window_size + 2000
          )
        )
      )
      reads.l = scanBam(RNA_filename, param = p)
      total_read_counts = countBam(RNA_filename)$records
      IP.gr = GRanges(
        seqnames = reads.l[[1]][["rname"]],
        ranges = IRanges(start = reads.l[[1]][["pos"]], width = 1),
        strand = reads.l[[1]][["strand"]]
      )
      idx = which(strand(IP.gr) == "-")
      ranges(IP.gr[idx]) =
        IRanges(start = start(IP.gr[idx]) + reads.l[[1]][["qwidth"]][idx] - 1,
                width = 1)
      pos.gr = IP.gr[strand(IP.gr) == '+']
      neg.gr = IP.gr[strand(IP.gr) == '-']
      ranges.gr = GRanges(seqnames = chr,
                          ranges = IRanges(
                            start = seq(pos - left_window_size, pos + right_window_size, 10),
                            width = 50
                          ))
      pos_RPKM = countOverlaps(ranges.gr, neg.gr, ignore.strand = T) * 10 ^
        9 / width(ranges.gr) / total_read_counts
      neg_RPKM = countOverlaps(ranges.gr, pos.gr, ignore.strand = T) * 10 ^
        9 / width(ranges.gr) / total_read_counts
      if (strand == '+') {
        plot(
          start(ranges.gr),
          pos_RPKM,
          type = 'h',
          col = rgb(1, 0, 0, 0.5),
          ylim = c(-100, 300),
          xlab = '',
          ylab = '',
          xaxs = 'i',
          yaxs = 'i',
          main = time_labels[i],
          xaxt = 'n'
        )
      }
      if (strand == '-') {
        plot(
          start(ranges.gr),
          pos_RPKM,
          type = 'h',
          col = rgb(1, 0, 0, 0.5),
          ylim = c(-3000, 200),
          xlab = '',
          ylab = '',
          xaxs = 'i',
          yaxs = 'i',
          main = time_labels[i],
          xaxt = 'n'
        )
      }
      points(start(ranges.gr),
             neg_RPKM * (-1),
             type = 'h',
             col = rgb(0, 0, 1, 0.5))
      
    }
    i = length(RNA_filenames)
    RNA_filename = RNA_filenames[i]
    p = ScanBamParam(
      what = c("rname", "pos", "strand", "qwidth"),
      which = GRanges(
        seqnames = chr,
        ranges = IRanges(
          start = pos - left_window_size - 2000,
          end = pos + right_window_size + 2000
        )
      )
    )
    reads.l = scanBam(RNA_filename, param = p)
    total_read_counts = countBam(RNA_filename)$records
    IP.gr = GRanges(
      seqnames = reads.l[[1]][["rname"]],
      ranges = IRanges(start = reads.l[[1]][["pos"]], width = 1),
      strand = reads.l[[1]][["strand"]]
    )
    idx = which(strand(IP.gr) == "-")
    ranges(IP.gr[idx]) =
      IRanges(start = start(IP.gr[idx]) + reads.l[[1]][["qwidth"]][idx] - 1,
              width = 1)
    pos.gr = IP.gr[strand(IP.gr) == '+']
    neg.gr = IP.gr[strand(IP.gr) == '-']
    ranges.gr = GRanges(seqnames = chr,
                        ranges = IRanges(
                          start = seq(pos - left_window_size, pos + right_window_size, 10),
                          width = 50
                        ))
    pos_RPKM = countOverlaps(ranges.gr, neg.gr, ignore.strand = T) * 10 ^
      9 / width(ranges.gr) / total_read_counts
    neg_RPKM = countOverlaps(ranges.gr, pos.gr, ignore.strand = T) * 10 ^
      9 / width(ranges.gr) / total_read_counts
    if (strand == '+') {
      plot(
        start(ranges.gr),
        pos_RPKM,
        type = 'h',
        col = rgb(1, 0, 0, 0.5),
        ylim = c(-100, 300),
        xlab = '',
        ylab = '',
        xaxs = 'i',
        yaxs = 'i',
        main = time_labels[i],
        xaxt = 'n'
      )
    }
    if (strand == '-') {
      plot(
        start(ranges.gr),
        pos_RPKM,
        type = 'h',
        col = rgb(1, 0, 0, 0.5),
        ylim = c(-3000, 200),
        xlab = '',
        ylab = '',
        xaxs = 'i',
        yaxs = 'i',
        main = time_labels[i],
        xaxt = 'n'
      )
    }
    axis(
      1,
      at = c(pos - left_window_size, pos, pos + 1000, pos + right_window_size),
      labels = c(-1 * left_window_size, 0, 1000, right_window_size)
    )
    points(start(ranges.gr),
           neg_RPKM * (-1),
           type = 'h',
           col = rgb(0, 0, 1, 0.5))
    make_gene_schematic(chr, pos - left_window_size, pos + right_window_size)
    
    for (i in 1:(length(MNase_filenames) - 1)) {
      filename = MNase_filenames[i]
      mid.df = get_midpoint_dataframe_from_gr(filename,
                                              chr,
                                              pos - left_window_size,
                                              pos + right_window_size)
      densCols = densColors(
        x = mid.df$mid,
        y = mid.df$length,
        nbin = 128,
        bandwidth = c(10, 15),
        transformation = function(x)
          x ^ 0.5,
        colramp = colorRampPalette(brewer.pal(9, "Reds"))
      )
      plot(
        mid.df,
        col = densCols,
        xaxs = 'i',
        yaxs = 'i',
        pch = '.',
        cex = 2,
        ylim = c(20, 250),
        xlim = c(pos - left_window_size, pos + right_window_size),
        main = time_labels[i],
        xlab = '',
        ylab = '',
        yaxt = 'n',
        xaxt = 'n'
      )
      axis(
        2,
        at = c(50, 150, 250),
        labels = c(50, 150, 250),
        tick = T
      )
    }
    
    i = length(MNase_filenames)
    filename = MNase_filenames[i]
    mid.df = get_midpoint_dataframe_from_gr(filename, chr, pos - left_window_size, pos +
                                              right_window_size)
    densCols = densColors(
      x = mid.df$mid,
      y = mid.df$length,
      nbin = 128,
      bandwidth = c(10, 15),
      transformation = function(x)
        x ^ 0.5,
      colramp = colorRampPalette(brewer.pal(9, "Reds"))
    )
    plot(
      mid.df,
      col = densCols,
      xaxs = 'i',
      yaxs = 'i',
      pch = '.',
      cex = 2,
      ylim = c(20, 250),
      xlim = c(pos - left_window_size, pos + right_window_size),
      main = time_labels[i],
      xlab = '',
      ylab = '',
      yaxt = 'n',
      xaxt = 'n'
    )
    axis(
      1,
      at = c(pos - left_window_size, pos, pos + 1000, pos + right_window_size),
      labels = c(-1 * left_window_size, 0, 1000, right_window_size)
    )
    axis(
      2,
      at = c(50, 150, 250),
      labels = c(50, 150, 250),
      tick = T
    )
  }
