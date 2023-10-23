densColors<-function(x, y = NULL, nbin = 128, bandwidth, transformation = function(x) x^1, colramp = colorRampPalette(blues9), z_factor = 1)
{
  library(RColorBrewer)
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
  #print(length(dens))
  colpal <- cut(dens, length(dens), labels = FALSE)
  colpal<-ceiling(as.integer(colpal/z_factor)+1)
  #print(length(colpal))
  cols <- rep(NA_character_, length(select))
  #print(length(cols))
  cols[select] <- colramp(length(dens))[colpal]
  #print(length(cols[select]))
  #print(length(colramp(length(dens))[colpal]))
  cols
}

