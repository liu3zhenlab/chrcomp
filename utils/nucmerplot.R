# R for plot Nucmer alignments

args <- commandArgs(trailingOnly=T)
nucmer_show <- args[1] # nucmer show-coordi output
band_color <- args[2] # color to be used for plotting
qryname <- args[3] # qry label
refname <- args[4] # ref label
outdir <- args[5] # PDF output directory
outpdf <- args[6] # PDF output file

###########################################################
#' module to determine xaxis
###########################################################
smartaxis <- function(maxnum) {
  numdigits <- nchar(maxnum)
  unit <- 10 ^ (numdigits - 1) / (2- round((maxnum / 10 ^ numdigits), 0)) # 1 or 5 e (numdigits - 1)
  subunit <- unit / 5 
  
  numsat <- unit * (0:10)
  numsat <- numsat[numsat < maxnum]
  
  if (numdigits >= 7) {
    numlabels <- numsat / 1000000
    label.scale <- "Mb"
  } else if (numdigits < 7 & numdigits >= 4) {
    numlabels <- numsat / 1000
    label.scale <- "kb"
  } else {
    numlabels <- numsat
    label.scale <- "bp"
  }
  
  subunits <- seq(0, maxnum, by = subunit)
  subunits <- subunits[!subunits %in% c(numsat, 0)] 
  # return
  list(numsat, numlabels, label.scale, subunits)
}

###########################################################
#' bandconnect
###########################################################
bandconnect <- function(topregion, bottomregion, topheight, bottomheight,
                        border=NULL, bandcol="grey") {
  
  topregion <- as.numeric(as.character(topregion))
  bottomregion <- as.numeric(as.character(bottomregion))
  
  ### module
  transform_curve <- function(startp, endp, npoint=1000) {
    ### computer transform_curve coordinates
    midp <- (startp + endp) / 2 
    beizer_value <- sqrt(1:npoint) / sqrt(npoint)  # sqrt as default
    
    if (startp[1] == endp[1]) {
      curve.x <- rep(startp[1], 2*npoint)
    } else {
      curve.x1 <- seq(startp[1], midp[1], by = (midp[1] - startp[1])/(npoint-1))
      curve.x2 <- seq(midp[1], endp[1], by = (endp[1] - midp[1])/(npoint-1))
      curve.x <- c(curve.x1, curve.x2)
    }
    
    if (startp[2] == endp[2]) {
      curve.y <- rep(startp[2], 2*npoint)
    } else {
      curve.y1 <- startp[2] - beizer_value * (startp[2] - midp[2])
      curve.y2 <- rev(endp[2] - beizer_value * (endp[2] - midp[2]))
      curve.y <- c(curve.y1, curve.y2)
    }
    list(x=curve.x, y=curve.y)
  }
  
  p1 <- transform_curve(c(topregion[1], topheight), c(bottomregion[1], bottomheight))
  p2 <- transform_curve(c(topregion[2], topheight), c(bottomregion[2], bottomheight))
  px <- c(p1$x, rev(p2$x))
  py <- c(p1$y, rev(p2$y))
  polygon(px, py, border=border, col=bandcol)
}

###########################################################
#' check if a color is valid
#' https://stackoverflow.com/questions/13289009/check-if-character-string-is-a-valid-color-representation
###########################################################
isColor <- function(col) {
  sapply(col, function(incol) {
    tryCatch(is.matrix(col2rgb(incol)), 
             error = function(e) FALSE)
  })
}

###########################################################
# nucmerplot
###########################################################
nucmerplot <- function(datafile, band_col="deepskyblue4", outpath=".", imageoutfile) {
  
  # check if the color is valid
  if (!isColor(band_col)) {
    cat(band_col, "is not valid. Plot using the default\n")
    band_col="deepskyblue4"
  }
  
  # data of alignments, ref and qry
  aln <- read.delim(datafile, stringsAsFactor=F)
  stopifnot(nrow(aln)>=1)
  
  rlen <- unique(as.numeric(as.character(aln$rlen)))
  qlen <- unique(as.numeric(as.character(aln$qlen)))
  ref <- unique(aln$ref)
  qry <- unique(aln$qry)

  ### plot
  outpdffile <- paste0(outpath, "/", imageoutfile)
  pdf(outpdffile, width=4.5, height=3)
  
  # plot setting
  par(mar=c(2.2, 0.2, 0.2, 0.2), mgp=c(1.5,0.1,0))
  
  # x-range
  max.xsize <- max(rlen, qlen)
  xrange <- c(-max.xsize / 50, max.xsize + max.xsize / 50)

  # connection colors
  high_identity_col <- band_col
  low_identity_col <- "grey96"
  colfunc <- colorRampPalette(c(low_identity_col, high_identity_col))
  color_num <- 20
  allcols <- colfunc(color_num)

  # identity
  identity_high <- max(aln$identity)
  identity_low <- min(aln$identity)
  if (identity_high - identity_low < 5) {
    identity_low <- identity_high - 5
  }

  identity20 <- seq(identity_low, identity_high, by=(identity_high - identity_low)/(color_num -1))

  # y-range
  yrange <- c(0, 1)

  # plot canvas
  plot(NULL, NULL, type="n", xlim=xrange, ylim=yrange,
       xlab="", ylab="", bty="n", xaxt="n", yaxt="n")
  
  # x-axis
  xaxis <- smartaxis(max.xsize)
  axis(1, at=xaxis[[1]], labels=xaxis[[2]], tick=F, col="gray50", cex.axis=0.8)
  mtext(paste0("position (", xaxis[[3]], ")"), side=1, line=1)
  
  # vertical lines
  #abline(v=xaxis[[1]], col="gray70", lwd=0.8)
  xaxis_thick <- xaxis[[1]]
  for (each_thick in xaxis_thick) {
  	lines(c(each_thick, each_thick), c(-0.025, 0.825), col="gray70", lwd=0.8, xpd=T) 
  }
  
  #abline(v=xaxis[[4]], col="gray90", lwd=0.5)
  xaxis_thin <- xaxis[[4]]
  for (each_thin in xaxis_thin) {
  	lines(c(each_thin, each_thin), c(-0.025, 0.825), col="gray90", lwd=0.5, xpd=T)
  }

  # region lines
  upper_ycenter <- 0.7
  upper_adj <- 0.005
  lower_ycenter <- 0.1
  lower_adj <- 0.005
  rect(0, upper_ycenter - upper_adj, qlen, upper_ycenter + upper_adj,
       col="gray70", border="gray70")
  rect(0, lower_ycenter - lower_adj, rlen, lower_ycenter + lower_adj,
       col="gray70", border="gray70")
  
  # connections
  for (i in 1:nrow(aln)) {
    band_gradient_col <- allcols[which.min(abs(identity20 - aln[i, "identity"]))]
    bandconnect(topregion=aln[i, 3:4], topheight=upper_ycenter - upper_adj - 0.005,
                bottomregion=aln[i, 1:2], bottomheight=lower_ycenter + lower_adj + 0.005,
                border=NA, bandcol=band_gradient_col)
  }
  
  # text labels
  upper_labels <- paste0(qryname, ":", qry)
  lower_labels <- paste0(refname, ":", ref)
  text(x= -max.xsize / 50, y=upper_ycenter + upper_adj + 0.04,
       labels=upper_labels, cex=0.8, xpd=T, pos=4)
  text(x= -max.xsize / 50, y=lower_ycenter - lower_adj - 0.055,
       labels=lower_labels, cex=0.8, xpd=T, pos=4)

  ### identity legends
  bar_ypos <- 0.96
  barlen <- max.xsize / 4
  barstep <- barlen / color_num
  barheight <- 0.04
  text(max.xsize - barlen/2, bar_ypos, labels = "identity", cex=0.8, pos=1, xpd=T)  ### plot LD name
  barlabels <- c(identity_low, identity_high)
  barlabels.num <- floor(barlabels * color_num)
  identity_low <- round(identity_low, 1)
  text(max.xsize-barlen, bar_ypos+0.01, labels=identity_low, cex=0.7, pos=1)  ### plot low identity
  identity_high <- round(identity_high, 1)
  if (identity_high == 100) {
  	identity_high <- round(identity_high, 0)
  }
  text(max.xsize, bar_ypos+0.01, labels=identity_high, cex=0.7, pos=1)  ### plot low identity
  
  # color bars for identity legends
  col_barpos <- max.xsize-barlen
  for (i in 1:color_num) {
    rect(col_barpos, bar_ypos,
	     col_barpos + barstep, bar_ypos + barheight,
         border=NA, col=allcols[i])
    col_barpos <- col_barpos + barstep
  }
  
  # close image output
  dev.off()
}

# execute plotting
nucmerplot(datafile=nucmer_show, band_col=band_color, outpath=outdir,
           imageoutfile=outpdf)

