#####################################################################################
#' module to determine xaxis
#####################################################################################
smart.axis <- function(maxnum) {
  numdigits <- nchar(maxnum)
  unit <- 10 ^ (numdigits - 1) / (2- round((maxnum / 10 ^ numdigits), 0)) # 1 or 5 e (numdigits - 1)
  subunit <- unit / 5
  
  numsat <- unit * (0:10)
  numsat <- numsat[numsat < maxnum]
  
  if (numdigits >= 6) {
    numlabels <- numsat / 1000000
    label.scale <- "Mb"
  } else if (numdigits < 6) {
    numlabels <- unit / 1000
    label.scale <- "kb"
  }
  
  subunits <- seq(0, maxnum, by = subunit)
  subunits <- subunits[!subunits %in% c(numsat, 0)]
  # return
  list(numsat, numlabels, label.scale, subunits)
}

#########################################################################################
#'@param sylist a tab-delimited file with seven columns with the header, required
#'              1. syri output file name
#'              2. Reference name
#'              3. file name of reference chr lengths (2 columsn: chr and length in bp)
#'              4. BED file name for reference highlights; set to NA if no inputs
#'              5. Query name
#'              6. file name of query chr lengths (2 columsn: chr and length in bp)
#'              7. BED file name for query highlights; set to NA if no inputs
#'                 1) chr 2) start 3) end 4) feature 5) color 6) proportion of the chr height
#'@param chr name of the select chromosome, required
#'@param add2existingplot logic value to specify whether drawing is added to an existing plot
#'@param min.syn.size minimum length of syntenic blocks to be plotted (10000)
#'                    either ref or qry block larger than min.syn.size will be plotted
#'@param min.others.size minimum length of other blocks to be plotted (50000)
#'                       either ref or qry block larger than min.syn.size will be plotted
#'@param genome.name.space space length used for labeling genome names in the proportion of plot width (0.1)
#'@param genome.label.col color for chromosome text ("palevioletred4")
#'@param chr.col chromosome background color ("gray80")
#'@param chr.height.prop height of each chromosome in the proportion of plot height (0.02)
#'@param main main text (NULL)
#'@param main.pos text coordinates (x, y) and a value (1, 2, 3, 4) to pass to "pos" to specify (below, left, above, right)
#'@param main.label.col color of the main text (palevioletred4)
#'@param band.legend.add logic value to indicate whether to add legend (TRUE)
#'@param band.legend.space.prop height of legend space in the proportion of plot height (0.1) 
#'@param band.legend.xpos start x-axis position for legends (0.5)
#'@param outpdf logic value to indicate whether a PDF will be output (FALSE)
#'@param outfile file name including path for the PDF output (NULL)
#'@param pdfwidth inch of width of output PDF figure (6)
#'@param pdfheight inch of height of output PDF figure (3.5)
#'@method sv types: synteny, duplication, translocation, inversion
#'                  synteny=SYN; duplication=DUP+INVDP; translocation=TRANS+INVTR; inversion=TDM
#'@author Sanzhen Liu, liu3zhen@gmail.com
#'@description 
#'@example 
#'@return
#'
mchrplot <- function(srlist, add2existingplot=F,
                   min.syn.size=10000, min.others.size=10000,
                   genome.name.space=0.1, genome.label.col="palevioletred4",
                   chr.col="gray80", chr.height.prop=0.02,
                   is.ygap.bw.comparison=F, ygap.bw.comparison=0.01,
                   main=NULL, main.pos=NULL, main.label.col="palevioletred4",
                   band.legend.add=T, band.legend.space.prop=0.1, band.legend.xpos=0.5,
                   high_chrcol="orange", low_chrcol="grey96",
                   outpdf=F, outfile=NULL, pdfwidth=6, pdfheight=4) {
  
  if (outpdf) {
    pdf(outfile, width=pdfwidth, height=pdfheight)
    stopifnot(!is.null(outfile))
  }
  
  # the canvas from 0 to 1 in x-axis
  xleft=0
  xright=1
  
  ### read data
  sr <- read.delim(srlist, stringsAsFactors=F, header=F)
  ncomparisons <- nrow(sr) # number of syri comparisons
  
  ### chromosome lengths and the maxial lengths from all 
  max.chrlen <- NULL
  for (i in 1:ncomparisons) {
    genome_name_ref <- sr[i, 2]
    genome_name_qry <- sr[i, 5]
    refchr <- sr[i, 8]
    qrychr <- sr[i, 9]
    chrlengths_ref <- read.delim(sr[i, 3], stringsAsFactors=F, header=F)
    chrlen_ref <- chrlengths_ref[chrlengths_ref[,1]==refchr, 2]
    chrlengths_qry <- read.delim(sr[i, 6], stringsAsFactors=F, header=F)
    chrlen_qry <- chrlengths_qry[chrlengths_qry[,1]==qrychr, 2]

    max.chrlen <- max(max.chrlen, chrlen_ref, chrlen_qry)
  }
  
  # canvasas (0, 1)
  par(mar=c(2, 0, 0, 0))
  
  # plot xrange and unit
  xaxtdata <- smart.axis(max.chrlen)
  x.unit <- max.chrlen / (xright - xleft)
  
  # coordinate conversion modules
  xpos.conversion <- function(pos) {
    pos / x.unit + xleft
  }
  
  #####################################################################################
  ### chromosome draw
  #####################################################################################
  chrdraw <- function(genome.name, chr, chrlen, ypos, height, chr.col,
                      text.add=T, text.yadj.value=0, text.col="palevioletred4") {
    rect(xpos.conversion(1), ypos-height/2,
         xpos.conversion(chrlen), ypos+height/2,
         col=chr.col, border=chr.col)
    chr.label <- paste(genome.name, chr)
    if (text.add) {
      text(0, ypos+text.yadj.value, labels=chr.label, pos=2, col=text.col)
    }
  }
  
  #####################################################################################
  #' syntenic link
  #####################################################################################
  syndraw <- function(alignment, ybottom, ytop, color="gray75") {
    # alignment is a vector with four column; first two from REF, 2nd two from QRY
    alignment <- as.numeric(as.character(alignment))
    alignment <- xpos.conversion(alignment)
	  polygon(alignment[c(1,3,4,2)], c(ybottom, ytop, ytop, ybottom), border=NA, col=color)
  }
  
  ###############################################################################
  ### module to draw bands
  ###############################################################################
  bandconnect <- function(qryregion, refregion, ytop, ybottom, inversion=F,
                          border=NA, bandcol="brown") {
    ### module
    transform_curve <- function(startp, endp, npoint=1000) {
      ### computer transform_curve coordinates
      midp <- (startp + endp) / 2
      beizer_value <- sqrt(0:(npoint-1)) / sqrt(npoint-1)  # sqrt as default
      
      curve.x1 <- seq(startp[1], midp[1], by=(midp[1] - startp[1])/(npoint-1))
      curve.y1 <- startp[2] - beizer_value * (startp[2] - midp[2])
      
      curve.x2 <- rev(seq(endp[1], midp[1], by=(midp[1] - endp[1])/(npoint-1)))
      curve.y2 <- rev(endp[2] - beizer_value * (endp[2] - midp[2]))
      
      curve.x <- c(curve.x1, curve.x2)
      curve.y <- c(curve.y1, curve.y2)
      list(x=curve.x, y=curve.y)
    }
    
    if (inversion) {
      qryregion <- rev(qryregion)
    }
    
    p1 <- transform_curve(c(qryregion[1], ytop), c(refregion[1], ybottom))
    p2 <- transform_curve(c(qryregion[2], ytop), c(refregion[2], ybottom))
    px <- c(p1$x, rev(p2$x))
    py <- c(p1$y, rev(p2$y))
    #lines(px, py)
    polygon(px, py, border=border, col=bandcol)
  }
  
  ###############################################################################
  ### module to draw bands for a certain sv type
  ###############################################################################
  type.band.dram <- function(syridata, svtype, min.len, bandcol, ybottom, ytop, inversion=F) {
    syri.type <- syridata[syridata$Type==svtype &
                            (abs(syridata[, 2] - syrichr[, 1]) >= min.len | 
                               abs(syrichr[, 4]- syrichr[, 3]) >= min.len), ]
    if (nrow(syri.type) > 0) {
      for (i in 1:nrow(syri.type)) {
        bandconnect(qryregion=as.numeric(xpos.conversion(syri.type[i, 3:4])),
                    refregion=as.numeric(xpos.conversion(syri.type[i, 1:2])),
                    ybottom=ybottom, ytop=ytop, inversion=inversion,
                    border=bandcol, bandcol=bandcol)
      }
    }
  }
  
  ###############################################################################
  ### module to highlight chromosome segments
  ###############################################################################
  seghighlight <- function(syridata, chr, svtype, genome=c("ref","qry"), min.len=5000, color, ypos, height) {
    # module to highlight regions of certain sv types on chromosome
    rectplot <- function(two_coordinates, ybottom, ytop, ...) {
      ### module to plot rect
      rect(xpos.conversion(two_coordinates[1]), ybottom,
           xpos.conversion(two_coordinates[2]), ytop, ...)
    }
    
    stopifnot(genome=="ref" | genome=="qry")
    if (genome == "ref") {
      col0 <- 1; col1 <- 2; col2 <- 3
    } else {
      col0 <- 6; col1 <- 7; col2 <- 8
    }
    # subset
    syridata <- syridata[syridata[, 11]==svtype & syridata[, col0]==chr, ]
    syridata[, col1] <- as.numeric(as.character(syridata[, col1]))
    syridata[, col2] <- as.numeric(as.character(syridata[, col2]))
    syri.highlight <- syridata[abs(syridata[, col2] - syridata[, col1]) >= min.len, ]
    if (nrow(syri.highlight) > 0) {
      apply(syri.highlight[, col1:col2], 1, rectplot, ybottom=ypos-height/2, ytop=ypos+height/2,
            col=color, border=NA)
    }
  }
  
  ###############################################################################
  ### module for chromosome highlights from an exterior input
  ###############################################################################
  chr.add.highlight <- function(bed, chr, ypos, chr.height) {
    features <- NULL
    feature.cols <- NULL
    feature.highlight <- read.delim(bed, header=F, comment.char="#", stringsAsFactors=F) # read bed file
    feature.highlight.chr <- feature.highlight[feature.highlight[,1]==chr, ]
    if (nrow(feature.highlight.chr) > 0) {
      for (i in 1:nrow(feature.highlight.chr)) {
        feature.highlight.pos <- xpos.conversion(as.numeric(feature.highlight.chr[i, 2:3]))
        feature.highlight.col <- feature.highlight.chr[i, 5]
        feature.highlight.chr.prop <- feature.highlight.chr[i, 6]
        feature.highlight.height <- chr.height * feature.highlight.chr.prop
        # draw
        rect(feature.highlight.pos[1], ypos-feature.highlight.height/2,
             feature.highlight.pos[2], ypos+feature.highlight.height/2,
             border=feature.highlight.col, col=feature.highlight.col, lwd=0.2)
        if (sum(features %in% feature.highlight.chr[i, 4]) == 0) {
          features <- c(features, feature.highlight.chr[i, 4])
          feature.cols <- c(feature.cols, feature.highlight.col)
        }
      }
    }
    
    col.out <- NULL
    if (!is.null(feature.cols)) {
      col.out <- data.frame(feature=features, col=feature.cols)
    }
    col.out
  }
  
  ########################################################################################
  ### plot preparation
  ########################################################################################
  ### y axis length
  ymax = ncomparisons / (1 - band.legend.space.prop)
  chr.height <- ymax * chr.height.prop # chromosome height

  if (is.ygap.bw.comparison) {
  	ygap.height <- ymax * ygap.bw.comparison
	ymax.adjust = ymax - ygap.height - chr.height
  } else {
  	ymax.adjust = ymax
  }
  
  ### color scheme
  band.types <- c("synteny", "duplication", "translocation", "inversion")
  band.cols <- c("gray75", "orchid4", "goldenrod3", "darkolivegreen4")
  names(band.cols) <- band.types
  
  ########################################################################################
  ### plot canvas and labels
  ########################################################################################
  genome.name.space.len <- (xleft - xright) * genome.name.space
  if (!add2existingplot) {
    plot(NULL, NULL, xlim=c(genome.name.space.len, xright), ylim=c(0, ymax.adjust), axes=F, xlab="", ylab="", main="")
  }
  
  ### xaxis
  xaxs.lines.at <- xpos.conversion(xaxtdata[[1]])
  axis(side=1, at=xaxs.lines.at, labels=xaxtdata[[2]])
  xaxs.lines.sub <- xpos.conversion(xaxtdata[[4]])
  
  #abline(v=xaxs.lines.at, col="lightsteelblue1", lwd=2)
  #abline(v=xaxs.lines.sub, col="gray95", lwd=1)
  axis_yvalue <- -0.03 * ymax.adjust
  
  for (tick in xaxs.lines.at) {
  	lines(c(tick, tick), c(axis_yvalue, ymax.adjust * 0.93),
	      col="lightsteelblue1", lwd=2)
  }
  
  for (sub.tick in xaxs.lines.sub) {
  	lines(c(sub.tick, sub.tick), c(axis_yvalue, ymax.adjust * 0.93),
	      col="lightsteelblue1", lwd=1)
  }

  if (is.null(main.pos)) {
    text(xleft, ymax.adjust * 0.96, pos=4, labels=main, cex=1.2, col=main.label.col, xpd=T)
  } else {
    text(main.pos[1], main.pos[2], pos=main.pos[3], labels=main, cex=1.2, col=main.label.col)
  }
  
  ### plot chr, highlights, and bands
  y.base <- 0
  previous.genome <- NULL
  for (i in 1:ncomparisons) {
    syri <- read.delim(sr[i, 1], header=F)
    ref.name <- sr[i, 2]
    refchr <- sr[i, 8]
	ref.ypos <- y.base
    qry.name <- sr[i, 5]
	qrychr <- sr[i, 9]
	
	chrlengths_ref <- read.delim(sr[i, 3], stringsAsFactors=F, header=F)
	chrlengths_qry <- read.delim(sr[i, 6], stringsAsFactors=F, header=F)

    if (is.ygap.bw.comparison) {
      qry.ypos <- y.base+1-ygap.height-chr.height
    } else {
      qry.ypos <- y.base+1
    }
    ### ref chromosome
    if (is.null(previous.genome) | is.ygap.bw.comparison) {
      ref.chrlen <- chrlengths_ref[chrlengths_ref[,1]==refchr, 2]
	  chrdraw(genome.name=ref.name, chr=refchr, chrlen=ref.chrlen, ypos=ref.ypos, chr.col=chr.col,
            height=chr.height, text.add=T, text.col=genome.label.col)
      ### ref chromosome highlights
      ref.bedfile <- sr[i, 4]
      if (!is.na(ref.bedfile)) {
        chr.add.highlight(bed=ref.bedfile, chr=refchr, ypos=ref.ypos, chr.height=chr.height)
      }
    } else if (!is.ygap.bw.comparison & ref.name != previous.genome) {
      stop("two neighboring SYRI results did not share a common genome\n")
    }
  
    ### qry chromosome
    qry.chrlen <- chrlengths_qry[chrlengths_qry[,1]==qrychr, 2]
	  chrdraw(genome.name=qry.name, chr=qrychr, chrlen=qry.chrlen, ypos=qry.ypos, chr.col=chr.col,
            height=chr.height, text.add=T, text.col=genome.label.col)
    
    ### ref chromosome highlights
    qry.bedfile <- sr[i, 7]
    if (!is.na(qry.bedfile)) {
      chr.add.highlight(bed=qry.bedfile, chr=qrychr, ypos=qry.ypos, chr.height=chr.height)
    }
    
    # memorize the query genome just plotted
    previous.genome <- qry.name
    
    ### extract SYRI data of a given chr
    syrichr <- syri[syri[, 1]==refchr & syri[, 6]==qrychr, c(2,3,7,8,9,10,11,12)]
    colnames(syrichr) <- c("refS", "refE", "qryS", "qryE", "ID", "ParentID", "Type", "Copy")
    for (i in 1:4) {
      syrichr[, i] <- as.numeric(syrichr[, i])
    }
    
    ### base to draw connection bands
    band.qrybase <- qry.ypos-chr.height/2
    band.refbase <- ref.ypos+chr.height/2
    
    ### syntenic block
    syn <- syrichr[syrichr$Type=="SYN" &
                  (abs(syrichr$refE - syrichr$refS) >= min.syn.size & 
                  abs(syrichr$qryE - syrichr$qryS) >= min.syn.size), ]
    stopifnot(nrow(syn)>0)
    
    # draw syn
    apply(syn[,1:4], 1, syndraw, ybottom=band.refbase, ytop=band.qrybase, color=band.cols["synteny"])
    
    ### duplication
    # type, minimum length, col
    type.band.dram(syridata=syrichr, svtype="DUP", min.len=min.others.size,
                   bandcol=band.cols["duplication"], ybottom=band.refbase, ytop=band.qrybase)
    
    ### inverted duplication
    type.band.dram(syridata=syrichr, svtype="INVDP", min.len=min.others.size,
                   bandcol=band.cols["duplication"], ybottom=band.refbase, ytop=band.qrybase)
    
    ### TRANS
    type.band.dram(syridata=syrichr, svtype="TRANS", min.len=min.others.size,
                   bandcol=band.cols["translocation"], ybottom=band.refbase, ytop=band.qrybase)
    
    ### INVTR
    type.band.dram(syridata=syrichr, svtype="INVTR", min.len=min.others.size,
                   bandcol=band.cols["translocation"], ybottom=band.refbase, ytop=band.qrybase)
    
    ### INV
    type.band.dram(syridata=syrichr, svtype="INV", min.len=min.others.size,
                   bandcol=band.cols["inversion"], ybottom=band.refbase, ytop=band.qrybase, inversion=T)
    
    y.base <- y.base + 1
  }
  
  ### legends
  if (band.legend.add) {
    legend(band.legend.xpos, ymax.adjust*1.03, legend=band.types, col=band.cols,
           lty=1, bty="n", lwd=3, ncol=2, xpd=T, cex=0.95)
  }
  
  ### pdf out
  if (outpdf) { dev.off() }
  
  ### return two values
  invisible(c(max.chrlen, x.unit, xleft))
}

