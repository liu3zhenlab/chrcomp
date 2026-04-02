#!/ucomp/bin/env Rscript

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
mchr_nucmerplot <- function(complist, 
                            min.aln=10000, min.identity=80,
							genome.name.space=0.15, genome.label.col="palevioletred4",
							chr.col="gray80", chr.height.prop=0.02,
							ygap.bw.comparison=0.03,
							main="", main.label.col="palevioletred4",
							inv.border.col="plum4",
							bandcol="burlywood3", col.legend.space.prop=0.15, col.legend.xpos=0.4,
							outfile=NULL, pdfwidth=6, pdfheight=5) {

  stopifnot(!is.null(outfile))
  pdf(outfile, width=pdfwidth, height=pdfheight)
  
  # the canvas from 0 to 1 in x-axis
  xleft=0
  xright=1
  
  ### read data
  comp <- read.delim(complist, stringsAsFactors=F, header=F)
  ncomparisons <- nrow(comp) # number of alignment comparisons
  
  ### chromosome lengths and the maxial lengths from all 
  ### range of alignment identity
  max.chrlen <- NULL
  identity_high <- NULL
  identity_low <- NULL
  for (i in 1:ncomparisons) {
    genome_name_ref <- comp[i, 2]
    genome_name_qry <- comp[i, 5]
    refchr <- comp[i, 8]
    qrychr <- comp[i, 9]
   
    # chr lengths 
	chrlengths.ref <- read.delim(comp[i, 3], stringsAsFactors=F, header=F)
	chrlen.ref <- chrlengths.ref[chrlengths.ref[,1]==refchr, 2]
	chrlengths.qry <- read.delim(comp[i, 6], stringsAsFactors=F, header=F)
	chrlen.qry <- chrlengths.qry[chrlengths.qry[,1]==qrychr, 2]
	
	max.chrlen <- max(max.chrlen, chrlen.ref, chrlen.qry)

	nuc_aln <- read.delim(comp[i, 1], stringsAsFactors=F)
    # filter alignments
	nuc_aln <- nuc_aln[(nuc_aln$rmatch >= min.aln | nuc_aln$qmatch >= min.aln) & (nuc_aln$identity >= min.identity), ]
	if (nrow(nuc_aln)>0) {
	  identity_high <- max(identity_high, max(nuc_aln$identity))
	  identity_low <- min(identity_low, min(nuc_aln$identity))
    }
  }
  
  ### check extracted data
  stopifnot(!is.null(max.chrlen))

  if (is.null(identity_high) | is.null(identity_low)) {
  	identity_high <- 100
	identity_low <- min.identity
  } else if (identity_high == identity_low) {
  	identity_low <- identity_high - 5
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
         col=chr.col, border=NA)
    chr.label <- paste(genome.name, chr)
    if (text.add) {
      text(0, ypos+text.yadj.value, labels=chr.label, pos=2, col=text.col, cex=0.9)
    }
  }
  
  ###############################################################################
  ### module to draw bands
  ###############################################################################
  bandconnect <- function(qryregion, refregion, ytop, ybottom,
                          inv.border.col="plum4", bandcol="grey") {
    
	qryregion <- as.numeric(as.character(qryregion))
	refregion <- as.numeric(as.character(refregion))
    
	### deal with inversion
	border=NA
    if ((qryregion[2] > qryregion[1] & refregion[2] < refregion[1]) |
	    (qryregion[2] < qryregion[1] & refregion[2] > refregion[1])) {
		# inversion
		border <- inv.border.col # draw border
	}

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
    
    p1 <- transform_curve(c(qryregion[1], ytop), c(refregion[1], ybottom))
    p2 <- transform_curve(c(qryregion[2], ytop), c(refregion[2], ybottom))
    px <- c(p1$x, rev(p2$x))
    py <- c(p1$y, rev(p2$y))
    polygon(px, py, border=border, lwd=0.1, col=bandcol)
  }
  
  ###########################################################
  #' check if a color is valid
  ###########################################################
  isColor <- function(col) {
    sapply(col, function(incol) {
      tryCatch(is.matrix(col2rgb(incol)), 
               error = function(e) FALSE)
    })
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
  ### y height
  ymax = ncomparisons / (1 - col.legend.space.prop)
  chr.height <- ymax * chr.height.prop # chromosome height
  ygap.height <- ymax * ygap.bw.comparison
  ymax.adjusted = ymax - ygap.height - chr.height
  
  # connection colors
  high_identity_col <- bandcol
  low_identity_col <- "grey96"
  colfunc <- colorRampPalette(c(low_identity_col, high_identity_col))
  color_num <- 20
  allcols <- colfunc(color_num)

  # identity
  if (identity_high - identity_low < 5) {
    identity_low <- identity_high - 5
  }

  identity20 <- seq(identity_low, identity_high, by=(identity_high - identity_low)/(color_num -1))

  ########################################################################################
  ### plot canvas and labels
  ########################################################################################
  genome.name.space.len <- (xleft - xright) * genome.name.space
    
  plot(NULL, NULL, xlim=c(genome.name.space.len, xright), ylim=c(0, ymax.adjusted), axes=F, xlab="", ylab="", main="")
  
  ### xaxis
  xaxs.lines.at <- xpos.conversion(xaxtdata[[1]])
  xaxt.labels <- xaxtdata[[2]]
  # add scale to the last element
  xaxt.labels[length(xaxt.labels)] <- paste(xaxt.labels[length(xaxt.labels)], xaxtdata[[3]])
  par(mgp = c(3, 0.5, 0))
  axis(side=1, at=xaxs.lines.at, labels=xaxt.labels)
  xaxs.lines.sub <- xpos.conversion(xaxtdata[[4]])

  axis_ylow <- -0.02 * ymax.adjusted

  axis_ytop <- ymax.adjusted - ymax * col.legend.space.prop + chr.height * 0.9

  for (tick in xaxs.lines.at) {
    lines(c(tick, tick), c(axis_ylow, axis_ytop), col="lightsteelblue2", lwd=2, xpd=T)
  }
  
  for (sub.tick in xaxs.lines.sub) {
    lines(c(sub.tick, sub.tick), c(axis_ylow, axis_ytop), col="lightsteelblue1", lwd=0.8, xpd=T)
  }

  ### main text
  if (main != "") {
    text(xleft, ymax.adjusted - ymax * col.legend.space.prop * 0.5,
	     pos=4, labels=main, cex=1.2, col=main.label.col, xpd=T)
  }

  ### plot chr, highlights, and bands
  y.base <- 0
  for (i in 1:ncomparisons) {
    ref.name <- comp[i, 2]
    refchr <- comp[i, 8]
	ref.ypos <- y.base
    qry.name <- comp[i, 5]
	qrychr <- comp[i, 9]
	
    qry.ypos <- y.base+1-ygap.height-chr.height

    ### ref chromosome
	chrlengths.ref <- read.delim(comp[i, 3], stringsAsFactors=F, header=F)
	chrlen.ref <- chrlengths.ref[chrlengths.ref[,1]==refchr, 2]
	chrdraw(genome.name=ref.name, chr=refchr, chrlen=chrlen.ref, ypos=ref.ypos, chr.col=chr.col,
            height=chr.height, text.add=T, text.col=genome.label.col)
    
	### ref chromosome highlights
    ref.bedfile <- comp[i, 4]
    if (!is.na(ref.bedfile)) {
      chr.add.highlight(bed=ref.bedfile, chr=refchr, ypos=ref.ypos, chr.height=chr.height)
    }
  
    ### qry chromosome
	chrlengths.qry <- read.delim(comp[i, 6], stringsAsFactors=F, header=F)
	chrlen.qry <- chrlengths.qry[chrlengths.qry[,1]==qrychr, 2]
	chrdraw(genome.name=qry.name, chr=qrychr, chrlen=chrlen.qry, ypos=qry.ypos, chr.col=chr.col,
            height=chr.height, text.add=T, text.col=genome.label.col)
    
    ### ref chromosome highlights
    qry.bedfile <- comp[i, 7]
    if (!is.na(qry.bedfile)) {
      chr.add.highlight(bed=qry.bedfile, chr=qrychr, ypos=qry.ypos, chr.height=chr.height)
    }
    
    ### base to draw connection bands
    band.qrybase <- qry.ypos-chr.height*0.6
    band.refbase <- ref.ypos+chr.height*0.6
    
    # connections
    nuc_aln <- read.delim(comp[i, 1])
	nuc_aln <- nuc_aln[(nuc_aln$rmatch >= min.aln | nuc_aln$qmatch >= min.aln) & (nuc_aln$identity >= min.identity), ]
    if (nrow(nuc_aln)>0) {
	  for (i in 1:nrow(nuc_aln)) {
        band_gradient_col <- allcols[which.min(abs(identity20 - nuc_aln[i, "identity"]))]
        qryregion_input <- xpos.conversion(as.numeric(as.character(nuc_aln[i, 3:4])))
	    refregion_input <- xpos.conversion(as.numeric(as.character(nuc_aln[i, 1:2])))
	    if (inv.border.col == "skip") {
		  inv.border.col <- NA
		}
		bandconnect(qryregion=qryregion_input, ytop=band.qrybase,
                  refregion=refregion_input, ybottom=band.refbase,
                  inv.border.col=inv.border.col, bandcol=band_gradient_col)
      }
    }
    # go to the next
    y.base <- y.base + 1
  }
  
  ### identity legends
  bar_ypos <- y.base + ymax * col.legend.space.prop * 0.4
  max.xsize <- 1
  barlen <- max.xsize / 4
  barstep <- barlen / color_num
  barheight <- ymax * col.legend.space.prop * 0.2
  text(max.xsize - barlen/2, bar_ypos, labels="identity", cex=0.8, pos=1, xpd=T)  ### plot LD name
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
         border=NA, col=allcols[i], xpd=T)
    col_barpos <- col_barpos + barstep
  }
 
  ### close pdf
  dev.off()
}

### run
library(optparse)

### input options
option_list <- list(
  # Required
  make_option("--complist", type="character", help="comparison list file (required)
                1. nucmer output file name, including the path
                2. ref name (e.g., species name)
                3. file name of ref chr lengths (2 columsn: chr and length in bp)
				4. BED file name for reference highlights; set to NA if no inputs
                   1) chr 2) start 3) end 4) feature 5) color 6) proportion of the chr height
                5. query name (e.g., species name)
                6. file name of query chr lengths (2 columsn: chr and length in bp)
				7. BED file name for query highlights; set to NA if no inputs
                8. ref chr name
                9. query chr name"),
 
  # alignment filter
  make_option("--min_aln", type="integer", default=10000,
              help="Minimum alignment length [default %default]"),
  make_option("--min_identity", type="double", default=80,
              help="Minimum identity of each alignment [default %default]"),

  # Genome label settings
  make_option("--genome_name_space", type="double", default=0.15,
              help="Genome name space proportion [default %default]"),
  make_option("--genome_label_col", type="character", default="palevioletred4",
              help="Genome label color [default %default]"),

  # Chromosome appearance
  make_option("--chr_col", type="character", default="gray80",
              help="Chromosome color [default %default]"),
  make_option("--chr_height_prop", type="double", default=0.02,
              help="Chromosome height proportion [default %default]"),

  # Y-gap control
  make_option("--ygap_bw_comparison", type="double", default=0.03,
              help="Gap proportion between comparisons [default %default]"),

  # Main title
  make_option("--main", type="character", default="", help="Main title text"),
  make_option("--main_label_col", type="character", default="palevioletred4",
              help="Main label color [default %default]"),

  # Legend
  make_option("--bandcol", type="character", default="burlywood3",
              help="color for alignmnet connection bands [default %default]"),
  make_option("--inv_border_col", type="character", default="plum4",
			  help="color of the border for alignmnet connection bands when the alignment is inverted:
			        skip - no border or a color [default %default]"),
  make_option("--col_legend_space_prop", type="double", default=0.15,
              help="Legend space proportion [default %default]"),
  make_option("--band_legend_xpos", type="double", default=0.4,
              help="Legend x position [default %default]"),

  # PDF output
  make_option("--outfile", type="character", default=NULL,
              help="Output PDF filename"),
  make_option("--pdfwidth", type="double", default=6,
              help="PDF width [default %default]"),
  make_option("--pdfheight", type="double", default=4,
              help="PDF height [default %default]")
)

opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$complist)) {
  stop("ERROR: --complist is required")
}

################################################################################
# Run function
################################################################################
mchr_nucmerplot(
  complist = opt$complist,
  min.aln = opt$min_aln,
  min.identity = opt$min_identity,
  genome.name.space = opt$genome_name_space,
  genome.label.col = opt$genome_label_col,
  chr.col = opt$chr_col,
  chr.height.prop = opt$chr_height_prop,
  ygap.bw.comparison = opt$ygap_bw_comparison,
  main = opt$main,
  main.label.col = opt$main_label_col,
  inv.border.col=opt$inv_border_col,
  bandcol=opt$bandcol,
  col.legend.space.prop = opt$col_legend_space_prop,
  col.legend.xpos = opt$col_legend_xpos,
  outfile = opt$outfile,
  pdfwidth = opt$pdfwidth,
  pdfheight = opt$pdfheight
)

