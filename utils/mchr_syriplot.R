#!/usr/bin/env Rscript

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
#'              8. Reference chr name
#'              9. Query chr name
#'@param chr name of the select chromosome, required
#'@param min.syn.size minimum length of syntenic blocks to be plotted (10000)
#'                    either ref or qry block larger than min.syn.size will be plotted
#'@param min.others.size minimum length of other blocks to be plotted (50000)
#'                       either ref or qry block larger than min.syn.size will be plotted
#'@param max.invseg.ratio maximal length ratio of two inversion segments (3)
#'@param minINSDEL minimal INSDEL in syntenic regions (5000)
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
#'@param outfile file name including path for the PDF output (NULL)
#'@param pdfwidth inch of width of output PDF figure (6)
#'@param pdfheight inch of height of output PDF figure (3.5)
#'@method sv types: synteny, duplication, translocation, invrsion
#'                  synteny=SYN; duplication=DUP+INVDP; translocation=TRANS+INVTR; inversion=TDM
#'@author Sanzhen Liu, liu3zhen@gmail.com
#'@description 
#'@example 
#'@return
#'
mchrplot <- function(srlist, min.syn.size=10000, minINSDEL=5000,
                   min.others.size=10000, max.invseg.ratio=3,
                   genome.name.space=0.15, genome.label.col="palevioletred4",
                   chr.col="gray80", chr.height.prop=0.02,
                   ygap.bw.comparison=0.03,
                   main="", main.pos=NULL, main.label.col="palevioletred4",
                   band.legend.add=T, band.legend.space.prop=0.15, band.legend.xpos=0.4,
                   outfile=NULL, pdfwidth=6, pdfheight=4) {

  stopifnot(!is.null(outfile))
  pdf(outfile, width=pdfwidth, height=pdfheight)
  
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
  
  synINSdraw <- function(alignment, ybottom, ytop, linecol="coral1", fillcol="#FFC0CB4D") {
    # alignment is a vector with four column; first two from REF, 2nd two from QRY
    alignment <- as.numeric(as.character(alignment))
    alignment <- xpos.conversion(alignment)
	polygon(alignment[c(1,3,4,2)], c(ybottom, ytop, ytop, ybottom), border=linecol, col=fillcol)
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
  #type.band.dram <- function(syridata, svtype, min.len, bandcol, ybottom, ytop, inversion=F) {
   # syri.type <- syridata[syridata$Type==svtype &
   #                         (abs(syridata[, 2] - syrichr[, 1]) >= min.len | 
    #                           abs(syrichr[, 4]- syrichr[, 3]) >= min.len), ]
    #if (nrow(syri.type) > 0) {
     #for (i in 1:nrow(syri.type)) {
      #  bandconnect(qryregion=as.numeric(xpos.conversion(syri.type[i, 3:4])),
       #             refregion=as.numeric(xpos.conversion(syri.type[i, 1:2])),
        #            ybottom=ybottom, ytop=ytop, inversion=inversion,
         #           border=bandcol, bandcol=bandcol)
      #}
   # }
  #}
  type.band.dram <- function(syridata, svtype, min.len, bandcol, ybottom, ytop, inversion=F, segmaxratio=10000) {
    syri.type <- syridata[syridata$Type==svtype &
                            (abs(syridata[, 2] - syrichr[, 1]) >= min.len | 
                               abs(syrichr[, 4]- syrichr[, 3]) >= min.len), ]

	if (nrow(syri.type) > 0) {
      for (i in 1:nrow(syri.type)) {
        seg_size_ratio = abs(syri.type[i, 2] - syri.type[i, 1]) / abs(syri.type[i, 4]- syri.type[i, 3])
		if (seg_size_ratio < 1) {
			seg_size_ratio <- 1 / seg_size_ratio;
		}
		if (seg_size_ratio <= segmaxratio) {
	      bandconnect(qryregion=as.numeric(xpos.conversion(syri.type[i, 3:4])),
                      refregion=as.numeric(xpos.conversion(syri.type[i, 1:2])),
                      ybottom=ybottom, ytop=ytop, inversion=inversion,
                      border=bandcol, bandcol=bandcol)
        }
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
  ### y height
  main.space.prop <- 0.1 # if no band legend is added
  main.space <- 0
  if (band.legend.add) {
    ymax = ncomparisons / (1 - band.legend.space.prop)
  } else if (main != "") {
    ymax = ncomparisons / (1 - main.space.prop)
	main.space <- ymax * main.space.prop
  } else {
    ymax = ncomparisons
  }

  chr.height <- ymax * chr.height.prop # chromosome height
  ygap.height <- ymax * ygap.bw.comparison
  ymax.adjusted = ymax - ygap.height - chr.height
  
  ### color scheme
  band.types <- c("synteny", "dup", "translocation", "inv")
  band.cols <- c("gray75", "orchid4", "goldenrod3", "darkolivegreen4")
  names(band.cols) <- band.types
  
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
  axis(side=1, at=xaxs.lines.at, labels=xaxt.labels)
  xaxs.lines.sub <- xpos.conversion(xaxtdata[[4]])

  axis_ylow <- -0.02 * ymax.adjusted

  if (band.legend.add) {
    axis_ytop <- ymax.adjusted - ymax * band.legend.space.prop + chr.height * 0.9
  } else {
    axis_ytop <- ymax.adjusted - main.space + chr.height * 0.9
  }

  for (tick in xaxs.lines.at) {
    lines(c(tick, tick), c(axis_ylow, axis_ytop), col="lightsteelblue2", lwd=2, xpd=T)
  }
  
  for (sub.tick in xaxs.lines.sub) {
    lines(c(sub.tick, sub.tick), c(axis_ylow, axis_ytop), col="lightsteelblue1", lwd=0.8, xpd=T)
  }

  ### main text
  if (main != "") {
    if (is.null(main.pos)) {
      if (band.legend.add) {
        text(xleft, ymax.adjusted - ymax * band.legend.space.prop * 0.5,
	         pos=4, labels=main, cex=1.2, col=main.label.col, xpd=T)
      } else {
	    text(xleft, ymax.adjusted - main.space * 0.4, pos=4, labels=main, cex=1.2, col=main.label.col, xpd=T)
      }
    } else {
      text(main.pos[1], main.pos[2], pos=main.pos[3], labels=main, cex=1.2, col=main.label.col)
    }
  }

  ### plot chr, highlights, and bands
  y.base <- 0
  for (i in 1:ncomparisons) {
    syri <- read.delim(sr[i, 1], header=F)
    ref.name <- sr[i, 2]
    refchr <- sr[i, 8]
	ref.ypos <- y.base
    qry.name <- sr[i, 5]
	qrychr <- sr[i, 9]
	
	chrlengths_ref <- read.delim(sr[i, 3], stringsAsFactors=F, header=F)
	chrlengths_qry <- read.delim(sr[i, 6], stringsAsFactors=F, header=F)

    qry.ypos <- y.base+1-ygap.height-chr.height

    ### ref chromosome
    ref.chrlen <- chrlengths_ref[chrlengths_ref[,1]==refchr, 2]
	chrdraw(genome.name=ref.name, chr=refchr, chrlen=ref.chrlen, ypos=ref.ypos, chr.col=chr.col,
            height=chr.height, text.add=T, text.col=genome.label.col)
    ### ref chromosome highlights
    ref.bedfile <- sr[i, 4]
    if (!is.na(ref.bedfile)) {
      chr.add.highlight(bed=ref.bedfile, chr=refchr, ypos=ref.ypos, chr.height=chr.height)
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
   
    # draw INSDEL on syn
    parentid <- syrichr$ParentID
    parentid <- gsub("[0-9]+$", "", parentid)
    synINS <- syrichr[(syrichr$Type=="INS" | syrichr$Type=="DEL") & parentid=="SYN" &
                   ((abs(syrichr$refE - syrichr$refS) >= minINSDEL ) |
				   (abs(syrichr$qryE - syrichr$qryS) >= minINSDEL)), ]
   
    if (nrow(synINS) > 0) {
  	  apply(synINS[,1:4], 1, synINSdraw, ybottom=band.refbase + 0.01, ytop=band.qrybase - 0.01)
    }
 
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
    type.band.dram(syridata=syrichr, svtype="INV", min.len=min.others.size, segmaxratio=max.invseg.ratio,
                   bandcol=band.cols["inversion"], ybottom=band.refbase, ytop=band.qrybase, inversion=T)
    
    y.base <- y.base + 1
  }
  
  ### legends
  if (band.legend.add) {
    legend(band.legend.xpos, ymax.adjusted*1.02, legend=band.types, col=band.cols,
           lty=1, bty="n", lwd=3, ncol=2, xpd=T, cex=0.9)
  }
  
  ### close pdf
  dev.off()
}

### run
library(optparse)

### input options
option_list <- list(
  # Required
  make_option("--srlist", type="character", help="SR list file (required)
                1. syri output file name
                2. Ref name (e.g., species name)
                3. file name of reference chr lengths (2 columsn: chr and length in bp)
                4. BED file name for reference highlights; set to NA if no inputs
                   1) chr 2) start 3) end 4) feature 5) color 6) proportion of the chr height
                5. Query name (e.g., species name)
                6. file name of query chr lengths (2 columsn: chr and length in bp)
                7. BED file name for query highlights; set to NA if no inputs
                8. Ref chr name
                9. Query chr name"),
 
  # Plot behavior
  make_option("--min_syn_size", type="integer", default=10000,
              help="Minimum synteny size [default %default]"),
  make_option("--minINSDEL", type="integer", default=5000,
              help="Minimum insertion or deletion in syntenic regions [default %default]"),
  make_option("--min_others_size", type="double", default=10000,
              help="Minimum other SV size [default %default]"),
  make_option("--max_invseg_ratio", type="double", default=3,
              help="Maximal length ratio of two inversion segments [default %default]"),

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
  make_option("--main_pos", type="character", default=NULL,
              help="Main position as comma-separated x,y,pos (e.g. 0.1,3,4)"),
  make_option("--main_label_col", type="character", default="palevioletred4",
              help="Main label color [default %default]"),

  # Legend
  make_option("--band_legend_add", type="logical", default=TRUE,
              help="Add legend [default %default]"),
  make_option("--band_legend_space_prop", type="double", default=0.15,
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

if (is.null(opt$srlist)) {
  stop("ERROR: --srlist is required")
}

################################################################################
# Convert main_pos if provided
################################################################################
main_pos_parsed <- NULL
if (!is.null(opt$main_pos)) {
  main_pos_parsed <- as.numeric(strsplit(opt$main_pos, ",")[[1]])
}

################################################################################
# Run function
################################################################################
mchrplot(
  srlist = opt$srlist,
  min.syn.size = opt$min_syn_size,
  minINSDEL = opt$minINSDEL,
  min.others.size = opt$min_others_size,
  max.invseg.ratio=opt$max_invseg_ratio,
  genome.name.space = opt$genome_name_space,
  genome.label.col = opt$genome_label_col,
  chr.col = opt$chr_col,
  chr.height.prop = opt$chr_height_prop,
  ygap.bw.comparison = opt$ygap_bw_comparison,
  main = opt$main,
  main.pos = main_pos_parsed,
  main.label.col = opt$main_label_col,
  band.legend.add = opt$band_legend_add,
  band.legend.space.prop = opt$band_legend_space_prop,
  band.legend.xpos = opt$band_legend_xpos,
  outfile = opt$outfile,
  pdfwidth = opt$pdfwidth,
  pdfheight = opt$pdfheight
)

