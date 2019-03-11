#!/usr/bin/env Rscript

# Make sure dependencies are installed and then load them
if (!require("conStruct")) stop("Error: The required package conStruct is not installed")
if (!require("optparse")) stop("Error: The required package optparse is not installed")

# Set command-line arguments
option_list <- list(make_option(c("-s", "--str"), 
                                type="character", 
                                default=NULL, 
                                help="Input structure filename; default = NULL", 
                                metavar="character"),
                    make_option(c("-p", "--popmap"),
                                 type="character",
                                 default=NULL,
                                 help="Input popmap filename; default = NULL",
                                 metavar="character"),
                    make_option(c("-g", "--geodist"),
                                type="character",
                                default=NULL,
                                help="Input geoDist matrix file (as CSV)",
                                metavar="character"),
                    make_option(c("-c", "--coords"),
                                type="character",
                                default=NULL,
                                help="Input coordinates file (as CSV)"),
                    make_option(c("-w", "--wd"),
                                type="character",
                                default="./",
                                metavar="character",
                                help="Set working directory; default = ./"))

opt_parser <- OptionParser(option_list=option_list,
                           description="Rscript to run conStruct")

opt <- parse_args(opt_parser)


#############################################################################

## FUNCTIONS ##

read.infile <- function(infile) {
  
  if(!is.null(infile)) {
    
    if (!file.exists(infile)) {
      stop("Error: Could not find required input file; aborting program")
    }
    
    file <- read.table(infile, 
                       header = FALSE, 
                       sep = "\t",
                       stringsAsFactors = FALSE)
    return(file)
  }
  
  else {
    return(NULL)
  }
}

#############################################################################

# Tell user what the WD will be
print(paste("Setting working directory to", getwd(), sep="... "))

# Set working directory
setwd(opt$wd)

# Read input files specified on command-line
str.file <- read.infile(opt$str)
popmap.file <- read.infile(opt$popmap)
geo.dist <- read.infile(opt$geodist)
#coords <- read.infile(opt$coords)

str.file <- read.csv("BOX_filteredPops_FINAL.str", 
         sep="\t", 
         header=FALSE,
         stringsAsFactors = FALSE)


pop.index <- str.file$V2

popmap.file <- read.csv("BOX_popmap_FINAL4construct.txt",
                        sep="\t",
                        header = FALSE,
                        stringsAsFactors = FALSE)
geo.dist <- read.csv("BOX_geodist.csv", 
                     header = FALSE,
                     stringsAsFactors = FALSE)

coords <- read.csv("coords.out.csv",
                   header = TRUE,
                   stringsAsFactors = FALSE)

# Formatting data
#vignette(topic="format-data",package="conStruct")

afreq <- structure2conStruct(infile = opt$str, 
                      onerowperind = FALSE, 
                      start.loci = 3, 
                      missing.datum = -9, 
                      outfile = "construct.out")

afreq <- structure2conStruct(infile = "BOX_filteredPops_FINAL.str", 
                             onerowperind = FALSE, 
                             start.loci = 3, 
                             missing.datum = -9, 
                             outfile = "construct.out")

pop.data.matrix <- matrix(NA,nrow=nrow(coords),ncol=ncol(afreq))
for(i in 1:nrow(pop.data.matrix)){
  pop.data.matrix[i,] <- colMeans(
    afreq[
      which(pop.index==i),,
      drop=FALSE
      ],na.rm=TRUE
  )
}

print(nrow())
print(pop.data.matrix)
#geodist <- conStruct.data$geoDist
#afreq <- conStruct.data$allele.frequencies
#coods <- conStruct.data$coords

#pop.data.matrix <- matrix(NA, nrow=4, ncol=ncol(afreq))


