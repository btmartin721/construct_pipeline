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

# Formatting data
#vignette(topic="format-data",package="conStruct")

geodist <- conStruct.data$geoDist
afreq <- conStruct.data$allele.frequencies
coods <- conStruct.data$coords

pop.data.matrix <- matrix(NA, nrow=4, ncol=ncol(afreq))


