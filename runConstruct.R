#!/usr/bin/env Rscript

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
                                help="Input geoDist matrix file (tab-delimited)",
                                metavar="character"),
                    make_option(c("-c", "--coords"),
                                type="character",
                                default=NULL,
                                help="Input coordinates file (tab-delimited)"),
                    make_option(c("-w", "--wd"),
                                type="character",
                                default="./",
                                metavar="character",
                                help="Set working directory; default = ./"),
                    make_option(c("--prefix"),
                                type="character",
                                default=NULL,
                                help="Specify prefix for output files in conStruct analysis",
                                metavar="character"),
                    make_option(c("--nchains"),
                                type="integer",
                                default=1,
                                help="Specify number of independend MCMC chains to run",
                                metavar="integer"),
                    make_option(c("-i", "--niter"),
                                type="integer",
                                default=1000,
                                help="Specify length of MCMC chain to run",
                                metavar="integer"),
                    make_option(c("-K", "--K"),
                                type="integer",
                                default=NULL,
                                help="Specify K value (Number of popualations) to run",
                                metavar="integer"),
                    make_option(c("-a", "--afreq"),
                                type="character",
                                default="construct.out",
                                help="Specify allele frequency file to output",
                                metavar="character"))

opt_parser <- OptionParser(option_list=option_list,
                           description="Rscript to run conStruct")

opt <- parse_args(opt_parser)


#############################################################################

## FUNCTIONS ##

read.infile <- function(infile, header) {
  
  if(!is.null(infile)) {
    
    if (!file.exists(infile)) {
      stop("Error: Could not find required input file; aborting program")
    }
    
    file <- read.table(infile, 
                       header = header, 
                       sep = "\t",
                       stringsAsFactors = FALSE)
    return(file)
  }
  
  else {
    return(NULL)
  }
}

required.args <- function(arg, argString) {
  
  if(is.null(arg)) {
    
    stop(paste0("\n\nError: ", argString, " is a required argument\n"))
  }
  
}

#############################################################################

# Make sure dependencies are installed and then load them
if (!require("conStruct")) stop("Error: The required package conStruct is not installed")

# Tell user what the WD will be
print(paste("Setting working directory to", getwd(), sep="... "))

# Set working directory
setwd(opt$wd)

required.args(opt$prefix, "--prefix")
required.args(opt$str, "--str")
required.args(opt$coords, "--coords")
required.args(opt$K, "--K")
required.args(opt$geodist, "--geodist")

# Read input files specified on command-line
str.file <- read.infile(opt$str, FALSE)
popmap.file <- read.infile(opt$popmap, FALSE)
coords <- read.infile(opt$coords, TRUE)
geo.dist <- read.infile(opt$geodist, FALSE)

# Append K value to afreq filename
opt$afreq <- paste0(opt$afreq, "_", opt$K)

# Get population ID vector from STRUCTURE file
# For use with getting population allele frequency means
oneInd.str <- str.file[!duplicated(str.file$V1), ]
pop.index <- oneInd.str$V2

pop.ids <- geo.dist[,1]

# Remove first column that contains population IDs
geo.dist <- geo.dist[,2:ncol(geo.dist)]

# Convert geo.dist and coords to matrix class; required by conStruct
geo.distMat <- data.matrix(geo.dist)
coordMat <- data.matrix(coords)

# Exception handling if allele frequency RData file already exists
if (file.exists(paste0(opt$afreq, ".RData"))) {
    file.remove(paste0(opt$afreq, ".RData"))
  }

# Convert STRUCTURE file to conStruct allele frequency format
afreq <- structure2conStruct(infile = opt$str, 
                      onerowperind = FALSE, 
                      start.loci = 3, 
                      missing.datum = -9, 
                      outfile = opt$afreq)


# Get population means of allele frequencies if popmap file is given
if (!is.null(opt$popmap)){
  pop.data.matrix <- matrix(NA,nrow=nrow(coords),ncol=ncol(afreq))
  for(i in 1:nrow(pop.data.matrix)){
    pop.data.matrix[i,] <- colMeans(
      afreq[
        which(pop.index==i),,
        drop=FALSE
        ],na.rm=TRUE
    )
  }
}

# Print parameter settings
print("\nArguments provided: \n")
print(paste0("--str = ", opt$str))
print(paste0("--popmap = ", opt$popmap))
print(paste0("--geodist = ", opt$geodist))
print(paste0("--coords = ", opt$coords))
print(paste0("--wd = ", opt$wd))
print(paste0("--prefix = ", opt$prefix))
print(paste0("--nchains = ", opt$nchains))
print(paste0("--niter = ", opt$niter))
print(paste0("--K = ", opt$K))
print(paste0("--afreq = ", opt$afreq))

# Set prefixes for spatial (sp) and nonspatial (nsp) models
sp.prefix <- paste0(opt$prefix, "_spK", opt$K)
nsp.prefix <- paste0(opt$prefix, "_nspK", opt$K)

# Spatial model
construct.sp <- conStruct(spatial = TRUE,
                             K = opt$K,
                             freqs = pop.data.matrix,
                             geoDist = geo.distMat,
                             coords = coordMat,
                             prefix = sp.prefix,
                             n.chains = opt$nchains,
                             n.iter = opt$niter,
                             make.figs = TRUE,
                             save.files = TRUE)
  
# Non-spatial model
construct.nsp <- conStruct(spatial = FALSE,
                               K = opt$K,
                               freqs = pop.data.matrix,
                               geoDist = NULL,
                               coords = coordMat,
                               prefix = nsp.prefix,
                               n.chains = opt$nchains,
                               n.iter = opt$niter,
                               make.figs = TRUE,
                               save.files = TRUE)


if(!file.exists("environment.RData")) {
  save(pop.data.matrix, geo.distMat, coordMat, sp.prefix, nsp.prefix, pop.ids, opt, file = "environment.RData")
}


print(paste0("conStruct analysis for K = ", opt$K, " finished!"))
