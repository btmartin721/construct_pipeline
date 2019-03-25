#!/usr/bin/env Rscript

if (!require("optparse")) stop("Error: The required package optparse is not installed")

# Set command-line arguments
option_list <- list(make_option(c("-s", "--str"), 
                                type="character", 
                                default=NULL, 
                                help="Required; Input structure filename; default = NULL", 
                                metavar="character"),
                    make_option(c("-p", "--popmap"),
                                 type="character",
                                 default=NULL,
                                 help="Required; Input popmap filename; default = NULL",
                                 metavar="character"),
                    make_option(c("-g", "--geodist"),
                                type="character",
                                default=NULL,
                                help="Required; Input geoDist matrix file (tab-delimited)",
                                metavar="character"),
                    make_option(c("-c", "--coords"),
                                type="character",
                                default=NULL,
                                help="Required; Input coordinates file (tab-delimited)"),
                    make_option(c("--prefix"),
                                type="character",
                                default=NULL,
                                help="Required; Specify prefix for output files in conStruct analysis",
                                metavar="character"),
                    make_option(c("-K", "--K"),
                                type="integer",
                                default=NULL,
                                help="Required; Specify K value (Number of popualations) to run",
                                metavar="integer"),
                    make_option(c("--spatial"),
                                action="store_true",
                                default=FALSE,
                                help="Boolean; Toggle on spatial model; --spatial and/or --nonspatial are required; default = FALSE"),
                    make_option(c("--nonspatial"),
                                action="store_true",
                                default=FALSE,
                                help="Boolean; Toggle on nonspatial model; --spatial and/or --nonspatial are required; default = FALSE"),
                    make_option(c("-w", "--wd"),
                                type="character",
                                default="./",
                                metavar="character",
                                help="Optional; Set working directory; default = ./"),
                    make_option(c("--nchains"),
                                type="integer",
                                default=1,
                                help="Optional; Specify number of independend MCMC chains to run",
                                metavar="integer"),
                    make_option(c("-i", "--niter"),
                                type="integer",
                                default=1000,
                                help="Optional; Specify length of MCMC chain to run",
                                metavar="integer"),
                    make_option(c("-a", "--afreq"),
                                type="character",
                                default="construct.out",
                                help="Optional; Specify allele frequency file to output",
                                metavar="character"),
                    make_option(c("-o", "--outdir"),
                                type="character",
                                default="./",
                                help="Optional; Specify directory to write output files",
                                metavar="character"),
                    make_option(c("-r", "--onerowperind"),
                                action="store_true",
                                default=FALSE,
                                help="Optional, Boolean; If toggled, specifies only one row per individual in structure file; default = FALSE"),
                    make_option(c("--data_column"),
                                type="integer",
                                default=3,
                                help="Optional; Specify column index (1-based) for first data column in structure file; default = 3"),
                    make_option(c("-m", "--missingval"),
                                type="integer",
                                default=-9,
                                help="Optional; Specify missing data value in structure file; default = -9"),
                    make_option(c("-D", "--adapt_delta"),
                                type="numeric",
                                default=0.8,
                                help="Optional; Set adapt_delta value if you are having mixing problems",
                                metavar="numeric"),
                    make_option(c("-d", "--max_treedepth"),
                                type="numeric",
                                default=10,
                                help="Optional; Set max_treedepth",
                                metavar="numeric"))
                    
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

# Set working directory
setwd(opt$wd)

# Tell user what the WD will be
print(paste("Setting working directory to", getwd(), sep="... "))

required.args(opt$prefix, "--prefix")
required.args(opt$str, "--str")
required.args(opt$coords, "--coords")
required.args(opt$K, "--K")
required.args(opt$geodist, "--geodist")

if (!opt$spatial && !opt$nonspatial) {
  stop("Error: One or both of the following boolean options must be toggled: --spatial and/or --nonspatial")
}

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
                      onerowperind = opt$onerowperind, 
                      start.loci = opt$data_column, 
                      missing.datum = opt$missingval, 
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
print(paste0("--onerowperind = ", opt$onerowperind))
print(paste0("--data_column = ", opt$data_column))
print(paste0("--missingval = ", opt$missingval))
print(paste0("--spatial = ", opt$spatial))
print(paste0("--nonspatial = ", opt$nonspatial))
print(paste0("--adapt_delta = ", opt$adapt_delta))
print(paste0("--max_treedepth = ", opt$max_treedepth))

# Set prefixes for spatial (sp) and nonspatial (nsp) models
sp.prefix <- paste0(opt$prefix, "_spK", opt$K)
nsp.prefix <- paste0(opt$prefix, "_nspK", opt$K)

dir.create(opt$outdir)
setwd(opt$outdir)

if (opt$K == 1) {
  
  if (file.exists(paste0(opt$prefix, "_environment.RData"))) {
    print(paste0("\n\nWarning: The file ", opt$prefix, "_environment.RData already exists; renaming it to ", opt$prefix, "_tempenv.RData\n"))
    file.rename(paste0(opt$prefix, "_environment.RData"), paste0(opt$prefix, "tempenv.RData"))
  }
}

if (opt$spatial){
  
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
                             save.files = TRUE,
                             adapt_delta = opt$adapt_delta,
                             max_treedepth = opt$max_treedepth)
}

if (opt$nonspatial) {
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
                               save.files = TRUE,
                               adapt_delta = opt$adapt_delta,
                               max_treedepth = opt$max_treedepth)
}

if (opt$K == 1 && opt$spatial) {
  save(pop.data.matrix, 
       geo.distMat, 
       coordMat, 
       sp.prefix, 
       nsp.prefix, 
       pop.ids, 
       file = paste0(opt$prefix, "_environment.RData"))
  saveRDS(opt, file = "arguments.RDS")
}



print(paste0("conStruct analysis for K = ", opt$K, " finished!"))
