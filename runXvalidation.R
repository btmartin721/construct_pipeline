#!/usr/bin/env Rscript

if (!require("optparse")) stop("Error: The required package optparse is not installed")
library("optparse")

# Set command-line arguments
option_list <- list(make_option(c("-K", "--maxK"),
                                type="integer",
                                default=NULL,
                                help="Required; Specify maximum K value (upper case K)",
                                metavar="integer"),
                    make_option(c("-r", "--nreps"),
                                type="integer",
                                default=NULL,
                                help="Required; Specify number of cross-validation replicates",
                                metavar="integer"),
                    make_option(c("-n", "--nodes"),
                                type="integer",
                                default=NULL,
                                help="Required; Specify number of CPU cores for parallelization",
                                metavar="integer"),
                    make_option(c("-t", "--trainProp"),
                                type="numeric",
                                default=0.9,
                                help="Optional; Specify training proportion for cross-validation",
                                metavar="numeric"),
                    make_option(c("-f", "--saveFiles"),
                                action="store_true",
                                default = FALSE,
                                help="Optional, Boolean; Toggle save Robj files from cross-validation; will save lots of files; default=FALSE"),
                    make_option(c("-F", "--saveFigs"),
                                action="store_false",
                                default=TRUE,
                                help="Optional, Boolean; Don't save figures from cross-validation; default=TRUE"),
                    make_option(c("-o", "--outdir"),
                                type="character",
                                default="./xvalResults",
                                help="Optional, Specify directory for output files; will be created if doesn't exist",
                                metavar="character"),
                    make_option(c("-s", "--str"), 
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
                    make_option(c("-w", "--wd"),
                                type="character",
                                default="./",
                                metavar="character"),
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
                    make_option(c("-R", "--onerowperind"),
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
                           description="Rscript to run cross-validation for conStruct")

opt2 <- parse_args(opt_parser)

######################################################################

### FUNCTIONS ###

# Stops if argument is NULL
required.args <- function(arg, argString) {
  
  if(is.null(arg)) {
    
    stop(paste0("\n\nError: ", argString, " is a required argument\n"))
  }
  
}

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


#####################################################################3
# Require these arguments or stop
required.args(opt2$maxK, "--maxK")
required.args(opt2$nreps, "--nreps")
required.args(opt2$nodes, "--nodes")
required.args(opt2$wd, "--wd")
required.args(opt2$prefix, "--prefix")
required.args(opt2$str, "--str")
required.args(opt2$coords, "--coords")
required.args(opt2$geodist, "--geodist")


library("conStruct")
library("parallel")
library("foreach")
library("doParallel")

# If more than one node is specified parallel = TRUE; otherwise FALSE
if (opt2$nodes > 1) {
  parallel <- TRUE
}else{
  parallel <- FALSE
}

# Print parameter settings
print("\nArguments provided: \n")
print(paste0("--maxK = ", opt2$maxK))
print(paste0("--nreps = ", opt2$nreps))
print(paste0("--nodes = ", opt2$nodes))
print(paste0("--trainProp = ", opt2$trainProp))
print(paste0("--saveFiles = ", opt2$saveFiles))
print(paste0("--saveFigs = ", opt2$saveFigs))
print(paste0("--outdir = ", opt2$outdir))
print(paste0("--str = ", opt2$str))
print(paste0("--popmap = ", opt2$popmap))
print(paste0("--geodist = ", opt2$geodist))
print(paste0("--coords = ", opt2$coords))
print(paste0("--wd = ", opt2$wd))
print(paste0("--prefix = ", opt2$prefix))
print(paste0("--niter = ", opt2$niter))
print(paste0("--afreq = ", opt2$afreq))
print(paste0("--onerowperind = ", opt2$onerowperind))
print(paste0("--data_column = ", opt2$data_column))
print(paste0("--missingval = ", opt2$missingval))
print(paste0("--adapt_delta = ", opt2$adapt_delta))
print(paste0("--max_treedepth = ", opt2$max_treedepth))

# Get path to script's directory
script.dir <- getwd()

print(paste("Setting working directory to", getwd(), sep="... "))

# Change to output directory from runConstruct.R
setwd(opt2$wd)

print("\n\nDONE!\n\n")

print("Reading input files...")

# Read input files specified on command-line
str.file <- read.infile(opt2$str, FALSE)
popmap.file <- read.infile(opt2$popmap, FALSE)
coords <- read.infile(opt2$coords, TRUE)
geo.dist <- read.infile(opt2$geodist, FALSE)

print("\n\nDONE!\n\n")

print("Processing and parsing input files...")

# Append K value to afreq filename

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
if (file.exists(paste0(opt2$afreq, ".RData"))) {
  file.remove(paste0(opt2$afreq, ".RData"))
}

print("\n\nDONE!\n\n")

print("Converting STRUCTURE file to conStruct allele frequency format...")

# Convert STRUCTURE file to conStruct allele frequency format
afreq <- structure2conStruct(infile = opt2$str, 
                             onerowperind = opt2$onerowperind, 
                             start.loci = opt2$data_column, 
                             missing.datum = opt2$missingval, 
                             outfile = opt2$afreq)


print("\n\nDONE!\n\n")

# Get population means of allele frequencies if popmap file is provided
if (!is.null(opt2$popmap)){
  print("Getting means of population allele frequencies...")
  pop.data.matrix <- matrix(NA,nrow=nrow(coords),ncol=ncol(afreq))
  for(i in 1:nrow(pop.data.matrix)){
    pop.data.matrix[i,] <- colMeans(
      afreq[
        which(pop.index==i),,
        drop=FALSE
        ],na.rm=TRUE
    )
  }
  print("\n\nDONE!\n\n")
}


# Prefix for cross-validation analyses
xval.prefix <- paste0(opt2$prefix, "_xval")

setwd(script.dir)

# Change to output directory
dir.create(opt2$outdir)
setwd(opt2$outdir)

if (parallel) {
  # For parallelization
  # Ran it like this because had issues with clean exit due to parallel not stopping
  cl <- makeCluster(opt2$nodes, type="FORK")
  registerDoParallel(cl)
}

print("Running Cross-validation analysis...")

# Run cross-validation
conStruct.xvals <- x.validation(train.prop = opt2$trainProp,
                                n.reps = opt2$nreps,
                                K = 1:opt2$maxK,
                                freqs = pop.data.matrix,
                                data.partitions = NULL,
                                geoDist = geo.distMat,
                                coords = coordMat,
                                prefix = xval.prefix,
                                n.iter = opt2$niter,
                                make.figs = opt2$saveFigs,
                                save.files = opt2$saveFiles,
                                parallel = parallel,
                                n.nodes = opt2$nodes,
                                control = setNames(list(opt2$adapt_delta, 
                                                        opt2$max_treedepth), 
                                                   c("adapt_delta", "max_treedepth")))
if(parallel) {
  # End parallelization
  stopCluster(cl)
}

print("\n\nDONE!\n\n")

print("Generating plots...")

# spatial results from cross-validation analysis
sp.results <- as.matrix(
                read.table(paste0(xval.prefix, "_sp_xval_results.txt"),
                                  header = TRUE,
                                  stringsAsFactors = FALSE))

# non-spatial results from cross-validation analysis
nsp.results <- as.matrix(
                read.table(paste0(xval.prefix, "_nsp_xval_results.txt"),
                           header = TRUE,
                           stringsAsFactors = FALSE))

# Generate 95% confidence intervals for spatial and non-spatial
sp.CIs <- apply(sp.results,1,function(x){mean(x) + c(-1.96,1.96) * sd(x)/length(x)})
nsp.CIs <- apply(nsp.results,1,function(x){mean(x) + c(-1.96,1.96) * sd(x)/length(x)})

# Write plots to PDF
pdf(file = paste0(opt2$prefix, "_xval.out.pdf"))

# Make plots for cross-validation results
# Across all K values for both spatial (blue) and non-spatial (green) analyses
par(mfrow=c(1,2))
plot(rowMeans(sp.results),
     pch=19,col="blue",
     ylab="predictive accuracy",xlab="values of K",
     ylim=range(sp.results,nsp.results),
     main="cross-validation results")
points(rowMeans(nsp.results),col="green",pch=19)

# finally, visualize results for the spatial model
# separately with its confidence interval bars
plot(rowMeans(sp.results),
     pch=19,col="blue",
     ylab="predictive accuracy",xlab="values of K",
     ylim=range(sp.CIs),
     main="spatial cross-validation results")
segments(x0 = 1:nrow(sp.results),
         y0 = sp.CIs[1,],
         x1 = 1:nrow(sp.results),
         y1 = sp.CIs[2,],
         col = "blue",lwd=2)

# Now make CI bars for non-spatial
plot(rowMeans(nsp.results),
     pch=19,col="green",
     ylab="predictive accuracy",xlab="values of K",
     ylim=range(nsp.CIs),
     main="non-spatial cross-validation results")
segments(x0 = 1:nrow(nsp.results),
         y0 = nsp.CIs[1,],
         x1 = 1:nrow(nsp.results),
         y1 = nsp.CIs[2,],
         col = "green",lwd=2)

dev.off()

# Make matrix for layer contributions validation
layer.contributions.sp <- matrix(NA,nrow=opt2$maxK,ncol=opt2$maxK)
layer.contributions.nsp <- matrix(NA,nrow=opt2$maxK,ncol=opt2$maxK)

# Get prefixes for Robj files between 1 and --maxK
all.sp <- character()
all.nsp <- character()
for (i in 1:opt2$maxK) {
  all.sp <- c(all.sp, paste0(opt2$prefix, "_xval_sp_rep1K", i))
  all.nsp <- c(all.nsp, paste0(opt2$prefix, "_xval_nsp_rep1K", i))
}

# Load K=1 conStruct.results and data.block Robj files
load(paste0(all.sp[1], "_data.block.Robj"))
load(paste0(all.sp[1], "_conStruct.results.Robj"))

# Calculate layer contributions for K=1
layer.contributions.sp[,1] <- c(calculate.layer.contribution(conStruct.results[[1]],data.block),rep(0,opt2$maxK-1))
tmp.sp <- conStruct.results[[1]]$MAP$admix.proportions

# load all other K values for spatial files
for (i in 2:opt2$maxK) {

  load(paste0(all.sp[i], "_data.block.Robj"))
  load(paste0(all.sp[i], "_conStruct.results.Robj"))
  
  # match layers up across runs to keep plotting colors consistent
  #   for the same layers in different runs
  tmp.order.sp <- match.layers.x.runs(tmp.sp,conStruct.results[[1]]$MAP$admix.proportions)  
  
  # calculate layer contributions
  layer.contributions.sp[,i] <- c(calculate.layer.contribution(conStruct.results=conStruct.results[[1]],
                                                            data.block=data.block,
                                                            layer.order=tmp.order.sp),
                                                             rep(0,opt2$maxK-i))
  
  tmp.sp <- conStruct.results[[1]]$MAP$admix.proportions[,tmp.order.sp]
}

# Load K=1 conStruct.resuls and data.block Robj files
load(paste0(all.nsp[1], "_data.block.Robj"))
load(paste0(all.nsp[1], "_conStruct.results.Robj"))

layer.contributions.nsp[,1] <- c(calculate.layer.contribution(conStruct.results[[1]],data.block),rep(0,opt2$maxK-1))
tmp.nsp <- conStruct.results[[1]]$MAP$admix.proportions

for (i in 2:opt2$maxK) {
# load all other K values for nonspatial files
  load(paste0(all.nsp[i], "_data.block.Robj"))
  load(paste0(all.nsp[i], "_conStruct.results.Robj"))
  
  # match layers up across runs to keep plotting colors consistent
  #   for the same layers in different runs
  tmp.order.nsp <- match.layers.x.runs(tmp.nsp,conStruct.results[[1]]$MAP$admix.proportions)  
  
  # calculate layer contributions
  layer.contributions.nsp[,i] <- c(calculate.layer.contribution(conStruct.results=conStruct.results[[1]],
                                                            data.block=data.block,
                                                            layer.order=tmp.order.nsp),
                               rep(0,opt2$maxK-i))
  
  tmp.nsp <- conStruct.results[[1]]$MAP$admix.proportions[,tmp.order.nsp]
}

# Change back to script directory
setwd(script.dir)

# Change to output directory for this script.
setwd(opt2$outdir)

# Write plots to PDF
pdf(file = paste0(opt2$prefix, "_layerContributions.out.pdf"))

# Make barplot for layer contributions
barplot(layer.contributions.sp,
        col=c("blue", "red", "goldenrod1", "forestgreen", "darkorchid1"),
        xlab="",
        ylab="layer contributions (spatial)",
        names.arg=paste0("K=",1:opt2$maxK))
barplot(layer.contributions.nsp,
        col=c("blue", "red", "goldenrod1", "forestgreen", "darkorchid1"),
        xlab="",
        ylab="layer contributions (non-spatial)",
        names.arg=paste0("K=",1:opt2$maxK))
dev.off()

print("\n\nDONE! Analysis completed!\n")

