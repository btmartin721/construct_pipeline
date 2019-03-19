#!/usr/bin/env Rscript

if (!require("optparse")) stop("Error: The required package optparse is not installed")
library("optparse")

# Set command-line arguments
option_list <- list(make_option(c("-K", "--maxK"),
                                type="integer",
                                default=NULL,
                                help="Specify maximum K value (upper case K)",
                                metavar="integer"),
                    make_option(c("-r", "--nreps"),
                                type="integer",
                                default=NULL,
                                help="Specify number of cross-validation replicates",
                                metavar="integer"),
                    make_option(c("-t", "--trainProp"),
                                type="numeric",
                                default=0.9,
                                help="Specify training proportion for cross-validation",
                                metavar="numeric"),
                    make_option(c("-n", "--nodes"),
                                type="integer",
                                default=NULL,
                                help="Specify number of CPU cores for parallelization",
                                metavar="integer"),
                    make_option(c("-f", "--saveFiles"),
                                action="store_true",
                                default = FALSE,
                                help="Boolean; Toggle save Robj files from cross-validation; will save lots of files; default=FALSE"),
                    make_option(c("-F", "--saveFigs"),
                                action="store_false",
                                default=TRUE,
                                help="Boolean; Don't save figures from cross-validation; default=TRUE"),
                    make_option(c("-o", "--outdir"),
                                type="character",
                                default="./",
                                help="Specify directory for output files; will be created if doesn't exist",
                                metavar="character"),
                    make_option(c("-w", "--wd"),
                                type="character",
                                default=NULL,
                                help="Specify directory containing output from runConstruct.R",
                                metavar="character"))
                    

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


#####################################################################3
# Require these arguments or stop
required.args(opt2$maxK, "--maxK")
required.args(opt2$nreps, "--nreps")
required.args(opt2$nodes, "--nodes")
required.args(opt2$wd, "--wd")

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

# Get path to script's directory
script.dir <- getwd()

# Change to output directory from runConstruct.R
setwd(opt2$wd)

if (!file.exists("arguments.RDS")) {
  stop(paste0("\n\nError: The file arguments.RDS does not exist in ", opt2$wd, "; this directory should contain all results from runConstruct.R\n"))
}

# Read command-line arguments from runConstruct.R
opt <- readRDS(file = "arguments.RDS")

# Make sure runConstruct.R has been run first to generate the environment.RData file
if (!file.exists(paste0(opt$prefix, "_environment.RData"))) {
  stop(paste0("\n\nError: The file ", opt$prefix, "_environment.RData could not be found in ", opt2$wd, "; this directory should contain all results from runConstruct.R\n"))
}

# load runConstruct.R environment
load(paste0(opt$prefix, "_environment.RData"))

# Prefix for cross-validation analyses
xval.prefix <- paste0(opt$prefix, "_xval")

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

# Run cross-validation
conStruct.xvals <- x.validation(train.prop = opt2$trainProp,
                                n.reps = opt2$nreps,
                                K = 1:opt2$maxK,
                                freqs = pop.data.matrix,
                                data.partitions = NULL,
                                geoDist = geo.distMat,
                                coords = coordMat,
                                prefix = xval.prefix,
                                n.iter = opt$niter,
                                make.figs = opt2$saveFigs,
                                save.files = opt2$saveFiles,
                                parallel = parallel,
                                n.nodes = opt2$nodes)
if(parallel) {
  # End parallelization
  stopCluster(cl)
}

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
pdf(file = paste0(opt$prefix, "_xval.out.pdf"))

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
  all.sp <- c(all.sp, paste0(opt$prefix, "_spK", i))
  all.nsp <- c(all.nsp, paste0(opt$prefix, "_nspK", i))
}

# Change to submission directory
setwd(script.dir)

# Change to runConstruct.R directory
setwd(opt2$wd)

# Load K=1 conStruct.resuls and data.block Robj files
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
pdf(file = paste0(opt$prefix, "_layerContributions.out.pdf"))

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

