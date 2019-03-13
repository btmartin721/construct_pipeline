#!/usr/bin/env Rscript

if (!require("optparse")) stop("Error: The required package optparse is not installed")
library("optparse")

# Set command-line arguments
option_list <- list(make_option(c("--minK"),
                                type="integer",
                                default=NULL,
                                help="Specify minimum K value",
                                metavar="integer"),
                    make_option(c("--maxK"),
                                type="integer",
                                default=NULL,
                                help="Specify maximum K value",
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
                                metavar="integer"))

opt_parser <- OptionParser(option_list=option_list,
                           description="Rscript to run conStruct")

opt2 <- parse_args(opt_parser)

######################################################################

library("conStruct")

# load runConstruct.R environment
load("environment.RData")

# load spatial files
load(paste0(sp.prefix, "_data.block.Robj"))
load(paste0(sp.prefix, "_conStruct.results.Robj"))

# load nonspatial files
load(paste0(nsp.prefix, "_data.block.Robj"))
load(paste0(nsp.prefix, "_conStruct.results.Robj"))