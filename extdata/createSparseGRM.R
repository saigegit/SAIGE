#!/usr/bin/env -S pixi run --manifest-path /app/pixi.toml Rscript

options(stringsAsFactors=F)

## load R libraries
#library(SAIGE, lib.loc="../../install_0.93")
library(SAIGE)
require(optparse) #install.packages("optparse")

print(sessionInfo())

## set list of cmd line arguments
option_list <- list(
  make_option("--plinkFile", type="character",default="",
    help="path to plink file for creating the genetic relationship matrix (GRM)"),
  make_option("--bedFile", type="character",default="",
    help="Path to bed file. If plinkFile is specified, 'plinkFile'.bed will be used"),
  make_option("--bimFile", type="character",default="",
    help="Path to bim file. If plinkFile is specified, 'plinkFile'.bim will be used"),
  make_option("--famFile", type="character",default="",
    help="Path to fam file. If plinkFile is specified, 'plinkFile'.fam will be used"),
  make_option("--nThreads", type="integer", default=16,
    help="Number of threads (CPUs) to use"),
  make_option("--memoryChunk", type="numeric", default=2,
   help="The size (Gb) for each memory chunk. By default, 2"),
  make_option("--outputPrefix", type="character", default="~/",
    help="path and prefix to the output files [default='~/']"),
  make_option("--numRandomMarkerforSparseKin", type="integer", default=2000,
    help="number of randomly selected markers to be used to identify related samples for sparse GRM [default=2000]"),
  make_option("--relatednessCutoff", type="numeric", default=0.125,
    help="The threshold to treat two samples as unrelated if IsSparseKin is TRUE [default=0.125]"),
  make_option("--isDiagofKinSetAsOne", type="logical", default=FALSE,
    help="Whether to set the diagnal elements in GRM to be 1 [default='FALSE']."),
  make_option("--minMAFforGRM", type="numeric", default=0.01,
    help="minimum MAF of markers used for GRM"),
  make_option("--maxMissingRateforGRM", type="numeric", default=0.15,
    help="Optional. Maximum missing rate of markers used for GRM")
)



## list of options
parser <- OptionParser(usage="%prog [options]", option_list=option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

set.seed(1)

if(opt$plinkFile != ""){
        bimFile = paste0(opt$plinkFile, ".bim")
        bedFile = paste0(opt$plinkFile, ".bed")
        famFile = paste0(opt$plinkFile, ".fam")	
}else{
	bimFile = opt$bimFile
	bedFile = opt$bedFile
	famFile = opt$famFile
}

createSparseGRM(bedFile = bedFile,
		bimFile = bimFile,
		famFile = famFile,
                outputPrefix=opt$outputPrefix,
                numRandomMarkerforSparseKin = opt$numRandomMarkerforSparseKin,
                relatednessCutoff = opt$relatednessCutoff,
                memoryChunk = opt$memoryChunk,
                isDiagofKinSetAsOne = opt$isDiagofKinSetAsOne,
                nThreads = opt$nThreads,
		minMAFforGRM = opt$minMAFforGRM,
		maxMissingRateforGRM = opt$maxMissingRateforGRM)
