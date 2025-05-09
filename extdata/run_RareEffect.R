# Read source and libraries
library(SAIGE)
library(data.table)
library(Matrix)
library(optparse)

option_list <- list(
    make_option(c("--rdaFile"), type="character", default="",
        help="Path to rda file from SAIGE step 1"),
    make_option(c("--chrom"), type="character", default="",
        help="Chromosome"),
    make_option(c("--geneName"), type="character", default="",
        help="Gene name to analyze"),
    make_option(c("--groupFile"), type="character", default="",
        help="Path to group file (containing functional annotation of variants)"),
    make_option(c("--traitType"), type="character", default="",
        help="Trait type (quantitative or binary)"),
    make_option(c("--bedFile"), type="character", default="",
        help="Path to bed file containing rare variants"),
    make_option(c("--bimFile"), type="character", default="",
        help="Path to bim file containing rare variants"),
    make_option(c("--famFile"), type="character", default="",
        help="Path to fam file containing rare variants"),
    make_option(c("--macThreshold"), type="integer", default=10,
        help="MAC threshold for ultra-rare variant collapsing (default: 10)"),
    make_option(c("--collapseLoF"), type="logical", default=FALSE,
        help="Collapse LoF variants into one super-variant (default: FALSE)"),
    make_option(c("--collapsemis"), type="logical", default=FALSE,
        help="Collapse missense variants into one super-variant (default: FALSE)"),
    make_option(c("--collapsesyn"), type="logical", default=FALSE,
        help="Collapse synonymous variants into one super-variant (default: FALSE)"),
    make_option(c("--apply_AR"), type="logical", default=FALSE,
        help="Apply adaptive ridge (approximated L0-regularization) when estimating effect size (defualt: FALSE)"),
    make_option(c("--outputPrefix"), type="character", default="",
        help="Path to save output (without extension and gene name)")
)

## list of options
parser <- OptionParser(usage="%prog [options]", option_list=option_list)
args <- parse_args(parser, positional_arguments = FALSE)
print(args)

run_RareEffect(
    rdaFile = args$rdaFile,
    chrom = args$chrom,
    geneName = args$geneName,
    groupFile = args$groupFile,
    traitType = args$traitType,
    bedFile = args$bedFile,
    bimFile = args$bimFile,
    famFile = args$famFile,
    macThreshold = args$macThreshold,
    collapseLoF = args$collapseLoF,
    collapsemis = args$collapsemis,
    collapsesyn = args$collapsesyn,
    apply_AR = args$apply_AR,
    outputPrefix = args$outputPrefix
)
