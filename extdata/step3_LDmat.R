#!/usr/bin/env -S pixi run --manifest-path /app/pixi.toml Rscript

#options(stringsAsFactors=F, scipen = 999)
options(stringsAsFactors=F)
library(SAIGE)
BLASctl_installed <- require(RhpcBLASctl)
library(optparse)
library(data.table)
library(methods)
print(sessionInfo())


option_list <- list(
  make_option("--vcfFile", type="character",default="",
    help="Path to vcf file."),
  make_option("--vcfFileIndex", type="character",default="",
    help="Path to vcf index file. Indexed by tabix. Path to index for vcf file by tabix, .csi file by tabix -p vcf csi file.vcf.gz"),
  make_option("--vcfField", type="character",default="DS",
    help="DS or GT, [default=DS]"),
  make_option("--savFile", type="character",default="",
    help="Path to the sav file."),
  make_option("--savFileIndex", type="character",default="",
    help="Path to the .s1r file (index of the sav file)."),
  make_option("--bgenFile", type="character",default="",
    help="Path to bgen file. Path to bgen file. Currently version 1.2 with 8 bit compression is supported"),
  make_option("--bgenFileIndex", type="character",default="",
    help="Path to the .bgi file (index of the bgen file)"),
  make_option("--sampleFile", type="character",default="",
    help="Path to the file that contains one column for IDs of samples in the dosage file. For version >= 0.38, this file is only needed for bgen files. "),
  make_option("--bedFile", type="character",default="",
    help="Path to bed file (PLINK)"),
  make_option("--bimFile", type="character",default="",
    help="Path to bim file (PLINK)"),
  make_option("--famFile", type="character",default="",
    help="Path to fam file (PLINK)"),
  make_option("--AlleleOrder", type="character",default="alt-first",
    help="alt-first or ref-first for bgen or PLINK files"),
  make_option("--chrom", type="character",default="",
    help="If LOCO is specified, chrom is required. chrom is also required for VCF/BCF/SAV input. Note: the string needs to exactly match the chromosome string in the vcf/sav file. For example, 1 does not match chr1."),

  make_option("--is_imputed_data", type="logical",default=FALSE,
    help="Whether the dosages/genotypes imputed are imputed. If TRUE, the program will output the imputed info score [default=FALSE]."),
  make_option("--minMAF", type="numeric", default=0,
    help="Minimum minor allele frequency for markers to be tested. The higher threshold between minMAC and minMAF will be used [default=0]."),
  make_option("--minMAC", type="numeric", default=0.5,
    help="Minimum minor allele count for markers to be tested. The higher threshold between minMAC and minMAF will be used [default=0.5]."),
  make_option("--minInfo", type="numeric", default=0,
    help="Minimum Info for markers to be tested if is_imputed_data=TRUE [default=0]"),
  make_option("--maxMissing", type="numeric", default=0.15,
    help="Maximum missing rate for markers to be tested [default=0.15]"),
  make_option("--impute_method", type="character",default="best_guess",
    help="Imputation method for missing dosages. best_guess, mean or minor. best_guess: missing dosages imputed as best guessed genotyes round(2*allele frequency). mean: missing dosages are imputed as mean (2*allele frequency). minor: missing dosages are imputed as minor allele homozygotes [default=minor]"),
  make_option("--SAIGEOutputFile", type="character", default="",
    help="Path to the output file containing assoc test results"),
  make_option("--groups_per_chunk", type="numeric",default=100,
    help="Number of groups/sets to be read in and tested in each chunk in the set-based assoc tests [default=100]"),
  make_option("--maxMAF_in_groupTest", type="numeric",default=0.5,
    help="Max MAF for markers included in the LD matrices, [default=0.5]"),
  make_option("--annotation_in_groupTest", type="character",default="missense;lof;synonymous",
    help="annotations of markers included in the LD matrices using ; to combine multiple annotations in the same test, e.g. missense;lof;synonymous will include missense+lof+synonymous variants. use ALL to include all markers in the genotype/dosage files.  default: missense;lof;synonymous"),
  make_option("--groupFile", type="character", default="",
    help="Path to the file containing the group information for gene-based tests. Each gene/set has 2 or 3 lines in the group file. The first element is the gene/set name. The second element in the first line is to indicate whether this line contains variant IDs (var), annotations (anno), or weights (weight). The line for weights is optional. If not specified, the default weights will be generated based on beta(MAF, 1, 25). Use --weights.beta to change the parameters for the Beta distribution. The variant ids must be in the format chr:pos_ref/alt. Elements are seperated by tab or space."),
  make_option("--markers_per_chunk_in_groupTest", type="numeric", default=100,
    help="Number of markers in each chunk when calculating the variance covariance matrix in the set/group-based tests  [default=100]."),
  make_option("--dosage_zerod_cutoff",type="numeric", default=0.2,
    help="If is_imputed_data = TRUE, For variants with MAC <= dosage_zerod_MAC_cutoff, dosages <= dosageZerodCutoff with be set to 0. [default=0.2]"),
  make_option("--dosage_zerod_MAC_cutoff", type="numeric", default=10,
   help="If is_imputed_data = TRUE, For variants with MAC <= dosage_zerod_MAC_cutoff, dosages <= dosageZerodCutoff with be set to 0. [default=10]"),
   make_option("--sample_include_inLDMat_File", type="character",default="",
    help="Path to the file that contains one column for IDs of samples who are included for LD estimation. If empty, all samples in the dosage/genotype file will be used. By default, empty"),
make_option("--is_overwrite_output", type="logical",default=TRUE,
    help="Whether to overwrite the output file if it exists. If FALSE, the program will continue the unfinished analysis instead of starting over from the beginining [default=TRUE]")  
)



parser <- OptionParser(usage="%prog [options]", option_list=option_list)

args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)


convertoNumeric = function(x,stringOutput){
        y= tryCatch(expr = as.numeric(x),warning = function(w) {return(NULL)})
        if(is.null(y)){
                stop(stringOutput, " is not numeric\n")
        }
        #else{
        #       cat(stringOutput, " is ", y, "\n")
        #}
        return(y)
}




##by Alex Petty @pettyalex
if (BLASctl_installed){
  # Set number of threads for BLAS to 1, this step does not benefit from multithreading or multiprocessing
  original_num_threads <- blas_get_num_procs()
  blas_set_num_threads(1)
}

if(opt$annotation_in_groupTest == "ALL"){
	opt$annotation_in_groupTest = NULL
}	


if(packageVersion("SAIGE")>="1.1.3"){
	generate_LDMat_forMataRegion(vcfFile=opt$vcfFile,
             vcfFileIndex=opt$vcfFileIndex,
             vcfField=opt$vcfField,
             savFile=opt$savFile,
             savFileIndex=opt$savFileIndex,
             bgenFile=opt$bgenFile,
             bgenFileIndex=opt$bgenFileIndex,
             sampleFile=opt$sampleFile,
             bedFile=opt$bedFile,
             bimFile=opt$bimFile,
             famFile=opt$famFile,
             AlleleOrder=opt$AlleleOrder,
	     chrom=opt$chrom,
	     is_imputed_data=opt$is_imputed_data,
             min_MAF = opt$minMAF,
             min_MAC = opt$minMAC,
             min_Info = opt$minInfo,
             max_missing = opt$maxMissing,
             impute_method = opt$impute_method,		     
	     SAIGEOutputFile=opt$SAIGEOutputFile,  			     
	     groups_per_chunk=opt$groups_per_chunk,
	     markers_per_chunk_in_groupTest=opt$markers_per_chunk_in_groupTest,	
             groupFile = opt$groupFile,
     	     dosage_zerod_cutoff = opt$dosage_zerod_cutoff,
             dosage_zerod_MAC_cutoff = opt$dosage_zerod_MAC_cutoff,
	     annotation_in_groupTest = opt$annotation_in_groupTest,
	     maxMAF_in_groupTest=opt$maxMAF_in_groupTest,
             is_overwrite_output = opt$is_overwrite_output,
             sampleFile_include_inLDMat = opt$sample_include_inLDMat_File)

}	
