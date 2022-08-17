#!/usr/bin/env Rscript

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
  make_option("--idstoIncludeFile", type="character",default="",
    help="Path to a file containing variant ids to be included from the dosage file. The file does not have a header and each line is for a marker ID. Marker ids are in the format chr:pos_ref/alt"),
  make_option("--rangestoIncludeFile", type="character",default="",
    help="Path to a file containing genome regions to be included from the dosage file. The file contains three columns for chromosome, start, and end respectively with no header. Note for vcf and sav files, only the first line in the file will be used."),
  make_option("--chrom", type="character",default="",
    help="If LOCO is specified, chrom is required. chrom is also required for VCF/BCF/SAV input. Note: the string needs to exactly match the chromosome string in the vcf/sav file. For example, 1 does not match chr1."),
  make_option("--is_imputed_data", type="logical",default=FALSE,
    help="Whether the dosages/genotypes imputed are imputed. If TRUE, the program will output the imputed info score [default=FALSE]."),
  make_option("--minMAF", type="numeric", default=0,
    help="Minimum minor allele frequency for markers to be tested. The higher threshold between minMAC and minMAF will be used [default=0]."),
  make_option("--minMAC", type="numeric", default=0.5,
    help="Minimum minor allele count for markers to be tested. The higher threshold between minMAC and minMAF will be used [default=0.5]."),
 make_option("--minGroupMAC_in_BurdenTest", type="numeric", default=5,
    help="Only applied when only Burden tests are performed (r.corr=1). Minimum minor allele count in the Burden test for the psueodo marker. [default=5]."),
  make_option("--minInfo", type="numeric", default=0,
    help="Minimum Info for markers to be tested if is_imputed_data=TRUE [default=0]"),
  make_option("--maxMissing", type="numeric", default=0.15,
    help="Maximum missing rate for markers to be tested [default=0.15]"),
  make_option("--impute_method", type="character",default="best_guess",
    help="Imputation method for missing dosages. best_guess, mean or minor. best_guess: missing dosages imputed as best guessed genotyes round(2*allele frequency). mean: missing dosages are imputed as mean (2*allele frequency). minor: missing dosages are imputed as minor allele homozygotes [default=minor]"),
  make_option("--LOCO", type="logical", default=TRUE,
    help="Whether to apply the leave-one-chromosome-out option. If TRUE, --chrom is required [default=FALSE] "),


  make_option("--GMMATmodelFile", type="character",default="",
    help="Path to the input file containing the glmm model, which is output from previous step. Will be used by load()"),
  make_option("--varianceRatioFile", type="character",default="",
    help="Path to the input file containing the variance ratio, which is output from the previous step"),
  make_option("--SAIGEOutputFile", type="character", default="",
    help="Path to the output file containing assoc test results"),
  make_option("--markers_per_chunk", type="numeric",default=10000,
    help="Number of markers to be tested and output in each chunk in the single-variant assoc tests [default=10000]"),
  make_option("--groups_per_chunk", type="numeric",default=100,
    help="Number of groups/sets to be read in and tested in each chunk in the set-based assoc tests [default=100]"),
  make_option("--is_output_moreDetails", type="logical",default=FALSE,
    help="Whether to output heterozygous and homozygous counts in cases and controls. By default, FALSE. If True, the columns homN_Allele2_cases, hetN_Allelelogical2_cases, homN_Allele2_ctrls, hetN_Allele2_ctrls will be output [default=FALSE]"),
  make_option("--is_overwrite_output", type="logical",default=TRUE,
    help="Whether to overwrite the output file if it exists. If FALSE, the program will continue the unfinished analysis instead of starting over from the beginining [default=TRUE]"),				
  make_option("--maxMAF_in_groupTest", type="character",default="0.0001,0.001,0.01",
    help="Max MAF for markers tested in group test seperated by comma. e.g. 0.0001,0.001,0.01, [default=0.01]"),
  make_option("--maxMAC_in_groupTest", type="character",default="0",
    help="Max MAC for markers tested in group test seperated by comma.The list will be combined with maxMAF_in_groupTest. e.g. 1,2 . By default, 0 and no maxMAC cutoff are applied. [default=0]"),
  make_option("--annotation_in_groupTest", type="character",default="lof,missense;lof,missense;lof;synonymous",
    help="annotations of markers to be tested in the set-based tests seperated by comma. using ; to combine multiple annotations in the same test, e.g. lof,missense;lof,missense;lof;synonymous will test lof variants only, missense+lof variants, and missense+lof+synonymous variants. default: lof,missense;lof,missense;lof;synonymous"),
  make_option("--groupFile", type="character", default="",
    help="Path to the file containing the group information for gene-based tests. Each gene/set has 2 or 3 lines in the group file. The first element is the gene/set name. The second element in the first line is to indicate whether this line contains variant IDs (var), annotations (anno), or weights (weight). The line for weights is optional. If not specified, the default weights will be generated based on beta(MAF, 1, 25). Use --weights.beta to change the parameters for the Beta distribution. The variant ids must be in the format chr:pos_ref/alt. Elements are seperated by tab or space."),
  make_option("--sparseGRMFile", type="character", default="",
   help="Path to the pre-calculated sparse GRM file that was used in Step 1"),
  make_option("--sparseGRMSampleIDFile", type="character", default="",
   help="Path to the sample ID file for the pre-calculated sparse GRM. No header is included. The order of sample IDs is corresponding to sample IDs in the sparse GRM"),
  make_option("--relatednessCutoff", type="numeric", default=0,
   help="Optional. Threshold (minimum realtedness coefficient) to treat two samples as unrelated when the sparse GRM is used [default=0]"), 
  make_option("--MACCutoff_to_CollapseUltraRare", type="numeric", default=10,
    help="MAC cutoff to collpase the ultra rare variants (<= MACCutoff_to_CollapseUltraRare) in the set-based association tests. By default, 10."),
  make_option("--cateVarRatioMinMACVecExclude",type="character", default="10,20.5",
    help="Optional. vector of float. Lower bound for MAC categories. The length equals to the number of MAC categories for variance ratio estimation. [default='10,20.5']"),
  make_option("--cateVarRatioMaxMACVecInclude",type="character", default="20.5",
    help="Optional. vector of float. Higher bound for MAC categories. The length equals to the number of MAC categories for variance ratio estimation minus 1. [default='20.5']"),
  make_option("--weights.beta", type="character", default="1,25",
    help="parameters for the beta distribution to weight genetic markers in gene-based tests."),
  make_option("--r.corr", type="numeric", default=0,
    help="If r.corr = 1, only Burden tests will be performed. If r.corr = 0, SKAT-O tests will be performed and results for Burden tests and SKAT tests will be output too. [default = 0]"),
  make_option("--markers_per_chunk_in_groupTest", type="numeric", default=100,
    help="Number of markers in each chunk when calculating the variance covariance matrix in the set/group-based tests  [default=100]."),
  make_option("--condition", type="character",default="",
    help="For conditional analysis. Variant ids are in the format chr:pos_ref/alt and seperated by by comma. e.g.chr3:101651171_C/T,chr3:101651186_G/A"),
  make_option("--weights_for_condition",type="character", default=NULL,
    help="vector of float. weights for conditioning markers for gene- or region-based tests. The length equals to the number of conditioning markers, delimited by comma. e.g. '1,2,3. If not specified, the default weights will be generated based on beta(MAF, 1, 25). Use --weights.beta to change the parameters for the Beta distribution."),
  make_option("--SPAcutoff", type="numeric", default=2,
    help=" If the test statistic lies within the standard deviation cutoff of the
mean, p-value based on traditional score test is returned. Default value is 2."),
  make_option("--dosage_zerod_cutoff",type="numeric", default=0.2,
    help="If is_imputed_data = TRUE, For variants with MAC <= dosage_zerod_MAC_cutoff, dosages <= dosageZerodCutoff with be set to 0. [default=0.2]"),
  make_option("--dosage_zerod_MAC_cutoff", type="numeric", default=10,
   help="If is_imputed_data = TRUE, For variants with MAC <= dosage_zerod_MAC_cutoff, dosages <= dosageZerodCutoff with be set to 0. [default=10]"),

  make_option("--is_single_in_groupTest", type="logical", default=FALSE,
    help="Whether to output single-variant assoc test results when perform group tests. Note, single-variant assoc test results will always be output when SKAT and SKAT-O tests are conducted with --r.corr=0. This parameter should only be used when only Burden tests are condcuted with --r.corr=1, [default=TRUE]"), 
  make_option("--is_no_weight_in_groupTest", type="logical", default=FALSE,
    help="Whether no weights are used in group Test. If FALSE, weights will be calcuated based on MAF from the Beta distribution with paraemters weights.beta or weights will be extracted from the group File if available [default=FALSE]"),
  make_option("--is_output_markerList_in_groupTest", type="logical", default=FALSE,
    help="Whether to output the marker lists included in the set-based tests for each mask.[default=TRUE]"),
  make_option("--is_Firth_beta", type="logical", default=FALSE,
    help="Whether to estimate effect sizes using approx Firth, only for binary traits [default=FALSE]"),
  make_option("--pCutoffforFirth", type="numeric", default=0.01,
    help="p-value cutoff to use approx Firth to estiamte the effect sizes. Only for binary traits. The effect sizes of markers with p-value <= pCutoffforFirth will be estimated using approx Firth [default=0.01]"),
  make_option("--is_fastTest", type="logical", default=FALSE,
    help="Whether to use the fast mode for tests"),
  make_option("--max_MAC_for_ER", type="numeric", default=4,
    help="p-values of genetic variants with MAC <= max_MAC_for_ER will be calculated via efficient resampling. [default=4]"),

 make_option("--subSampleFile", type="character",default="",
    help="Path to the file that contains one column for IDs of samples that are included in Step 1 and will be also included in Step 2. This option is used when any sample included in Step 1 but does not have dosages/genotypes for Step 2. Please make sure it contains one column of sample IDs that will be used for subsetting samples from the Step 1 results for Step 2 jobs. Note: Thi option has not been fully evaluated. If more than 5% of the samples in Step 1 are missing in Step 2, please consider re-run Step 1. ")
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
	#	cat(stringOutput, " is ", y, "\n")
	#}
	return(y)	
}


#weights.beta.rare <- as.numeric(strsplit(opt$weights.beta.rare,",")[[1]])
#weights.beta.rare <- convertoNumeric(x=strsplit(opt$weights.beta.rare,",")[[1]], "weights.beta.rare")
#weights.beta.common <- convertoNumeric(x=strsplit(opt$weights.beta.common,",")[[1]], "weights.beta.common")
#if(sum(weights.beta.common!=weights.beta.rare) > 0){stop("weights.beta.common option is not functioning, so weights.beta.common needs to be equal to weights.beta.rare")}
weights.beta <- convertoNumeric(x=strsplit(opt$weights.beta,",")[[1]], "weights.beta")

cateVarRatioMinMACVecExclude <- convertoNumeric(x=strsplit(opt$cateVarRatioMinMACVecExclude,",")[[1]], "cateVarRatioMinMACVecExclude")
cateVarRatioMaxMACVecInclude <- convertoNumeric(x=strsplit(opt$cateVarRatioMaxMACVecInclude,",")[[1]], "cateVarRatioMaxMACVecInclude")
if(is.null(opt$weights_for_condition)){
	weights_for_condition=NULL
}else{
	weights_for_condition <- convertoNumeric(x=strsplit(opt$weights_for_condition,",")[[1]], "weights_for_condition")
}

maxMAF_in_groupTest = convertoNumeric(x=strsplit(opt$maxMAF_in_groupTest,",")[[1]], "maxMAF_in_groupTest")
maxMAC_in_groupTest = convertoNumeric(x=strsplit(opt$maxMAC_in_groupTest,",")[[1]], "maxMAC_in_groupTest")

annotation_in_groupTest = gsub(":",";",opt$annotation_in_groupTest)
annotation_in_groupTest = unlist(strsplit(annotation_in_groupTest,",")[[1]])

#try(if(length(which(opt == "")) > 0) stop("Missing arguments"))

##by Alex Petty @pettyalex
if (BLASctl_installed){
  # Set number of threads for BLAS to 1, this step does not benefit from multithreading or multiprocessing
  original_num_threads <- blas_get_num_procs()
  blas_set_num_threads(1)
}

print("opt$r.corr")
print(opt$r.corr)

if(packageVersion("SAIGE")<"1.1.3"){


  SPAGMMATtest(vcfFile=opt$vcfFile,
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
	     idstoIncludeFile = opt$idstoIncludeFile,
	     rangestoIncludeFile = opt$rangestoIncludeFile,
	     chrom=opt$chrom,
             is_imputed_data=opt$is_imputed_data,
             min_MAF = opt$minMAF,
             min_MAC = opt$minMAC,
             min_Info = opt$minInfo,
             max_missing = opt$maxMissing,	
	     impute_method = opt$impute_method,
	     LOCO=opt$LOCO,
             GMMATmodelFile=opt$GMMATmodelFile,
             varianceRatioFile=opt$varianceRatioFile,
             SAIGEOutputFile=opt$SAIGEOutputFile,	
	     markers_per_chunk=opt$markers_per_chunk,
	     groups_per_chunk=opt$groups_per_chunk,
             markers_per_chunk_in_groupTest=opt$markers_per_chunk_in_groupTest,
	     is_output_moreDetails =opt$is_output_moreDetails,
	     is_overwrite_output = opt$is_overwrite_output,
	     maxMAF_in_groupTest = maxMAF_in_groupTest,
	     maxMAC_in_groupTest = maxMAC_in_groupTest,
	     minGroupMAC_in_BurdenTest = opt$minGroupMAC_in_BurdenTest,
	     annotation_in_groupTest = annotation_in_groupTest,
	     groupFile = opt$groupFile,
	     sparseGRMFile=opt$sparseGRMFile,
             sparseGRMSampleIDFile=opt$sparseGRMSampleIDFile,
	     relatednessCutoff=opt$relatednessCutoff,	
	     MACCutoff_to_CollapseUltraRare = opt$MACCutoff_to_CollapseUltraRare,	
	     cateVarRatioMinMACVecExclude = cateVarRatioMinMACVecExclude,
             cateVarRatioMaxMACVecInclude = cateVarRatioMaxMACVecInclude,
	     weights.beta = weights.beta,
	     r.corr = opt$r.corr,
	     condition = opt$condition,
	     weights_for_condition = weights_for_condition, 
	     SPAcutoff = opt$SPAcutoff,
	     dosage_zerod_cutoff = opt$dosage_zerod_cutoff,
	     dosage_zerod_MAC_cutoff = opt$dosage_zerod_MAC_cutoff,
	     is_Firth_beta = opt$is_Firth_beta,
	     pCutoffforFirth = opt$pCutoffforFirth,
	     is_single_in_groupTest = opt$is_single_in_groupTest,
             is_no_weight_in_groupTest = opt$is_no_weight_in_groupTest,
	     is_output_markerList_in_groupTest = opt$is_output_markerList_in_groupTest,
	     is_fastTest = opt$is_fastTest
)
}else{

if(packageVersion("SAIGE")>"1.1.4"){	
  SPAGMMATtest(vcfFile=opt$vcfFile,
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
             idstoIncludeFile = opt$idstoIncludeFile,
             rangestoIncludeFile = opt$rangestoIncludeFile,
             chrom=opt$chrom,
             is_imputed_data=opt$is_imputed_data,
             min_MAF = opt$minMAF,
             min_MAC = opt$minMAC,
             min_Info = opt$minInfo,
             max_missing = opt$maxMissing,
             impute_method = opt$impute_method,
             LOCO=opt$LOCO,
             GMMATmodelFile=opt$GMMATmodelFile,
             varianceRatioFile=opt$varianceRatioFile,
             SAIGEOutputFile=opt$SAIGEOutputFile,
             markers_per_chunk=opt$markers_per_chunk,
             groups_per_chunk=opt$groups_per_chunk,
             markers_per_chunk_in_groupTest=opt$markers_per_chunk_in_groupTest,
             is_output_moreDetails =opt$is_output_moreDetails,
             is_overwrite_output = opt$is_overwrite_output,
             maxMAF_in_groupTest = maxMAF_in_groupTest,
             maxMAC_in_groupTest = maxMAC_in_groupTest,
             minGroupMAC_in_BurdenTest = opt$minGroupMAC_in_BurdenTest,
             annotation_in_groupTest = annotation_in_groupTest,
             groupFile = opt$groupFile,
             sparseGRMFile=opt$sparseGRMFile,
             sparseGRMSampleIDFile=opt$sparseGRMSampleIDFile,
             relatednessCutoff=opt$relatednessCutoff,
             MACCutoff_to_CollapseUltraRare = opt$MACCutoff_to_CollapseUltraRare,
             cateVarRatioMinMACVecExclude = cateVarRatioMinMACVecExclude,
             cateVarRatioMaxMACVecInclude = cateVarRatioMaxMACVecInclude,
             weights.beta = weights.beta,
             r.corr = opt$r.corr,
             condition = opt$condition,
             weights_for_condition = weights_for_condition,
             SPAcutoff = opt$SPAcutoff,
             dosage_zerod_cutoff = opt$dosage_zerod_cutoff,
             dosage_zerod_MAC_cutoff = opt$dosage_zerod_MAC_cutoff,
             is_Firth_beta = opt$is_Firth_beta,
             pCutoffforFirth = opt$pCutoffforFirth,
             is_single_in_groupTest = opt$is_single_in_groupTest,
             is_no_weight_in_groupTest = opt$is_no_weight_in_groupTest,
             is_output_markerList_in_groupTest = opt$is_output_markerList_in_groupTest,
             is_fastTest = opt$is_fastTest,
	     max_MAC_use_ER = opt$max_MAC_for_ER,
	     subSampleFile = opt$subSampleFile
)
  }else{

	SPAGMMATtest(vcfFile=opt$vcfFile,
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
             idstoIncludeFile = opt$idstoIncludeFile,
             rangestoIncludeFile = opt$rangestoIncludeFile,
             chrom=opt$chrom,
             is_imputed_data=opt$is_imputed_data,
             min_MAF = opt$minMAF,
             min_MAC = opt$minMAC,
             min_Info = opt$minInfo,
             max_missing = opt$maxMissing,
             impute_method = opt$impute_method,
             LOCO=opt$LOCO,
             GMMATmodelFile=opt$GMMATmodelFile,
             varianceRatioFile=opt$varianceRatioFile,
             SAIGEOutputFile=opt$SAIGEOutputFile,
             markers_per_chunk=opt$markers_per_chunk,
             groups_per_chunk=opt$groups_per_chunk,
             markers_per_chunk_in_groupTest=opt$markers_per_chunk_in_groupTest,
             is_output_moreDetails =opt$is_output_moreDetails,
             is_overwrite_output = opt$is_overwrite_output,
             maxMAF_in_groupTest = maxMAF_in_groupTest,
             maxMAC_in_groupTest = maxMAC_in_groupTest,
             minGroupMAC_in_BurdenTest = opt$minGroupMAC_in_BurdenTest,
             annotation_in_groupTest = annotation_in_groupTest,
             groupFile = opt$groupFile,
             sparseGRMFile=opt$sparseGRMFile,
             sparseGRMSampleIDFile=opt$sparseGRMSampleIDFile,
             relatednessCutoff=opt$relatednessCutoff,
             MACCutoff_to_CollapseUltraRare = opt$MACCutoff_to_CollapseUltraRare,
             cateVarRatioMinMACVecExclude = cateVarRatioMinMACVecExclude,
             cateVarRatioMaxMACVecInclude = cateVarRatioMaxMACVecInclude,
             weights.beta = weights.beta,
             r.corr = opt$r.corr,
             condition = opt$condition,
             weights_for_condition = weights_for_condition,
             SPAcutoff = opt$SPAcutoff,
             dosage_zerod_cutoff = opt$dosage_zerod_cutoff,
             dosage_zerod_MAC_cutoff = opt$dosage_zerod_MAC_cutoff,
             is_Firth_beta = opt$is_Firth_beta,
             pCutoffforFirth = opt$pCutoffforFirth,
             is_single_in_groupTest = opt$is_single_in_groupTest,
             is_no_weight_in_groupTest = opt$is_no_weight_in_groupTest,
             is_output_markerList_in_groupTest = opt$is_output_markerList_in_groupTest,
             is_fastTest = opt$is_fastTest,
             max_MAC_use_ER = opt$max_MAC_for_ER
)



}	

}	
if(BLASctl_installed){
  # Restore originally configured BLAS thread count
  blas_set_num_threads(original_num_threads)
}
