getChromNumber = function(chrom = ""){
  if(chrom == ""){
    stop("chrom is not specified\n")
  }else{	  
    chrom_v2 = as.character(chrom)
    chrom_v2 = gsub("CHR", "", chrom_v2, ignore.case = T)
    chrom_v3 = as.numeric(gsub("[^0-9.]", "",chrom_v2))

   if(!is.na(chrom_v3)){

   
    if(chrom_v3 <= 22 &  chrom_v3 >= 1) {
      #stop("chromosome ", chrom, " is out of the range of null model LOCO results\n")
    #}else {
      cat("Leave chromosome ", chrom_v3, " out will be applied\n")
    }else{
      chrom_v3 = chrom	
    }
   }else{
	chrom_v3 = chrom
   }
  }
  
  return(chrom_v3)
}	

removeLOCOResult = function(chromList, obj.glmm.null){
  for (chr in 1:22) {
    if(chr %in% chromList){    
      obj.glmm.null$LOCOResult[chr] = list(NULL)
      cat("chromosome ", chr, " model results are removed to save memory\n")
      gc()
    }
  }
  return(obj.glmm.null)
}	


ReadModel <- function(GMMATmodelFile = "", chrom = "", LOCO = TRUE, is_Firth_beta = FALSE, is_EmpSPA = FALSE, espa_nt = 9999, espa_range = c(-20, 20)) {
  # Check file existence
  Check_File_Exist(GMMATmodelFile, "GMMATmodelFile")
  if(!LOCO %in% c(TRUE, FALSE))
    stop("LOCO should be TRUE or FALSE.")
  # load GMMATmodelFile
  load(GMMATmodelFile)
  modglmm$Y = NULL
  modglmm$linear.predictors = NULL
  modglmm$coefficients = NULL
  modglmm$cov = NULL
  gc()
  chrom_v3=NULL

  if(!LOCO) {
    print("Leave-one-chromosome-out is not applied")
    if(modglmm$LOCO) {
      for (chr in 1:22) {
        modglmm$LOCOResult[chr] = list(NULL)
        cat("chromosome ", chr, " model results are removed to save memory\n")
        gc()
      }
    }
  }else{
    if (!modglmm$LOCO){
      stop("LOCO is TRUE but the null model file .rda does not contain LOCO results. In order to apply Leave-one-chromosome-out, please run Step 1 using LOCO. Otherwise, please set LOCO=FALSE in this step (Step 2).\n")
    }else{
        if(chrom == ""){
          stop("chrom needs to be specified in order to apply Leave-one-chromosome-out on gene- or region-based tests")
        }else{
          chrom_v3 = getChromNumber(chrom)
        }
   }

 if(is.numeric(chrom_v3)){
  if(chrom_v3 >= 1 & chrom_v3 <= 22){

   chromList = c(1:22)
   chromList = chromList[which(chromList != chrom_v3)]
   modglmm = removeLOCOResult(chromList, modglmm)
   if(!is.null(modglmm$LOCOResult[[chrom_v3]])){
   modglmm$fitted.values = modglmm$LOCOResult[[chrom_v3]]$fitted.values
   modglmm$residuals = modglmm$LOCOResult[[chrom_v3]]$residuals
   modglmm$obj.noK = modglmm$LOCOResult[[chrom_v3]]$obj.noK
   if(is_Firth_beta){
     if(!is.null(modglmm$LOCOResult[[chrom_v3]]$offset)){
       modglmm$offset = modglmm$LOCOResult[[chrom_v3]]$offset
     }
   }
  }
   #modglmm$offset = modglmm$LOCOResult[[chrom_v3]]$linear.predictors -  modglmm$LOCOResult[[chrom_v3]]$coefficients[1]
   modglmm$LOCOResult[chrom_v3] = list(NULL)
  }else{ #if(chrom_v3 >= 1 & chrom_v3 <= 22){
     chromList = c(1:22)
     modglmm = removeLOCOResult(chromList, modglmm)
  }
 }else{
        chromList = c(1:22)
        modglmm = removeLOCOResult(chromList, modglmm)
  }
   gc()

 }

 modglmm$mu = as.vector(modglmm$fitted.values)
 modglmm$fitted.values <- NULL
 gc()
 tau = modglmm$theta
 N = length(modglmm$mu)
     if(modglmm$traitType == "binary"){
             modglmm$mu2 = (modglmm$mu) *(1-modglmm$mu)
	     modglmm$obj_cc = SKAT::SKAT_Null_Model(modglmm$y ~ modglmm$X-1, out_type="D", Adjustment = FALSE)
	     modglmm$obj_cc$mu = modglmm$mu
	     modglmm$obj_cc$res = modglmm$res
	     modglmm$obj_cc$pi_1 = modglmm$mu2
     }else if(modglmm$traitType == "quantitative"){
             modglmm$mu2 = (1/tau[1])*rep(1,N)
     }else if(modglmm$traitType == "survival"){
	      modglmm$mu2 = modglmm$mu
	      modglmm$obj.noK$XVX = modglmm$obj.noK$XVX_fg
	      modglmm$obj.noK$XVX_fg = NULL
	      modglmm$obj.noK$XXVX_inv = modglmm$obj.noK$XXVX_inv_fg
	      modglmm$obj.noK$XXVX_inv_fg = NULL
	      modglmm$obj.noK$XV = modglmm$obj.noK$XV_fg
	      modglmm$obj.noK$XV_fg = NULL
	      modglmm$obj.noK$XVX_inv_XV = modglmm$obj.noK$XVX_inv_XV_fg
	      modglmm$obj.noK$XVX_inv_XV_fg = NULL
	      modglmm$X = modglmm$obj.noK$X1_fg
	      modglmm$obj.noK$X1_fg = NULL
	      gc()
     }
 #if(FALSE){

 if(is_Firth_beta){
        if(modglmm$traitType == "binary"){
                if(is.null(modglmm$offset)){
                        #if(FALSE){
                        #print("WARNING. is_Firth_beta = TRUE, but offset was not computed in Step 1. Please re-run Step 1 using the more rencet version of SAIGE/SAIGE-GENE.")
                        cat("Applying is_Firth_beta = TRUE and note the results would be more accurate withe Step 1 results using the more rencet version of SAIGE/SAIGE-GENE. \n")

                        if(ncol(modglmm$X) == 1){
                                covoffset = rep(0,nrow(modglmm$X))
                        }else{
                                formulastr = paste0("y ~ ", paste(colnames(modglmm$X)[-1], collapse="+"))
                                modglmm$X = cbind(modglmm$X, as.vector(modglmm$y))
                                colnames(modglmm$X)[ncol(modglmm$X)] = "y"
                                modglmm$X = as.data.frame(modglmm$X)
                                formula.new = as.formula(formulastr)
                                modwitcov = glm(formula.new, data = modglmm$X, family = binomial)
                                modglmm$X = modglmm$X[,-ncol(modglmm$X)]
                                covoffset = as.matrix(modglmm$X[,-1]) %*%  modwitcov$coefficients[-1]
                                modglmm$X = as.matrix(modglmm$X)
                        }
                        modglmm$offset = covoffset
                        #}
                }
        }
  }


 if(is.null(modglmm$offset)){
         modglmm$offset = rep(0,nrow(modglmm$X))
  }


 if(is.null(modglmm$Sigma_iXXSigma_iX)){
        modglmm$Sigma_iXXSigma_iX = matrix(0, nrow=1, ncol=1)
 }
 return(modglmm)
}




ReadModel_subsample = function(GMMATmodelFile = "", chrom="", LOCO=TRUE, is_Firth_beta=FALSE, subSampleFile=""){
  # Check file existence
  Check_File_Exist(GMMATmodelFile, "GMMATmodelFile")
  if(!LOCO %in% c(TRUE, FALSE)){
    stop("LOCO should be TRUE or FALSE.")
  }

   if(!file.exists(subSampleFile)){
		stop("ERROR! ", subSampleFile, " does not exist\n")
   }else{
		cat("subSampleFile is specified. This option is used when any sample included in Step 1 but does not have dosages/genotypes for Step 2. Please make sure it contains one column of sample IDs that will be used for subsetting samples from the Step 1 results for Step 2 jobs\n")
	subSampleID0 = data.table::fread(subSampleFile, header = F, colClasses = c("character"), data.table=F)
   	subSampleID0 = subSampleID0[,1]
   }

  load(GMMATmodelFile)
  
  includeIndex = which(modglmm$sampleID %in% subSampleID0)
  if(length(includeIndex) == 0){
	stop("No samples analyzed in Step 1 are included in ", subSampleFile, "\n")
  }else{
	cat(length(includeIndex), " samples analyzed in Step 1 will be included for Step 2\n")

  	percNum = length(includeIndex) / (length(modglmm$sampleID))
	percNum = percNum * 100
	if(percNum < 90){
		cat("WARNING: ", percNum, "% samples analyzed in Step 1 will be included for Step 2. This proportion is small and using subsetted Step 1 results could lead to incorrect Step 2 results. Please re-run Step 1 for the sample subset\n")
	}	
  }	  


  modglmm$Y = NULL
  #modglmm$offset = modglmm$linear.predictors - modglmm$coefficients[1]
  modglmm$linear.predictors = NULL
  modglmm$coefficients = NULL
  modglmm$cov = NULL
  modglmm$obj.glm.null = NULL
  modglmm$fitted.values = modglmm$fitted.values[includeIndex]
  modglmm$Y = modglmm$Y[includeIndex]
  modglmm$residuals = modglmm$residuals[includeIndex]
  modglmm$sampleID = modglmm$sampleID[includeIndex]
  modglmm$y = modglmm$y[includeIndex]
  modglmm$X = modglmm$X[includeIndex,,drop=F]




  #rm(modglmm)
  gc()
  #traitType = modglmm$traitType
  #y = modglmm$y
  #X = modglmm$X
  #N = length(y)
  #tauVec = modglmm$theta
  #indChromCheck = FALSE
  chrom_v3=NULL

  if(!LOCO) {
    print("Leave-one-chromosome-out is not applied")
    if(modglmm$LOCO) {
      for (chr in 1:22) {
        modglmm$LOCOResult[chr] = list(NULL)
        cat("chromosome ", chr, " model results are removed to save memory\n")
        gc()
      }
    }
  }else{
    if (!modglmm$LOCO){
      stop("LOCO is TRUE but the null model file .rda does not contain LOCO results. In order to apply Leave-one-chromosome-out, please run Step 1 using LOCO. Otherwise, please set LOCO=FALSE in this step (Step 2).\n")
    }else{
        if(chrom == ""){
          stop("chrom needs to be specified in order to apply Leave-one-chromosome-out on gene- or region-based tests")
        }else{
          chrom_v3 = getChromNumber(chrom)
        }
   }

 if(is.numeric(chrom_v3)){
  if(chrom_v3 >= 1 & chrom_v3 <= 22){

   chromList = c(1:22)
   chromList = chromList[which(chromList != chrom_v3)]
   modglmm = removeLOCOResult(chromList, modglmm)
   if(!is.null(modglmm$LOCOResult[[chrom_v3]])){
   modglmm$fitted.values = modglmm$LOCOResult[[chrom_v3]]$fitted.values
   modglmm$residuals = modglmm$LOCOResult[[chrom_v3]]$residuals
   modglmm$residuals = modglmm$residuals[includeIndex]
   modglmm$obj.noK = modglmm$LOCOResult[[chrom_v3]]$obj.noK


   
   if(is_Firth_beta){
     if(!is.null(modglmm$LOCOResult[[chrom_v3]]$offset)){
       modglmm$offset = modglmm$LOCOResult[[chrom_v3]]$offset
     }
   }
  }
   #modglmm$offset = modglmm$LOCOResult[[chrom_v3]]$linear.predictors -  modglmm$LOCOResult[[chrom_v3]]$coefficients[1]
   modglmm$LOCOResult[chrom_v3] = list(NULL)
  }else{ #if(chrom_v3 >= 1 & chrom_v3 <= 22){
     chromList = c(1:22)
     modglmm = removeLOCOResult(chromList, modglmm)
  }
 }else{
        chromList = c(1:22)
        modglmm = removeLOCOResult(chromList, modglmm)
  }
   gc()

  }

 modglmm$mu = as.vector(modglmm$fitted.values)
 modglmm$mu = modglmm$mu[includeIndex]


 tau = modglmm$theta
 N = length(modglmm$mu)
     if(modglmm$traitType == "binary"){
             modglmm$mu2 = (modglmm$mu) *(1-modglmm$mu)
             modglmm$obj_cc = SKAT::SKAT_Null_Model(modglmm$y ~ modglmm$X-1, out_type="D", Adjustment = FALSE)
             modglmm$obj_cc$mu = modglmm$mu

	     if(!is.null(modglmm$res)){
		modglmm$res = modglmm$res[includeIndex]
	     }	     
             modglmm$obj_cc$res = modglmm$res
             modglmm$obj_cc$pi_1 = modglmm$mu2
    }else if(modglmm$traitType == "quantitative"){
             modglmm$mu2 = (1/tau[1])*rep(1,N)
     }
 #if(FALSE){

 if(is_Firth_beta){
        if(modglmm$traitType == "binary"){
                if(is.null(modglmm$offset)){
                        #if(FALSE){
                        #print("WARNING. is_Firth_beta = TRUE, but offset was not computed in Step 1. Please re-run Step 1 using the more rencet version of SAIGE/SAIGE-GENE.")
                        cat("Applying is_Firth_beta = TRUE and note the results would be more accurate withe Step 1 results using the more rencet version of SAIGE/SAIGE-GENE. \n")

                        if(ncol(modglmm$X) == 1){
                                covoffset = rep(0,nrow(modglmm$X))
                        }else{
                                formulastr = paste0("y ~ ", paste(colnames(modglmm$X)[-1], collapse="+"))
                                modglmm$X = cbind(modglmm$X, as.vector(modglmm$y))
                                colnames(modglmm$X)[ncol(modglmm$X)] = "y"
                                modglmm$X = as.data.frame(modglmm$X)
                                formula.new = as.formula(formulastr)
                                modwitcov = glm(formula.new, data = modglmm$X, family = binomial)
                                modglmm$X = modglmm$X[,-ncol(modglmm$X)]
                                covoffset = as.matrix(modglmm$X[,-1]) %*%  modwitcov$coefficients[-1]
                                modglmm$X = as.matrix(modglmm$X)
                        }
                        modglmm$offset = covoffset
                        #}
                }
        }
  }
 
 
obj.noK = ScoreTest_NULL_Model(modglmm$mu, modglmm$mu2, modglmm$y, modglmm$X)
modglmm$obj.noK = obj.noK
rm(obj.noK)	  


 if(is.null(modglmm$offset)){
         modglmm$offset = rep(0,nrow(modglmm$X))
  }else{
	  if(length(modglmm$offset) > N){
	    modglmm$offset = modglmm$offset[includeIndex]
	  }
  }	  


 if(is.null(modglmm$Sigma_iXXSigma_iX)){
        modglmm$Sigma_iXXSigma_iX = matrix(0, nrow=1, ncol=1)
 }else{
	modglmm$Sigma_iXXSigma_iX = modglmm$Sigma_iXXSigma_iX[includeIndex,, drop=F]
  }	 
 return(modglmm)
}





Get_Variance_Ratio<-function(varianceRatioFile, cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude, isGroupTest, isSparseGRM, useSparseGRMtoFitNULL){

    iscateVR = FALSE
    # check variance ratio
    if (!file.exists(varianceRatioFile)) {
	if(varianceRatioFile != ""){
	    stop("varianceRatioFile is specified but the file ", varianceRatioFile, " does not exist\n")
	}else{	
            cat("varianceRatioFile is not specified so variance ratio won't be used\n")
	    
	}
	if(isSparseGRM){
          cat("WARNING: Sparse GRM is specified. Please make sure the null model was fit using the same sparse GRM in Step 1.\n")
          ratioVec_sparse = c(1)
	  ratioVec_null = c(-1)
	}else{
          ratioVec_sparse = c(-1)		
	  ratioVec_null = c(1)
	}
    }else{
        varRatioData = data.frame(data.table:::fread(varianceRatioFile, header = F, stringsAsFactors = FALSE))
	if(ncol(varRatioData) == 3){
	    spindex = which(varRatioData[,2] == "sparse")
	    if(length(spindex) > 0){
	        ratioVec_sparse = varRatioData[which(varRatioData[,2] == "sparse"),1]
		ratioVec_sparse = as.numeric(ratioVec_sparse)
		#if(!isSparseGRM & sum(ratioVec_sparse > 1.0001 | ratioVec_sparse < 0.9999) > 0){
		#       	stop("sparse GRM is not specified but it was used for estimating variance ratios in Step 1. Please specify --sparseGRMFile and --sparseGRMSampleIDFile\n")
		#}	
	    }else{    
		ratioVec_sparse = c(-1)
	    	if(isSparseGRM){
			stop("sparse GRM is specified but the variance ratio for sparse GRM was not estimatedin Step 1. Pleae remove --sparseGRMFile and --sparseGRMSampleIDFile\n")
		}	
	    }
	    ratioVec_null = varRatioData[which(varRatioData[,2] == "null"),1]
	    ratioVec_null = as.numeric(ratioVec_null)
            cat("variance Ratio null is ", ratioVec_null, "\n")
	    if(length(ratioVec_null) > 1){
		iscateVR = TRUE
		nrv = length(ratioVec_null)
	    }else{
		nrv = 1
	    }
	}else{
	    cat("Variance ratios were estimated with version < 1.0.6\n")	
	    if(isSparseGRM){
	        ratioVec_sparse = varRatioData[,1]
	        ratioVec_sparse = as.numeric(ratioVec_sparse)
        	cat("variance Ratio is ", ratioVec_sparse, "\n")
		ratioVec_null = rep(-1, length(varRatioData[,1]))
		cat("WARNING: Sparse GRM is specified and the variance ratio(s) were specified. Please make sure the variance ratios were estimated using a full GRM and a sparse GRM.")
	    }else{
	        ratioVec_null = varRatioData[,1]
	        ratioVec_null = as.numeric(ratioVec_null)
                cat("variance Ratio is ", ratioVec_null, "\n")
		cat("WARNING: Sparse GRM is not specified and the variance ratio(s) were specified. Please make sure that in Step 1, 1. the null model was fit using a full GRM (--useSparseGRMtoFitNULL=FALSE) and the variacne ratio was NOT estimated with the sparse GRM (useSparseGRMforVarRatio=FALSE) or 2. the null model was fit using a sparse GRM (--useSparseGRMtoFitNULL=TRUE) and the variacne ratio was estiamted with the sparse GRM and null --skipVarianceRatioEstimation=FALSE\n")
		ratioVec_sparse = c(-1)
	    }
	    if(length(varRatioData[,1]) > 1){
		iscateVR = TRUE
		nrv = length(varRatioData[,1])
	    }else{
		nrv = 1
	    }
	}


	if(iscateVR){
            ln = length(cateVarRatioMinMACVecExclude)
            hn = length(cateVarRatioMaxMACVecInclude)

            if (nrv != ln) {
                stop("ERROR! The number of variance ratios are different from the length of cateVarRatioMinMACVecExclude\n")
            }
            if (ln != (hn + 1)) {
                stop("ERROR! The length of cateVarRatioMaxMACVecInclude does not match with the lenght of cateVarRatioMinMACVecExclude (-1)\n")
            }
        }
	

    }

    print("ratioVec_sparse")
    print(ratioVec_sparse)
    print("ratioVec_null")
    print(ratioVec_null)
    return(list(ratioVec_sparse = ratioVec_sparse, ratioVec_null = ratioVec_null))
}

ReadModel_multiTrait <- function(GMMATmodelFileList = "", chrom = "", LOCO = TRUE, is_Firth_beta = FALSE, is_EmpSPA = FALSE, espa_nt = 9999, espa_range = c(-20, 20)) {
  # Check file existence
  GMMATmodelFileVec <- unlist(strsplit(GMMATmodelFileList, split = ","))
  modglmmList <- list()
  sampleIDList <- list()
  if (require("furrr")) {
    future::plan("multicore")
    modglmmList <- furrr::future_map(GMMATmodelFileVec, function(GMMATmodelFile) {
      ReadModel(GMMATmodelFile, chrom, LOCO, is_Firth_beta, is_EmpSPA, espa_nt, espa_range)
    })
    future::plan("sequential")
  } else {
    for (gm in 1:length(GMMATmodelFileVec)) {
      GMMATmodelFile <- GMMATmodelFileVec[gm]
      modglmmList[[gm]] <- ReadModel(GMMATmodelFile, chrom, LOCO, is_Firth_beta, is_EmpSPA, espa_nt, espa_range)
    }
  }
  union_vector <- c()

  for (gm in 1:length(GMMATmodelFileVec)) {
    vec <- modglmmList[[gm]]$sampleID
    union_vector <- unique(c(union_vector, vec))
  }
  nx = 0
  upperx = 0
  for (gm in 1:length(GMMATmodelFileVec)) {
    vec <- modglmmList[[gm]]$sampleID 
    modglmmList[[gm]]$sampleIndices <- match(vec, union_vector)
    modglmmList[[gm]]$sampleID = NULL
    nx = nx + ncol(modglmmList[[gm]]$X)
    upperx = max(upperx, ncol(modglmmList[[gm]]$X))
  cat("upperx ", upperx, "\n")
  #print(ncol(modglmmList[[gm]]$obj.noK$X1)+1)
  print("ncol(modglmmList[[gm]]$X)")
  print(ncol(modglmmList[[gm]]$X))
  }

  cat("upperx ", upperx, "\n")

  modglmmList[[1]]$sampleID = union_vector
  modglmmList[[1]]$nx = nx
  modglmmList[[1]]$upperx = upperx
  gc()
  return(modglmmList)
}

Get_Variance_Ratio_multiTrait <- function(varianceRatioFileList, cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude, isGroupTest, isSparseGRM, useSparseGRMtoFitNULL) {
  varianceRatioFileVec <- unlist(strsplit(varianceRatioFileList, split = ","))
  # varianceRatioList = list()
  ratioVec_sparse <- NULL
  ratioVec_null <- NULL
  #ratioVec_null_sample <- NULL
  ratioVec_null_noXadj <- NULL
  #ratioVec_null_eg <- NULL
  #ratioVec_sparse_eg <- NULL
  for (vrl in 1:length(varianceRatioFileVec)) {
    varianceRatioFile <- varianceRatioFileVec[vrl]
    ratioVecList <- Get_Variance_Ratio(varianceRatioFile, cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude, isGroupTest, isSparseGRM)
    print(ratioVecList)
    ratioVec_sparse <- cbind(ratioVec_sparse, ratioVecList$ratioVec_sparse)
    ratioVec_null <- cbind(ratioVec_null, ratioVecList$ratioVec_null)
    #ratioVec_null_sample <- cbind(ratioVec_null_sample, ratioVecList$ratioVec_null_sample)
    ratioVec_null_noXadj <- cbind(ratioVec_null_noXadj, ratioVecList$ratioVec_null_noXadj)
    #ratioVec_null_eg <- cbind(ratioVec_null_eg, ratioVecList$ratioVec_null_eg)
    #ratioVec_sparse_eg <- cbind(ratioVec_sparse_eg, ratioVecList$ratioVec_sparse_eg)
  }
  #return(list(ratioVec_sparse = ratioVec_sparse, ratioVec_null = ratioVec_null, ratioVec_null_sample = ratioVec_null_sample, ratioVec_null_noXadj = ratioVec_null_noXadj, ratioVec_null_eg = ratioVec_null_eg, ratioVec_sparse_eg = ratioVec_sparse_eg))
  return(list(ratioVec_sparse = ratioVec_sparse, ratioVec_null = ratioVec_null))
}
