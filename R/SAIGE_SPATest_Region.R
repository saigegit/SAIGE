setSparseSigma = function(sparseSigmaFile){

  Check_File_Exist(sparseSigmaFile, "sparseSigmaFile")
  sparseSigma = Matrix:::readMM(sparseSigmaFile)
  locations = rbind(sparseSigma@i, sparseSigma@j)
  values = sparseSigma@x
  nSubj = dim(sparseSigma)[1]
  sigmaMatListR = list(locations = locations,
                     values = values,
                     nSubj = nSubj)
  return(sigmaMatListR)	
}	


setSparseSigma_new = function(sparseGRMFile, sparseGRMSampleIDFile, relatednessCutoff, sampleIDInModel, tauVec, W, traitType){

  Check_File_Exist(sparseGRMFile, "sparseGRMFile")
  Check_File_Exist(sparseGRMSampleIDFile, "sparseGRMSampleIDFile")

  sparseGRM = Matrix:::readMM(sparseGRMFile)
  sparseGRMSampleID = data.frame(data.table:::fread(sparseGRMSampleIDFile, header=F, stringsAsFactors=FALSE, colClasses=c("character")))
  colnames(sparseGRMSampleID) = "sampleID"
  sparseGRMSampleID$IndexGRM = c(1:nrow(sparseGRMSampleID))
  cat("length(sparseGRMSampleID$IndexGRM): ", length(sparseGRMSampleID$IndexGRM), "\n")
  cat("nrow(sparseGRMSampleID): ", nrow(sparseGRMSampleID), "\n")
  sampleInModel = NULL
  sampleInModel$IID = sampleIDInModel
  #rm(sampleIDInModel)
  sampleInModel = data.frame(sampleInModel)
  sampleInModel$IndexInModel = seq(1,length(sampleInModel$IID), by=1)
  cat(nrow(sampleInModel), " samples have been used to fit the glmm null model\n")
  mergeID = merge(sampleInModel, sparseGRMSampleID, by.x="IID", by.y = "sampleID")
  mergeID = mergeID[with(mergeID, order(IndexInModel)), ]
  print(dim(mergeID))
  print(head(mergeID))
  indexIDofGRM=mergeID$IndexGRM
  sparseGRM = sparseGRM[indexIDofGRM, indexIDofGRM]
  if(length(indexIDofGRM) < nrow(sampleInModel)){
    stop(nrow(sampleInModel)-length(indexIDofGRM), " samples were not found in the sparse GRM\n")
  }else{
    print("Subsetting GRM")
  }
  removeIndex = which(sparseGRM@x < relatednessCutoff)
  if(length(removeIndex) > 0){
	cat("Removing ", length(removeIndex), " elements in the sparse GRM < ", relatednessCutoff, ".\n")  
  	sparseGRM@x = sparseGRM@x[-removeIndex]
	sparseGRM@i = sparseGRM@i[-removeIndex]
	sparseGRM@j = sparseGRM@j[-removeIndex]
  }

  sparseGRM@x = sparseGRM@x * tauVec[2]  
  #sparseSigma = sparseGRM * tauVec[2]
  if(traitType == "binary"){
	sparseGRM@x[which(sparseGRM@i == sparseGRM@j)] = sparseGRM@x[which(sparseGRM@i == sparseGRM@j)] + 1/W 
   #diag(sparseSigma) = W + diag(sparseSigma)
  }else if(traitType == "quantitative"){
	sparseGRM@x[which(sparseGRM@i == sparseGRM@j)] = tauVec[1] + sparseGRM@x[which(sparseGRM@i == sparseGRM@j)]
  }
   #diag(sparseSigma) = tauVec[1] + diag(sparseSigma)

  locations = rbind(sparseGRM@i, sparseGRM@j)
  values = sparseGRM@x
  nSubj = dim(sparseGRM)[1]
  sigmaMatListR = list(locations = locations,
                     values = values,
                     nSubj = nSubj)
  return(sigmaMatListR)
}



                        #DosageCutoff_for_UltraRarePresence, 
			#method_to_CollapseUltraRare,
##working
SAIGE.Region = function(mu,
			OutputFile,
                        MACCutoff_to_CollapseUltraRare,
			groupFile, 
			annolist, 
			maxMAFlist,
		        markers_per_chunk_in_groupTest,	
			genoType, 
			markerInfo,
			traitType,
			isImputation,
			isCondition,
			weight_cond,
			groups_per_chunk,
			r.corr,
			isOverWriteOutput,
			is_single_in_groupTest,
			is_no_weight_in_groupTest,
			is_output_markerList_in_groupTest,
			chrom){
  OutputFileIndex = NULL	
  if(is.null(OutputFileIndex))
    OutputFileIndex = paste0(OutputFile, ".index")

  outList = checkOutputFile(OutputFile, OutputFileIndex, "Region", 1, isOverWriteOutput) # Check 'Util.R'

  indexChunk = outList$indexChunk
  Start = outList$Start
  End = outList$End	

  cat("Start ", Start, "\n")
  cat("End ", End, "\n")


  if(End)
  {
    message = paste0("The analysis has been completed in earlier analysis. Results are saved in '", OutputFile, "'. ",
                     "If you want to change parameters and restart the analysis, please use another 'OutputFile'.")
    return(message)
  }

 isappend=FALSE
 if(!Start){ 
  isappend=TRUE
 }

  n = length(mu) #sample size 


  ## annotation in region
  ##need to revise
  if(r.corr==0){
    out.method = SKAT:::SKAT_Check_Method(method="optimal.adj", r.corr=0)
    method=out.method$method
    r.corr=out.method$r.corr
    cat("SKAT-O test will be performed. P-values for BURDEN and SKAT tests will also be output\n")
    regionTestType = "SKAT-O"
    is_single_in_groupTest = TRUE
    #cat("is_single_in_groupTest = TRUE. Single-variant assoc tests results will be output\n")
  }else if(r.corr == 1){	  
    method = NULL
    cat("BURDEN test will be performed\n")
    regionTestType = "BURDEN"

    #output the result from Rcpp
    cat("isappend ", isappend, "\n")
    isOpenOutFile = openOutfile(traitType, isappend)
    if(!isOpenOutFile){
	stop("Output file ", OutputFile, " can't be opened\n")
    }

  }else{
    stop("r.corr needs to be either 1 (BURDEN test) or 0 (SKAT-O test)\n")
  }	  

  if(is_single_in_groupTest){
      cat("is_single_in_groupTest = TRUE. Single-variant assoc tests results will be output\n")
    isOpenOutFile_singleinGroup = openOutfile_singleinGroup(traitType, isImputation, isappend)
    if(!isOpenOutFile_singleinGroup){
        stop("Output file ", OutputFile, ".singleAssoc.txt can't be opened\n")
    }

  }else{
      cat("is_single_in_groupTest = FALSE. Single-variant assoc tests results will not be output\n")
  }


##check group file
  region_list = checkGroupFile(groupFile)
  nRegions = region_list$nRegions
  is_weight_included = region_list$is_weight_included
  if(is_no_weight_in_groupTest & is_weight_included){
    stop("is_no_weight_in_groupTest = TRUE but weights are found in the group file.\n")
  }

  if(is_no_weight_in_groupTest){
     cat("No weights are used in the group test\n")
  }

  if(is_weight_included){
    nline_per_gene = 3
  }else{
    nline_per_gene = 2
  }	  

  gf = file(groupFile, "r")

  cat("indexChunk is ", indexChunk, "\n")

  skipline = indexChunk*nline_per_gene
  if(indexChunk > 0 & indexChunk < nRegions){
    for(k in 1:skipline){
    marker_group_line_temp = readLines(gf, n = 1) 
    rm(marker_group_line_temp)
    }
  }

  if(regionTestType != "BURDEN"){
  	P1Mat = matrix(0, markers_per_chunk_in_groupTest, n)
  	P2Mat = matrix(0, n, markers_per_chunk_in_groupTest)
  }else{
	P1Mat = matrix(0, 1, 1)
	P2Mat = matrix(0, 1, 1)  
  }	  

  chrom1 = "FakeCHR";

  gc()
  num_region = 0
  mth = 0

  numberRegionsInChunk = 0
  cat("indexChunk ", indexChunk, "\n")
  cat("nRegions ", nRegions, "\n")
  pval.Region.all = NULL
  OutList.all = NULL
  Output_MarkerList.all = NULL
  cth_chunk_to_output=1

  i = indexChunk+1
  while(i <= nRegions){
  #for(i in (indexChunk+1):nRegions){


   if(mth ==  numberRegionsInChunk){
      if(i + groups_per_chunk > nRegions){
  	      nregions_ro_read = nRegions - i + 1	      
      }else{
	      nregions_ro_read = groups_per_chunk
      }
      nlinetoread = nregions_ro_read * nline_per_gene
      marker_group_line = readLines(gf, n = nlinetoread)
      RegionList = SAIGE.getRegionList_new(marker_group_line, nline_per_gene, annolist, markerInfo, chrom)
      cat("Read in ", nregions_ro_read, " region(s) from the group file.\n")
      mth = 0
      #numberRegionsInChunk = length(RegionList)
      numberRegionsInChunk = nregions_ro_read
    }

   mth = mth + 1
   if(!is.null(RegionList)){
    pval.Region = NULL
    region = RegionList[[mth]]

    annolist = region$annoVec 
    regionName = names(RegionList)[mth]
    i = i + 1
    if(!is.null(region$SNP) & length(annolist) > 0){

      SNP = region$SNP
      if(genoType == "vcf"){
        SNPlist = paste(c(regionName, SNP), collapse = "\t") 
        set_iterator_inVcf(SNPlist, chrom1, 1, 200000000)
        isVcfEnd =  check_Vcf_end()
    	if(!isVcfEnd){
		region$genoIndex = rep("0", length(SNP))
		region$genoIndex_prev = rep("0", length(SNP))
    	}else{
        	warning("No markers in region ", regionName, " are found in the VCF file")
		next
    	}	
      }

      if(is_no_weight_in_groupTest){
	WEIGHT = rep(1, length(SNP))
      }else{	      
        WEIGHT = region$WEIGHT
      }

      annoIndicatorMat = region$annoIndicatorMat

      #chrom = region$chrom

      print(paste0("Analyzing Region ", regionName, " (",i-1,"/",nRegions,")."))
      #tp1 = proc.time()
      #gc()

      #time_mainRegionInCPP = system.time({outList = mainRegionInCPP(genoType, genoIndex, annoIndicatorMat, maxMAFlist, OutputFile, traitType, n, P1Mat, P2Mat, regionTestType, isImputation, WEIGHT, weight_cond, is_single_in_groupTest, is_output_markerList_in_groupTest)})
      outList = mainRegionInCPP(genoType, region$genoIndex_prev, region$genoIndex, annoIndicatorMat, maxMAFlist, OutputFile, traitType, n, P1Mat, P2Mat, regionTestType, isImputation, WEIGHT, weight_cond, is_single_in_groupTest, is_output_markerList_in_groupTest, annolist, regionName)
#print("time_mainRegionInCPP")
#print(time_mainRegionInCPP)
      rm(region)
      #rm(genoIndex)
      #gc()
    #tb0 = proc.time()
if(is_single_in_groupTest){
      #OutList = as.data.frame(outList$OUT_DF)
      noNAIndices = which(!is.na(outList$pvalVec))
      if(sum(WEIGHT) > 0){
        AnnoWeights = c(WEIGHT, rep(1, outList$numofUR))
      }
}



annoMAFIndicatorMat = outList$annoMAFIndicatorMat

if((sum(outList$NumUltraRare_GroupVec) + sum(outList$NumRare_GroupVec)) > 0){

if(regionTestType != "BURDEN"){
     #is_single_in_groupTest is TRUE
 	#ta0 = proc.time()
       if(traitType == "binary"){	     
         outList$gyVec = outList$gyVec[noNAIndices]
       }

       if(isCondition){
         outList$VarMatAdjCond = outList$VarMatAdjCond[noNAIndices,noNAIndices]
         outList$TstatAdjCond = outList$TstatAdjCond[noNAIndices]
         outList$G1tilde_P_G2tilde_Weighted_Mat = outList$G1tilde_P_G2tilde_Weighted_Mat[noNAIndices,,drop=F]
         weightMat_G2_G2 = outList$G2_Weight_cond %*% t(outList$G2_Weight_cond)
       }	  
        
       ### Get annotation maf indicators

       ### 3. Adjust for saddlepoint approximation
       #notNAindice = which(!is.na(outList$TstatVec_flip))
       #print("notNAindice")
       #print(notNAindice)
       StatVec = outList$TstatVec_flip[noNAIndices]
       VarSVec = diag(outList$VarMat)
       VarSVec = VarSVec[!is.na(VarSVec)]
       adjPVec = outList$pvalVec[!is.na(outList$pvalVec)]
		
       #varTestedIndices = which(apply(annoMAFIndicatorMat, 1, isContainValue, val=1))
       #print("varTestedIndices")
       #print(varTestedIndices)
       #varTestedIndices = which(rowSums(annoMAFIndicatorMat) > 0)
       #annoMAFIndicatorMat = annoMAFIndicatorMat[varTestedIndices, , drop=F]
       #MAFVec = outList$MAFVec[varTestedIndices]
       annoMAFIndicatorMat = annoMAFIndicatorMat[noNAIndices, , drop=F]
       MAFVec = outList$MAFVec[noNAIndices]
       #AnnoWeights = dbeta(MAFVec,1,25)

       if(sum(WEIGHT) > 0){
 	 #AnnoWeights = AnnoWeights[varTestedIndices] 
 	 AnnoWeights = AnnoWeights[noNAIndices] 
       }else{
         AnnoWeights = dbeta(MAFVec,1,25)
       }
       weightMat = AnnoWeights %*% t(AnnoWeights)


       if(isCondition){    
    	 weightMat_G1_G2 = AnnoWeights %*% t(outList$G2_Weight_cond)
       }	
       wStatVec = StatVec * AnnoWeights

      wadjVarSMat = outList$VarMat * weightMat

	if(isCondition){
	  wStatVec_cond = wStatVec - outList$TstatAdjCond
          wadjVarSMat_cond = wadjVarSMat - outList$VarMatAdjCond	
      }

    #gc()


    annoMAFIndVec = c()
    for(j in 1:length(annolist)){
	AnnoName = annolist[j]
	maxMAF0 = outList$q_maf_for_annoVec[j]
	isPolyRegion = TRUE
	for(m in 1:length(maxMAFlist)){
		jm = (j-1)*(length(maxMAFlist)) + m
		maxMAFName = maxMAFlist[m]
	    if(m <= maxMAF0){
		tempPos = which(annoMAFIndicatorMat[,jm] == 1)
	       if(length(tempPos) > 0){
	       isPolyRegion = TRUE
		annoMAFIndVec = c(annoMAFIndVec, jm)
		Phi = wadjVarSMat[tempPos, tempPos, drop=F]
		Score = wStatVec[tempPos]
		if(traitType == "binary"){
			p.new = adjPVec[tempPos]
			g.sum = outList$genoSumMat[,jm]
			q.sum<-sum(outList$gyVec[tempPos] * AnnoWeights[tempPos])
			mu.a = mu
			re_phi = get_newPhi_scaleFactor(q.sum, mu.a, g.sum, p.new, Score, Phi, regionTestType)
		        Phi = re_phi$val
                }
		groupOutList = get_SKAT_pvalue(Score, Phi, r.corr, regionTestType)

		resultDF = data.frame(Region = regionName,
                                                    Group = AnnoName,
                                                    max_MAF = maxMAFName,
                                                    Pvalue = groupOutList$Pvalue_SKATO,
                                                    Pvalue_Burden = groupOutList$Pvalue_Burden,
                                                    Pvalue_SKAT = groupOutList$Pvalue_SKAT,
                                                    BETA_Burden = groupOutList$BETA_Burden,
                                                    SE_Burden = groupOutList$SE_Burden)
	      if(isCondition){
		if(traitType == "binary"){
			G1tilde_P_G2tilde_Mat_scaled = t(t((outList$G1tilde_P_G2tilde_Weighted_Mat[tempPos,,drop=F]) * sqrt(as.vector(re_phi$scaleFactor))) * sqrt(as.vector(outList$scalefactor_G2_cond)))
#t(t(b * sqrt(a1)) * sqrt(a2))
		        adjCondTemp = G1tilde_P_G2tilde_Mat_scaled %*% outList$VarInvMat_G2_cond_scaled	
			VarMatAdjCond = adjCondTemp %*% t(G1tilde_P_G2tilde_Mat_scaled)
			TstatAdjCond = adjCondTemp %*% (outList$Tstat_G2_cond * outList$G2_Weight_cond)
			Phi_cond = re_phi$val - diag(VarMatAdjCond)
			Score_cond = Score - TstatAdjCond
			
		}else{
			Score_cond = wStatVec_cond[tempPos]
			Phi_cond = wadjVarSMat_cond[tempPos, tempPos]
		}

		groupOutList_cond = get_SKAT_pvalue(Score_cond, Phi_cond, r.corr, regionTestType)

		resultDF$Pvalue_cond = groupOutList_cond$Pvalue_SKATO
		resultDF$Pvalue_Burden_cond = groupOutList_cond$Pvalue_Burden
		resultDF$Pvalue_SKAT_cond = groupOutList_cond$Pvalue_SKAT
		resultDF$BETA_Burden_cond = groupOutList_cond$BETA_Burden
		resultDF$SE_Burden_cond = groupOutList_cond$SE_Burden
	     }#if(isCondition){
		pval.Region = rbind.data.frame(pval.Region, resultDF)

	   }else{#if(length(tempPos) > 0){
		isPolyRegion = FALSE
	   }	

	}else{ #if(m <= maxMAF0){
	   if(isPolyRegion){
		annoMAFIndVec = c(annoMAFIndVec, jm)   
		resultDF$Region = regionName
		resultDF$Group = AnnoName
		resultDF$max_MAF = maxMAFName	
		pval.Region = rbind.data.frame(pval.Region, resultDF)
	   }			   
	}#if(m <= maxMAF0){

    }#for(m in 1:length(maxMAFlist)){
}#for(j in 1:length(annolist)){


gc()
}



if(regionTestType != "BURDEN"){

    if(length(annoMAFIndVec) > 0){
      pval.Region$MAC = outList$MAC_GroupVec[annoMAFIndVec]
      if(traitType == "binary"){
    	pval.Region$MAC_case = outList$MACCase_GroupVec[annoMAFIndVec]
    	pval.Region$MAC_control = outList$MACCtrl_GroupVec[annoMAFIndVec]
      }
      pval.Region$Number_rare = outList$NumRare_GroupVec[annoMAFIndVec]
      pval.Region$Number_ultra_rare = outList$NumUltraRare_GroupVec[annoMAFIndVec]
    }

if(length(annolist) > 1 | length(maxMAFlist) > 1){

   cctpval_Burden = get_CCT_pvalue(pval.Region$Pvalue_Burden)
   if(regionTestType != "BURDEN"){
     cctpval = get_CCT_pvalue(pval.Region$Pvalue)
     cctpval_SKAT = get_CCT_pvalue(pval.Region$Pvalue_SKAT)
     cctVec = c(regionName, "Cauchy", NA, cctpval, cctpval_Burden, cctpval_SKAT, NA, NA)
   }else{
	cctVec = c(regionName, "Cauchy", NA, cctpval_Burden, NA, NA)
   } 	   
   if(isCondition){
   	cctpval_Burden = get_CCT_pvalue(pval.Region$Pvalue_Burden_cond)

   	if(regionTestType != "BURDEN"){
	  cctpval = get_CCT_pvalue(pval.Region$Pvalue_cond)
	  cctpval_SKAT = get_CCT_pvalue(pval.Region$Pvalue_SKAT_cond)
	  cctVec = c(cctVec, cctpval, cctpval_Burden, cctpval_SKAT, NA, NA)
	}else{
	  cctVec = c(cctVec, cctpval_Burden, NA, NA)

	}	
   }
   cctVec = c(cctVec, NA)
   if(traitType == "binary"){
     cctVec = c(cctVec, NA, NA)
   }
   cctVec = c(cctVec, NA, NA)

   pval.Region = rbind(pval.Region, cctVec)
}

#ta1 = proc.time()
#print("ta0 - tb0")
#print(ta0 - tb0)
#print("ta1 - ta0")
#print(ta1 - ta0)
}#if(regionTestType != "BURDEN"){

  Output_MarkerList = NULL
  if(is_output_markerList_in_groupTest){
    for(j in 1:length(annolist)){
        AnnoName = annolist[j]
        for(m in 1:length(maxMAFlist)){
                jm = (j-1)*(length(maxMAFlist)) + m
                maxMAFName = maxMAFlist[m]
                tempPos = which(outList$annoMAFIndicatorMat[,jm] > 0)
                marker_rare_pos = which(outList$markerIndcatorVec == 1)
                marker_ultrarare_pos = which(outList$markerIndcatorVec == 2)
                if(length(tempPos) > 0){
                        if(length(marker_rare_pos) > 0){
				markerind_b = which(marker_rare_pos %in% tempPos)
				if(length(markerind_b) > 0){
					markerind = marker_rare_pos[markerind_b]
					SNPlist_rare = paste(SNP[markerind], collapse=",")
				}else{
					SNPlist_rare = ""
				}	
                        }else{
                                SNPlist_rare = ""
                        }
                        if(length(marker_ultrarare_pos) > 0){
				markerind_UR_b = which(marker_ultrarare_pos %in% tempPos)
				if(length(markerind_UR_b) > 0){
					markerindUR = marker_ultrarare_pos[markerind_UR_b]
					SNPlist_Ultra_rare = paste(SNP[markerindUR], collapse=",")
				}else{
					SNPlist_Ultra_rare = ""
				}	
                        }else{
                                SNPlist_Ultra_rare = ""
                        }
                        Output_MarkerList = rbind(Output_MarkerList, c(regionName, AnnoName, maxMAFName, SNPlist_rare, SNPlist_Ultra_rare))
                }

        }
    }
  }

#ta2 = proc.time()	
#print("ta2 - ta1")
#print(ta2 - ta1)
   if(is_output_markerList_in_groupTest){
	colnames(Output_MarkerList) = c("Region", "Group", "max_MAF", "Rare_Variants", "Ultra_Rare_Variants")
	Output_MarkerList.all = rbind(Output_MarkerList.all, Output_MarkerList)   
   }else{
	Output_MarkerList.all = NULL
   }	   

  indexChunk = i
  #Start = (i==1)
  Start = (cth_chunk_to_output==1)
  End = (i==nRegions)
  AnalysisType = "Region"
  nEachChunk = 1

if(regionTestType != "BURDEN"){  
    pval.Region.all = rbind(pval.Region.all, pval.Region)
}


# output

if(mth ==  numberRegionsInChunk){
  message1 = "This is the output index file for SAIGE package to record the end point in case users want to restart the analysis. Please do not modify this file."
  message2 = paste("This is a", AnalysisType, "level analysis.")
  message3 = paste("nEachChunk =", nEachChunk)
  message4 = paste("Have completed the analysis of chunk", indexChunk)
  message5 = "Have completed the analyses of all chunks."
  #n1 = length(Output)
  #n2 = length(OutputFile)
  print("write to output")
  #cat("n1 is ", n1, "\n")
  #cat("n2 is ", n2, "\n")
  if(regionTestType != "BURDEN"){  
      if(Start){
        if(!is.null(pval.Region.all)){
          fwrite(pval.Region.all, OutputFile, quote = F, sep = "\t", append = F, col.names = T, row.names = F, na="NA")
        }
      }else{
        if(!is.null(pval.Region.all)){
          fwrite(pval.Region.all, OutputFile, quote = F, sep = "\t", append = T, col.names = F, row.names = F, na="NA")
        }
        #write.table(Output, OutputFile, quote = F, sep = "\t", append = T, col.names = F, row.names = F)
      }
  }
      if(is_output_markerList_in_groupTest){ 
        if(Start){
        if(!is.null(Output_MarkerList.all)){
          fwrite(Output_MarkerList.all, paste0(OutputFile, ".markerList.txt"), quote = F, sep = "\t", append = F, col.names = T, row.names = F, na="NA")
        }
      }else{
        if(!is.null(Output_MarkerList.all)){
          fwrite(Output_MarkerList.all, paste0(OutputFile, ".markerList.txt"), quote = F, sep = "\t", append = T, col.names = F, row.names = F, na="NA")
        }
        #write.table(Output, OutputFile, quote = F, sep = "\t", append = T, col.names = F, row.names = F)
      }
    }



#if(FALSE){
  #print("write Output 2")
  if(Start){
    write.table(c(message1, message2, message3), OutputFileIndex,
                quote = F, sep = "\t", append = F, col.names = F, row.names = F)
  }
  write.table(message4, OutputFileIndex, quote = F, sep = "\t", append = T, col.names = F, row.names = F)

  if(End){
    write.table(message5, OutputFileIndex, quote = F, sep = "\t", append = T, col.names = F, row.names = F)
  }
#}#if(FALSE)

pval.Region.all = NULL
OutList.all = NULL
Output_MarkerList.all = NULL
cth_chunk_to_output = cth_chunk_to_output + 1
gc()
}

#tp2 = proc.time()
#print("tp2 - tp1")
#print(tp2 - tp1)
#print("tp2 - ta2")
#print(tp2 - ta2)
   #if(is_single_in_groupTest){
   #  rm(OutList)
   # }

   if(is_output_markerList_in_groupTest){
     rm(Output_MarkerList)
   } 	   

   rm(outList)
   rm(pval.Region)
   if(regionTestType != "BURDEN"){
     rm(resultDF)
   } 
   gc()
  
 }#if(length(noNAIndices) > 0){ 
  }else{#if(!is.null(region)){
    cat(regionName, " is empty.\n")
  }

   }else{#if(!is.null(RegionList)){
     cat("The chunk is empty\n")	   
     #mth = 0
     mth = numberRegionsInChunk
     i = i + numberRegionsInChunk
     pval.Region = NULL
   }
	   
}

  message = paste0("Analysis done! The results have been saved to '", OutputFile,"' and '",
                   paste0(OutputFile, ".markerInfo"),"'.")

  message = "Analysis done!"
  message = paste0(message, " The set-based tests results have been saved to '", OutputFile, "'.")
  if(is_output_markerList_in_groupTest){
	message = paste0(message, " The marker lists have been saved to '", OutputFile, ".markerList.txt'.")
  }	  
  if(is_single_in_groupTest){
  	message = paste0(message, " The single-variant association tests results have been saved to '", OutputFile, ".singleAssoc.txt'.")

  }

  return(message)

}



SAIGE.getRegionList_new = function(marker_group_line,
			nline_per_gene,	   
                         annoVec, #c("lof","lof;missense"
                         markerInfo, 
			 chrom="")
{

  # read group file
  ngroup<-length(marker_group_line)/nline_per_gene
  #cat("ngroup is ", ngroup, "\n")
  RegionData = NULL
  geneList = c() 
  for(i in 1:ngroup){
	  marker_group_line_list = strsplit(marker_group_line[1+(i-1)*nline_per_gene], split=c(" +", "\t"))[[1]]
          gene=marker_group_line_list[1]
          var=marker_group_line_list[3:length(marker_group_line_list)]
          marker_group_line_list_anno = strsplit(marker_group_line[2+(i-1)*nline_per_gene], split=c(" +", "\t"))[[1]]
	  anno=marker_group_line_list_anno[3:length(marker_group_line_list_anno)]
	  if(nline_per_gene == 3){
              marker_group_line_list_weight = strsplit(marker_group_line[3+(i-1)*nline_per_gene], split=c(" +", "\t"))[[1]]
              weight=marker_group_line_list_weight[3:length(marker_group_line_list_weight)]
	      RegionData = rbind(RegionData, cbind(rep(gene, length(var)), var, anno, weight))
          }else if(nline_per_gene == 2){
	      RegionData = rbind(RegionData, cbind(rep(gene, length(var)), var, anno))
	  }

	  if(gene %in% geneList){
		stop(gene, " is duplicated in the group File\n")
          }else{		  
	  	geneList = c(geneList, gene)
	  }
  }
    if(nline_per_gene == 2){
    	colnames(RegionData) = c("REGION", "SNP", "ANNO")
    }else if(nline_per_gene == 3){
	colnames(RegionData) = c("REGION", "SNP", "ANNO", "WEIGHT")
    }	    
    RegionData = as.data.frame(RegionData)
    setDT(RegionData)
    uRegion0 = unique(RegionData$REGION)    
    if(chrom != "" & is.null(markerInfo)){
      RegionData[, c("chr") := tstrsplit(RegionData$SNP, ":")[[1]] ]
      setkey(RegionData, "chr")
      RegionData = RegionData[chr == as.numeric(chrom)]
      RegionData[,chr:=NULL]
    }

if(nrow(RegionData) != 0){
if(!is.null(markerInfo)){
  setkey(RegionData, "SNP")
  RegionData = merge(RegionData, markerInfo, by.x = "SNP", by.y = "ID", all.x = T, sort = F) 

  if(!is.null(markerInfo$ID2)){
        RegionData = merge(RegionData, markerInfo, by.x = "SNP", by.y = "ID2", all.x = T, sort = F)
	#SNP REGION     ANNO CHROM.x POS.x genoIndex2.x  ID2 genoIndex_prev.x CHROM.y POS.y    ID genoIndex2.y genoIndex_prev.y
        setnames(RegionData, "genoIndex.x", "genoIndex")
        setnames(RegionData, "genoIndex.y", "genoIndex2")
        #setnames(RegionData, "genoIndex.y", "genoIndex")
        setnames(RegionData, "CHROM.x", "CHROM")
        setnames(RegionData, "CHROM.y", "CHROM2")
        setnames(RegionData, "POS.x", "POS")
        setnames(RegionData, "POS.y", "POS2")

	if(!is.null(RegionData$genoIndex_prev.y)){
		setnames(RegionData, "genoIndex_prev.x", "genoIndex_prev")
		setnames(RegionData, "genoIndex_prev.y", "genoIndex_prev2")
		#markerInfo[,genoIndex_prev.y:=NULL]

	}
	posNA = which(is.na(RegionData$genoIndex) & !is.na(RegionData$genoIndex2))

	if(length(posNA) != 0){
		RegionData$genoIndex[posNA] = RegionData$genoIndex2[posNA]
		RegionData$CHROM[posNA] = RegionData$CHROM2[posNA]
		RegionData$POS[posNA] = RegionData$POS2[posNA]
		if(!is.null(RegionData$genoIndex_prev)){
			RegionData$genoIndex_prev[posNA] = RegionData$genoIndex_prev2[posNA]			
		}
	}
	#RegionData$genoIndex2 = NULL	
  }
  posNA = which(is.na(RegionData$genoIndex))
  
    if(length(posNA) != 0){
      RegionData = RegionData[which(!is.na(RegionData$genoIndex))]
      cat(length(posNA)," markers in 'RegionFile' are not in 'GenoFile'.\n")
    }
   setorderv(RegionData, col=c("CHROM", "POS"))
}

}


if(nrow(RegionData) != 0){

  
  #HeaderInRegionData = colnames(RegionData)
  HeaderInRegionData = unique(RegionData$ANNO)
  RegionAnnoHeaderList = list()
  if(length(annoVec) == 0){
        stop("At least one annotation is required\n")
  }
  for(q in 1:length(annoVec)){
        RegionAnnoHeaderList[[q]] = strsplit(annoVec[q],";")[[1]]
  }

  RegionList = list()
  uRegion = unique(RegionData$REGION)
  #RegionData = as.data.frame(RegionData)


  for(r in uRegion0){
             #print(paste0("Analyzing region ",r,"...."))
    #print(RegionData$REGION)
    #print(r)
    #which(as.vector(RegionData$REGION) == r)
   
    posSNP = which(RegionData$REGION == r)
    if(length(posSNP) > 0){
    SNP = RegionData$SNP[posSNP]

    if(nline_per_gene == 3){	
      WEIGHT = RegionData$WEIGHT[posSNP]
    }

    if(any(duplicated(SNP)))
      stop("Please check RegionFile: in region ", r,": duplicated SNPs exist.")

    if(!is.null(markerInfo)){
      genoIndex = as.numeric(RegionData$genoIndex[posSNP])
      if(!is.null(RegionData$genoIndex_prev)){
		genoIndex_prev = as.numeric(RegionData$genoIndex_prev[posSNP])	
      }
      chrom = RegionData$CHROM[posSNP]
      #uchrom = unique(chrom)
      #if(length(uchrom) != 1)
      #  stop("In region ",r,", markers are from multiple chromosomes.")
    }

    annoIndicatorMat = matrix(0, nrow=length(posSNP), ncol=length(annoVec))
    annoVecNew = c()
        for(q in 1:length(annoVec)){
        indiceVec = which(RegionData$ANNO[posSNP] %in% RegionAnnoHeaderList[[q]])
        if(length(indiceVec) > 0){
                annoVecNew = c(annoVecNew, annoVec[q])
                annoIndicatorMat[indiceVec, q] = 1
                #annoIndicatorMat[which(RegionData$ANNO[posSNP] %in% RegionAnnoHeaderList[[q]]), q] = 1
        }
    }

  RegionAnnoHeaderListNew = list()
  if(length(annoVecNew) == 0){
        warning("No markers are found for at least one annotation, so region ", r, " is skipped\n")
        #stop("At least one annotation is required\n")
  }else{
   if(length(annoVecNew) < length(annoVec)){
        annoIndicatorMat = matrix(0, nrow=length(posSNP), ncol=length(annoVecNew))
    for(q in 1:length(annoVecNew)){
        RegionAnnoHeaderListNew[[q]] = strsplit(annoVecNew[q],";")[[1]]
        indiceVec = which(RegionData$ANNO[posSNP] %in% RegionAnnoHeaderListNew[[q]])
               #if(length(indiceVec) > 0){
        annoIndicatorMat[indiceVec, q] = 1
                #annoIndicatorMat[which(RegionData$ANNO[posSNP] %in% RegionAnnoHeaderList[[q]]), q] = 1
    }
   }else{
        annoVecNew = annoVec
   }


  }

  annoIndicatorMat_rmind = which(rowSums(annoIndicatorMat) == 0)
  if(length(annoIndicatorMat_rmind) > 0){
    SNP = SNP[-annoIndicatorMat_rmind]

    if(nline_per_gene == 3){
      WEIGHT = WEIGHT[-annoIndicatorMat_rmind]
    }

    annoIndicatorMat = annoIndicatorMat[-annoIndicatorMat_rmind,,drop=F]
    if(!is.null(markerInfo)){
      genoIndex = genoIndex[-annoIndicatorMat_rmind]
     if(!is.null(RegionData$genoIndex_prev)){
                genoIndex_prev = genoIndex_prev[-annoIndicatorMat_rmind]
      }

    }
  }

  if(nline_per_gene != 3){
	WEIGHT = c(0)
  }	  

   if(!is.null(markerInfo)){
    RegionList[[r]] = list(SNP = SNP,
			   WEIGHT=WEIGHT,
                           annoIndicatorMat = annoIndicatorMat,
                           genoIndex =  as.character(format(genoIndex, scientific = FALSE)),
#                           chrom = uchrom,
                           annoVec = annoVecNew)
    if(!is.null(RegionData$genoIndex_prev)){
	RegionList[[r]]$genoIndex_prev = as.character(format(genoIndex_prev, scientific = FALSE)) 
     }else{
	RegionList[[r]]$genoIndex_prev = c("-1")	
     }

   }else{
    RegionList[[r]] = list(SNP = SNP,
			   WEIGHT=WEIGHT,
                           annoIndicatorMat = annoIndicatorMat,
 #                          chrom = uchrom,
                           annoVec = annoVecNew)
    }
  
   }else{ #if(length(posSNP) > 0){
     RegionList[[r]] = list(SNP = NULL)

   }
}

}else{#if(nrow(RegionData) == 0){
	RegionList = NULL
}


  return(RegionList)
}





mainRegionURV = function(NullModelClass = "SAIGE_NULL_Model",
                         genoType,
                         genoIndex,
                         n)
{
  if(NullModelClass == "SAIGE_NULL_Model")
    obj.mainRegionURV = mainRegionURVInCPP("SAIGE", genoType, genoIndex, n)

  return(obj.mainRegionURV)
}

