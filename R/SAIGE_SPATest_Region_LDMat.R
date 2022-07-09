generate_LDMat_forMataRegion = function(bgenFile = "",
                 bgenFileIndex = "",
                 sampleFile = "",
                 vcfFile = "",
                 vcfFileIndex = "",
                 vcfField = "DS",
                 savFile = "",
                 savFileIndex = "",
                 bedFile="",
                 bimFile="",
                 famFile="",
                 AlleleOrder = "alt-first", #new
                 idstoIncludeFile = "",
                 rangestoIncludeFile = "",
                 chrom = "", #for vcf file
                 max_missing = 0.15,  #new
                 impute_method = "best_guess",  #"mean", "minor", "best_guess"     #new
                 min_MAC = 0.5,
                 min_MAF = 0,
                 min_Info = 0,
                 is_imputed_data = FALSE, #new
                 SAIGEOutputFile = "",
                 groups_per_chunk = 100,
                 markers_per_chunk_in_groupTest = 100, #new
                 groupFile="",
                 dosage_zerod_cutoff = 0.2,
                 dosage_zerod_MAC_cutoff = 10,
                 annotation_in_groupTest ="missense;lof;synonymous",  #new
                 maxMAF_in_groupTest = 0.5,
                 is_overwrite_output = TRUE, 
		 sampleFile_include_inLDMat = "")
{
    OutputFile = SAIGEOutputFile
    OutputFileIndex=NULL
    if(is.null(OutputFileIndex)){OutputFileIndex = paste0(OutputFile, ".index")}

    maxMAFLimit = max(maxMAF_in_groupTest)

    ##check the variance ratio file and extract the variance ratio vector
    setGlobalVarsInCPP_LDmat(impute_method,
			dosage_zerod_cutoff,
                        dosage_zerod_MAC_cutoff,
                        max_missing,
			maxMAFLimit,
                        min_MAF,
                        min_MAC,
                        min_Info,
			markers_per_chunk_in_groupTest,
			SAIGEOutputFile)

   if(sampleFile_include_inLDMat != ""){
     Check_File_Exist(sampleFile_include_inLDMat, "sampleFile_include_inLDMat")
     sampleInModel = data.table:::fread(sampleFile_include_inLDMat, header = F, colClasses = list(character = 1), select=c(1), data.table=F)[,1]
     cat(length(sampleInModel), " samples were found in sampleFile_include_inLDMat: ", sampleFile_include_inLDMat, "\n")

   }else{
     sampleInModel = NULL
     cat("sampleFile_include_inLDMat is not specified, so all samples in the genotype/dosage file are included for LD matrices\n")
   }	   


   objGeno = setGenoInput_LDmat(bgenFile = bgenFile,
                 bgenFileIndex = bgenFileIndex,
                 vcfFile = vcfFile,   #not activate yet
                 vcfFileIndex = vcfFileIndex,
                 vcfField = vcfField,
                 savFile = savFile,
                 savFileIndex = savFileIndex,
                 sampleFile = sampleFile,
                 bedFile=bedFile,
                 bimFile=bimFile,
                 famFile=famFile,
                 idstoIncludeFile = idstoIncludeFile,
                 rangestoIncludeFile = rangestoIncludeFile,
                 chrom = chrom,
                 AlleleOrder = AlleleOrder,
		 sampleInModel = sampleInModel)

	SAIGE.Region.LDmat(OutputFile,
			   groupFile,
			   annotation_in_groupTest,
			   markers_per_chunk_in_groupTest,
			   objGeno$genoType,
                           objGeno$markerInfo,
		           is_imputed_data,
		           groups_per_chunk,
		           chrom,
		           is_overwrite_output)

        closeOutfile_single_LDmat();

}


#	SMat.list: list object of G^TG matrices from each cohorts

#	SMat_Info.list: list object of dataframe with 

# SNPID, 	MajorAllele, MinorAllele, MAC (or MAF), N
SAIGE.Region.LDmat = function(
                        OutputFile,
                        groupFile,
                        annolist,
                        markers_per_chunk_in_groupTest,
                        genoType,
                        markerInfo,
                        isImputation,
                        groups_per_chunk,
                        chrom,
			isOverWriteOutput)
{
  OutputFileIndex = NULL
  if(is.null(OutputFileIndex)){
    OutputFileIndex = paste0(OutputFile, ".index")
  }
  outList = checkOutputFile(paste0(OutputFile,".marker_info.txt"), OutputFileIndex, "Region", 1, isOverWriteOutput) # Check 'Util.R'

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

 cat("isappend ", isappend, "\n")
 isOpenOutFile = openOutfile_single_LDmat(isappend)
 if(!isOpenOutFile){
   stop(OutputFile, ".marker_info.txt can't be opened to output marker information.\n") 
 } 	 

  #n = length(mu) #sample size
  ##check group file
  region_list = checkGroupFile(groupFile)
  nRegions = region_list$nRegions
  is_weight_included = region_list$is_weight_included

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


    chrom1 = "FakeCHR";

  gc()
  num_region = 0
  mth = 0

  numberRegionsInChunk = 0
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
      
      RegionList = SAIGE.getRegionList_forLDmat(marker_group_line, nline_per_gene, annolist, markerInfo, chrom)  ##
      #print(RegionList)
      cat("Read in ", nregions_ro_read, " region(s) from the group file.\n")
      mth = 0
      #numberRegionsInChunk = length(RegionList)
      numberRegionsInChunk = nregions_ro_read
    }

   mth = mth + 1




   if(!is.null(RegionList)){
    pval.Region = NULL
    region = RegionList[[mth]]

    annolistsub = region$annoVec ##annolist contains all markers
    regionName = names(RegionList)[mth]
    i = i + 1
    if(!is.null(region$SNP) & length(annolistsub) > 0){

      SNP = region$SNP
      if(genoType == "vcf"){
        SNPlist = paste(c(regionName, SNP), collapse = "\t")
        if(length(SNP) == 1){
                 fakem = strsplit(SNP, split=":")[[1]]
                 fakemb = paste(c(fakem[1], as.numeric(fakem[2
])+2, "N", "N"), collapse=":")
                 SNPlisttemp = paste(c(SNPlist, fakemb), collapse = "\t")
                set_iterator_inVcf(SNPlisttemp, chrom1, 1, 250000000)
        }else{
                set_iterator_inVcf(SNPlist, chrom1, 1, 250000000)
        }

         isVcfEnd =  check_Vcf_end()



        if(!isVcfEnd){
                region$genoIndex = rep("0", length(SNP))
                region$genoIndex_prev = rep("0", length(SNP))
        }else{
                warning("No markers in region ", regionName, " are found in the VCF file")
                next
        }
      }

      annoIndicatorMat = region$annoIndicatorMat

      #chrom = region$chrom
      print(paste0("Analyzing Region ", regionName, " (",i-1,"/",nRegions,")."))

      #LDmatInCPP(genoType, region$genoIndex_prev, region$genoIndex, annoIndicatorMat, maxMAFlist, OutputFile, traitType, n, P1Mat, P2Mat, regionTestType, isImputation, WEIGHT, weight_cond, is_single_in_groupTest, is_output_markerList_in_groupTest, annolistsub, regionName, is_fastTest, is_output_moreDetails)
      n=Unified_getSampleSizeinAnalysis(genoType)
      LDmatRegionInCPP(genoType, region$genoIndex_prev, region$genoIndex, annoIndicatorMat, OutputFile, n, isImputation, annolistsub, regionName)

    }else{#if(!is.null(region$SNP) & length(annolistsub) > 0){
      cat(regionName, " is empty.\n")	
   }	   

  }else{#if(!is.null(RegionList)){
     cat("The chunk is empty\n")
     #mth = 0
     mth = numberRegionsInChunk
     i = i + numberRegionsInChunk
     pval.Region = NULL
   }    
  indexChunk = i
  #Start = (i==1)
  Start = (cth_chunk_to_output==1)
  End = (i==nRegions)
  nEachChunk = 1

  if(mth ==  numberRegionsInChunk){
  message1 = "This is the output index file for SAIGE package to record the end point in case users want to restart the analysis. Please do not modify this file."
  message2 = paste("This is to generate LD matrices for rare-variant meta-analyses.")
  message3 = paste("nEachChunk =", nEachChunk)
  message4 = paste("Have completed the analysis of chunk", indexChunk-1)
  message5 = "Have completed the analyses of all chunks."
  cat("write to output\n")
  if(Start){
    write.table(c(message1, message2, message3), OutputFileIndex,
                quote = F, sep = "\t", append = F, col.names = F, row.names = F)
  }
    write.table(message4, OutputFileIndex, quote = F, sep = "\t", append = T, col.names = F, row.names = F)

  if(End){
    write.table(message5, OutputFileIndex, quote = F, sep = "\t", append = T, col.names = F, row.names = F)
  }

}  
  }
gc()
close(gf)
}



SAIGE.getRegionList_forLDmat = function(marker_group_line,
                         nline_per_gene,
                         annoVec, #c("lof","lof;missense"
                         markerInfo,
                         chrom="")
{

  chrom_nochr = gsub("CHR", "", chrom, ignore.case = T)
  # read group file
  ngroup<-length(marker_group_line)/nline_per_gene
  #cat("ngroup is ", ngroup, "\n")
  RegionData = NULL
  geneList = c()
  for(i in 1:ngroup){
          marker_group_line_list = strsplit(marker_group_line[1+(i-1)*nline_per_gene], split="[\ \t]+")[[1]]
          gene=marker_group_line_list[1]
          var=marker_group_line_list[3:length(marker_group_line_list)]
          marker_group_line_list_anno = strsplit(marker_group_line[2+(i-1)*nline_per_gene], split="[\ \t]+")[[1]]
          anno=marker_group_line_list_anno[3:length(marker_group_line_list_anno)]
          if(nline_per_gene == 3){
              marker_group_line_list_weight = strsplit(marker_group_line[3+(i-1)*nline_per_gene], split="[\ \t]+")[[1]]
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
      RegionData = RegionData[chr == chrom | chr == chrom_nochr]
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
        #stop("At least one annotation is required\n")
        cat("annotation_in_groupTest is ALL so all markers will be included\n")
        annoVec = c("ALL")
        RegionAnnoHeaderList[[1]] = "ALL" 
  }else{  
  	for(q in 1:length(annoVec)){
        	RegionAnnoHeaderList[[q]] = strsplit(annoVec[q],";")[[1]]
		cat("Markers with annotations in ", annoVec[q], " will be included\n")
  	}
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


    if(length(annoVec) == 1 & annoVec[1] == "ALL"){
		annoIndicatorMat[,1] = 1
		annoVecNew = annoVec
	}else{	
          for(q in 1:length(annoVec)){
            indiceVec = which(RegionData$ANNO[posSNP] %in% RegionAnnoHeaderList[[q]])
            if(length(indiceVec) > 0){
                annoVecNew = c(annoVecNew, annoVec[q])
                annoIndicatorMat[indiceVec, q] = 1
                #annoIndicatorMat[which(RegionData$ANNO[posSNP] %in% RegionAnnoHeaderList[[q]]), q] = 1
            }
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
    ##VCF
    if(length(SNP) > 0){
        orderposind = order(as.numeric(tstrsplit(SNP, ":")[[2]]))
        SNP = SNP[orderposind]
        annoIndicatorMat = annoIndicatorMat[orderposind,, drop=F]
        if(length(WEIGHT) == length(orderposind)){
            WEIGHT = WEIGHT[orderposind]
        }
    }
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





