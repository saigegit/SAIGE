#' Check if sample identifiers are stored in a BGEN file (for BGEN v1.2)
#'
#' Check if sample identifiers are stored in a BGEN file (for BGEN v1.2)
#'
#' @param bgenFile a character of BGEN file. Sometimes, BGEN file does not include sample IDs. This information can be extracted from BGEN file. Please refer to https://www.well.ox.ac.uk/~gav/bgen_format/spec/v1.2.html for more details.
#' @examples
#'
#' BGENFile = system.file("extdata", "example_bgen_1.2_8bits.bgen", package = "SAIGE")
#' checkIfSampleIDsExist(BGENFile)
#' @export
checkIfSampleIDsExist = function(bgenFile)
{
  con = file(bgenFile, "rb")
  seek(con, 4)
  LH = readBin(con, n = 1, what = "integer", size = 4)
  seek(con, 4 + LH - 4)
  header = rawToBits(readBin(con, n = 4, what = "raw", size = 1, signed = FALSE))
  close(con)
  return(header[32] == 01)
}

checkGenoInput = function(bgenFile = "",
                 bgenFileIndex = "",
                 vcfFile = "",
		 vcfFileIndex = "",
                 vcfField = "DS",
                 savFile = "",
		 savFileIndex = "",
                 sampleFile = "",
                 bedFile="",
                 bimFile="",
                 famFile="", 
		 sampleInModel = NULL){
   
    if(is.null(sampleInModel)){
    	stop("sampleInModel is not specified.")
    }	    
    # genotype data
    if (bgenFile != "") {
        Check_File_Exist(bgenFile, "bgenFile")
        Check_File_Exist(bgenFileIndex, "bgenFileIndex")
        dosageFileType = "bgen"
    }else if (vcfFile != "") {
        Check_File_Exist(vcfFile, "vcfFile")

	vcfFileIndex_require = paste(vcfFile, ".csi", sep = "")
        if(vcfFileIndex != ""){
                if(vcfFileIndex != vcfFileIndex_require){
                        stop("ERROR! vcfFileIndex is specfied. Please note it must be ", vcfFileIndex_require, "\n")
                }else{
			cat("Please note the argument vcfFileIndex will not be used in future versions because the vcf index file must has the name 'vcfFile'.csi\n")
		}	
        }else{
                cat("vcfFileIndex is not specified. ", vcfFileIndex_require, " will be used\n")
        }
    	vcfFileIndex = vcfFileIndex_require

        if (!file.exists(paste(vcfFile, ".csi", sep = ""))) {
            stop("ERROR! vcfFileIndex ", vcfFileIndex, " does not exist\n")
        }
        dosageFileType = "vcf"

    }else if (savFile != "") {
        Check_File_Exist(savFile, "savFile")
	#Check_File_Exist(savFileIndex, "savFileIndex")
	savFileIndex_require = paste(savFile, ".s1r", sep = "")
	if(savFileIndex != ""){
		if(savFileIndex != paste(savFile, ".s1r", sep = "")){
			stop("ERROR! savFileIndex is specfied. Please note it must be ", savFileIndex_require, "\n") 
		}else{
			cat("Please note the argument savFileIndex will not be used in future versions because the sav index file must has the name 'savFile'.s1r\n")
		}	
	}else{
		cat("savFileIndex is not specified. ", savFileIndex_require, " will be used\n")
	}	
	savFileIndex = savFileIndex_require
	if (!file.exists(savFileIndex)) {
            stop("ERROR! savFileIndex ", savFileIndex, " does not exist\n")
        }
	vcfFile = savFile
	vcfFileIndex = savFileIndex
        dosageFileType = "vcf"

    }else if(bedFile != ""){
	Check_File_Exist(bedFile, "bedFile")
	if(bimFile == ""){
		bimFile = gsub("bed$", "bim", bedFile)
	}
	if(famFile == ""){
		famFile = gsub("bed$", "fam", bedFile)
	}	

	Check_File_Exist(bimFile, "bimFile")
	Check_File_Exist(famFile, "famFile")
	dosageFileType = "plink"
    }	    
    
    cat("dosageFile type is ", dosageFileType, "\n")

    if(dosageFileType == "vcf"){
        #if (chrom == "") {
        #    stop("ERROR! chrom needs to be specified for the vcf file\n")
        #}
	if(vcfField != "DS" & vcfField != "GT"){
	    stop("vcfField has to be DS or GT\n")	
	}
    }

    return(dosageFileType)
}	


splitreformatMarkerIDinBgen = function(x){
	a = strsplit(x, split=":")[[1]]
	b=x
	if(length(a) == 4){
		b = paste0(a[1],":",a[2],"_",a[3],"/",a[4])
	}
	return(b)
}

setGenoInput = function(bgenFile = "",
                 bgenFileIndex = "",
                 vcfFile = "",
                 vcfFileIndex = "",
                 vcfField = "DS",
                 savFile = "",
                 savFileIndex = "",
                 sampleFile = "",
                 bedFile="",
                 bimFile="",
                 famFile="",
                 idstoIncludeFile = "",
                 rangestoIncludeFile = "",
                 chrom = "",
		 AlleleOrder = NULL,
		 sampleInModel = NULL)
{

  dosageFileType = checkGenoInput(bgenFile = bgenFile,
                 bgenFileIndex = bgenFileIndex,
                 vcfFile = vcfFile,
                 vcfFileIndex = vcfFileIndex,
                 vcfField = vcfField,
                 savFile = savFile,
                 savFileIndex = savFileIndex,
                 bedFile = bedFile,
                 bimFile = bimFile,
                 famFile = famFile,
		 sampleInModel = sampleInModel)

  ########## ----------  Plink format ---------- ##########
  
  if(dosageFileType == "plink"){
    if(is.null(AlleleOrder)) AlleleOrder = "alt-first"

    cat("allele order in the plink file is ", AlleleOrder, ".\n")

    if(bimFile == ""){
    	bimFile = gsub("bed$", "bim", bedFile)
    }
    if(famFile == ""){
    	famFile = gsub("bed$", "fam", bedFile)
    } 
    markerInfo = data.table::fread(bimFile, header = F, select = c(1, 2, 4, 5, 6))
    #markerInfo = as.data.frame(markerInfo)
    #markerInfo is a data.table 
    if(AlleleOrder == "alt-first")
      names(markerInfo) = c("CHROM", "ID", "POS", "ALT", "REF")
      #markerInfo = markerInfo[,c(1,2,3,5,4)]  # https://www.cog-genomics.org/plink/2.0/formats#bim
    if(AlleleOrder == "ref-first")
      names(markerInfo) = c("CHROM", "ID", "POS", "REF", "ALT")
      #markerInfo = markerInfo[,c(1,2,3,4,5)]  # https://www.cog-genomics.org/plink/2.0/formats#bim
    #colnames(markerInfo) = c("CHROM", "ID", "POS", "REF", "ALT")
    #colnames(markerInfo) = c("CHROM", "POS", "ID", "REF", "ALT")
    #markerInfo$ID = paste0(markerInfo$CHROM,":",markerInfo$POS, "_", markerInfo$REF,"/", markerInfo$ALT)
    markerInfo$genoIndex = 1:nrow(markerInfo) - 1  # -1 is to convert 'R' to 'C++' 
    markerInfo$genoIndex_prev = rep(0, length(markerInfo$genoIndex))
    if(chrom != ""){
      markerInfo = markerInfo[which(markerInfo[,1] == chrom), ]
    }
    #markerInfo$ID2 = lapply(markerInfo$ID, splitreformatMarkerIDinBgen)    
    markerInfo$ID2 = paste0(markerInfo$CHROM,":",markerInfo$POS, ":", markerInfo$REF,":", markerInfo$ALT)
#    markerInfo[,POS:=NULL]
    #print(is.data.table(markerInfo))
    markerInfo[,REF:=NULL]
    markerInfo[,ALT:=NULL]
    setkeyv(markerInfo, c("ID","ID2"))
    #markerInfo$genoIndex_prev = NULL
    #sampleInfo = data.table::fread(famFile, select = c(2), data.table=F)
    #samplesInGeno = sampleInfo[,1]
    #SampleIDs = updateSampleIDs(SampleIDs, samplesInGeno)
    #markerInfo$ID = paste0(markerInfo$CHROM,":", markerInfo$POS ,"_", markerInfo$REF, "/", markerInfo$ALT) 
    setPLINKobjInCPP(bimFile, famFile, bedFile, sampleInModel, AlleleOrder)
  }
  
  ########## ----------  BGEN format ---------- ##########
  
  if(dosageFileType == "bgen"){
    
    
    if(is.null(AlleleOrder)) AlleleOrder = "ref-first"
    
    if(sampleFile != "" | !checkIfSampleIDsExist(bgenFile)){
	print("Sample IDs were not found in the bgen file.")
	Check_File_Exist(sampleFile)
	sf = file(sampleFile, "r")
	first_sample_line = readLines(sf, n = 1)
	close(sf)
	first_sample_line_list = strsplit(first_sample_line, split="[\ \t]+")[[1]]
        		
	if(first_sample_line == "ID_1 ID_2 missing sex" | first_sample_line == "ID_1 ID_1 0 0" | length(first_sample_line_list) == 4){
	    cat("sample file is in the bgenix format\n")
	    sampleData = data.table::fread(sampleFile, header=F, colClasses = rep("character", 4), data.table=F, skip=2)
	    samplesInGeno = as.character(sampleData[,2])
	}else{
	    if(length(first_sample_line_list) == 1){
	        cat("sample file only has one column and has no header\n")
	        sampleData = data.table::fread(sampleFile, header=F, colClasses = c("character"), data.table=F)
                samplesInGeno = as.character(sampleData[,1])
	    }else{
		stop("sample file should either has one column and no header or in the bgenix format with four columns\n")	
	    }
	}
    }else{
	samplesInGeno = getSampleIDsFromBGEN(bgenFile)
 	#print(samplesInGeno[1:100])		    
    } 


    db_con <- RSQLite::dbConnect(RSQLite::SQLite(), bgenFileIndex)
    on.exit(RSQLite::dbDisconnect(db_con), add = TRUE)
    markerInfo = dplyr::tbl(db_con, "Variant")
    #markerInfo = as.data.frame(markerInfo)
    #markerInfo = as.data.frame(markerInfo)
    markerInfo = as.data.table(markerInfo, keep.rownames = FALSE)
    rmcol = setdiff(names(markerInfo), c("chromosome", "position", "rsid", "allele1", "allele2", 'file_start_position', 'size_in_bytes'))
    markerInfo = markerInfo[, !..rmcol]
    if(AlleleOrder == "alt-first")
      names(markerInfo) = c("CHROM", "POS", "ID", "ALT", "REF", "genoIndex", "size_in_bytes")
      #markerInfo = markerInfo[,c(1,2,3,6,5,7)]  # https://www.well.ox.ac.uk/~gav/bgen_format/spec/v1.2.html
    if(AlleleOrder == "ref-first")
      names(markerInfo) = c("CHROM", "POS", "ID", "REF", "ALT", "genoIndex", "size_in_bytes")

      #markerInfo = markerInfo[,c(1,2,3,5,6,7)]  # https://www.well.ox.ac.uk/~gav/bgen_format/spec/v1.2.html
    #markerInfo$ID2 = lapply(markerInfo$ID, splitreformatMarkerIDinBgen)
    markerInfo$ID2 = paste0(markerInfo$CHROM,":", markerInfo$POS ,":", markerInfo$REF, ":", markerInfo$ALT)
    markerInfo$genoIndex_prev = markerInfo$genoIndex + markerInfo$size_in_bytes
#    markerInfo[,POS:=NULL]
    markerInfo[,REF:=NULL]
    markerInfo[,ALT:=NULL]
    markerInfo[,size_in_bytes:=NULL]
    
    if(chrom != ""){
      markerInfo = markerInfo[which(markerInfo[,1] == chrom), ]
    }  
    	
    setkeyv(markerInfo, c("ID","ID2"))
    #markerInfo$ID2 = paste0(markerInfo$CHROM,":", markerInfo$POS ,"_", markerInfo$ALT, "/", markerInfo$REF)
    setBGENobjInCPP(bgenFile, bgenFileIndex, t_SampleInBgen = samplesInGeno, t_SampleInModel = sampleInModel, AlleleOrder)
  }
  

  if(dosageFileType == "vcf"){
    if(idstoIncludeFile != "" & rangestoIncludeFile != ""){
      stop("We currently do not support both 'idstoIncludeFile' and 'rangestoIncludeFile' at the same time for vcf files\n")
    }
    #if(chrom==""){
    #  stop("chrom needs to be specified for VCF/BCF/SAV input\n")
    #}
    #markerInfo = NULL
  }


  #Files = c("idstoIncludeFile", "idstoExcludeFile", "rangestoIncludeFile", "rangestoExcludeFile")

  anyInclude = FALSE
  #anyExclude = FALSE

  markersInclude = c()
  #markersExclude = c()
  IDsToInclude = NULL
  RangesToInclude = NULL
  if(idstoIncludeFile != ""){
    IDsToInclude = data.table::fread(idstoIncludeFile, header = F, colClasses = c("character"), data.table=F)
    if(ncol(IDsToInclude) != 1){
      stop("'idstoIncludeFile' of ", idstoIncludeFile, " should only include one column.")
    }else{
      IDsToInclude = IDsToInclude[,1]
      cat(length(IDsToInclude), " marker IDs are found in the idstoIncludeFile file\n")
    }

  if(dosageFileType != "vcf"){
    if(!is.null(markerInfo$ID2)){

    	posRows = which((markerInfo$ID %in% IDsToInclude) | (markerInfo$ID2 %in% IDsToInclude))
    }else{
	posRows = which((markerInfo$ID %in% IDsToInclude))

    }	    
    if(length(posRows) != 0){
      cat(length(posRows), " markers in idstoIncludeFile are found in the geno/dosage file.\n") 
      markersInclude = c(markersInclude, markerInfo$ID[posRows])
    }  
 
  }

    anyInclude = TRUE
  }

  if(rangestoIncludeFile != ""){
    RangesToInclude = data.table::fread(rangestoIncludeFile, header = F, colClasses = c("character", "numeric", "numeric"), data.table=F)
    if(ncol(RangesToInclude) != 3){
      stop("rangestoIncludeFile should only include three columns.")
    }else{
      cat(nrow(RangesToInclude), " ranges are in rangestoIncludeFile\n")
    }

    colnames(RangesToInclude) = c("CHROM", "START", "END")
 if(dosageFileType != "vcf"){ 
   if(nrow(RangesToInclude) > 0){
    for(i in 1:nrow(RangesToInclude)){
      CHROM1 = RangesToInclude$CHROM[i]
      START = RangesToInclude$START[i]
      END = RangesToInclude$END[i]
      posRows = with(markerInfo, which(CHROM == CHROM1 & POS >= START & POS <= END))
      if(length(posRows) != 0)
        markersInclude = c(markersInclude, markerInfo$ID[posRows])
	cat(length(markerInfo$ID[posRows]), " markers in the range ", i, " are found in the geno/dosage file.\n")
    }
   }
  } 
    anyInclude = TRUE
  }

if(FALSE){

  if(!is.null(idstoExcludeFile)){
    if(anyInclude)
      stop("We currently do not support both 'IncludeFile' and 'ExcludeFile' at the same time.")
    IDsToExclude = data.table::fread(idstoExcludeFile, header = F, colClasses = c("character"))
    if(ncol(IDsToExclude) != 1)
      stop("idstoExcludeFile should only include one column.")
    IDsToExclude = IDsToExclude[,1]

    posRows = which(markerInfo$ID %in% IDsToExclude)
    if(length(posRows) != 0)
      markersExclude = c(markersExclude, markerInfo$ID[posRows])
    anyExclude = TRUE
  }

  if(!is.null(rangestoExcludeFile)){
    if(anyInclude)
      stop("We currently do not support both 'IncludeFile' and 'ExcludeFile' at the same time.")
    RangesToExclude = data.table::fread(rangestoExcludeFile, header = F, colClasses = c("character", "numeric", "numeric"))
    if(ncol(RangesToExclude) != 3)
      stop("rangestoExcludeFile should only include three columns.")

    colnames(RangesToExclude) = c("CHROM", "START", "END")

    for(i in 1:nrow(RangesToExclude)){
      CHROM1 = RangesToExclude$CHROM[i]
      START = RangesToExclude$START[i]
      END = RangesToExclude$END[i]
      posRows = with(markerInfo, which(CHROM == CHROM1 & POS >= START & POS <= END))
      if(length(posRows) != 0)
        markersExclude = c(markersExclude, markerInfo$ID[posRows])
    }
    anyExclude = TRUE
  }
}

  markersInclude = unique(markersInclude)
#  markersExclude = unique(markersExclude)

  # return genotype
  #cat("Based on the 'GenoFile' and 'GenoFileIndex',", genoType, "format is used for genotype data.\n")

  if(anyInclude){
   if(dosageFileType != "vcf"){	  
    markerInfo = subset(markerInfo, ID %in% markersInclude)
   }
  }
#  if(anyExclude)
#    markerInfo = subset(markerInfo, !ID %in% markersExclude)

#  anyQueue = anyInclude | anyExclude

  if(dosageFileType == "vcf"){

    if(vcfFile != ""){
	vcfFileIndex = paste(vcfFile, ".csi", sep = "") 
    }else if(savFile != ""){	    
	vcfFile = savFile
        vcfFileIndex = paste(vcfFile, ".s1r", sep = "")
    }
	    
    setVCFobjInCPP(vcfFile, vcfFileIndex, vcfField, t_SampleInModel = sampleInModel)
    if(!is.null(IDsToInclude)){
      SNPlist = paste(c("set1", IDsToInclude), collapse = "\t")
      in_chrom="fake_chrom"
      in_beg_pd=1
      in_end_pd=250000000
      set_iterator_inVcf(SNPlist, in_chrom, in_beg_pd, in_end_pd)
      anyInclude = TRUE
    }

    if(!is.null(RangesToInclude)){
      CHROM1 = RangesToInclude$CHROM
      START = RangesToInclude$START
      END = RangesToInclude$END
      if(length(CHROM1) > 1){
        stop("We do not support query with multiple regions for vcf file. Please only include one region in the ", rangestoIncludeFile, "\n")
      }else{
	inSNPlist=""
        in_chrom=CHROM1[1]
  	in_beg_pd=START[1]
	in_end_pd=END[1]
	cat("Tests will be performed for chromosome: ",in_chrom, ", start: ", in_beg_pd," end: ", in_end_pd,"\n")	
        set_iterator_inVcf(inSNPlist, in_chrom, in_beg_pd, in_end_pd)
	anyInclude = TRUE
      }
    } 
    #markerInfo = data.table(CHROM=NULL, POS=NULL, genoIndex_prev=NULL, genoIndex=NULL)
    markerInfo = NULL
  }
  #genoList = list(genoType = genoType, markerInfo = markerInfo, SampleIDs = SampleIDs, AlleleOrder = AlleleOrder, GenoFile = GenoFile, GenoFileIndex = GenoFileIndex, anyQueue = anyQueue)
  #genoList = list(dosageFileType = dosageFileType, markerInfo = markerInfo, anyQueue = anyQueue, genoType = dosageFileType)
  #genoList = list(dosageFileType = dosageFileType, markerInfo = markerInfo, genoType = dosageFileType)
 # if(nrow(markerInfo) != 0){
 #   markerInfo[,POS:=NULL]
 # }

  genoList = list(markerInfo = markerInfo, genoType = dosageFileType, anyInclude = anyInclude)
  return(genoList)
}




updateSampleIDs = function(SampleIDs, samplesInGeno)
{
  if(is.null(SampleIDs)){
    print("Since 'SampleIDs' not specified, we use all samples in 'GenoFile'.")
    SampleIDs = samplesInGeno;
  }
  
  if(any(!SampleIDs %in% samplesInGeno))
    stop("At least one sample from 'SampleIDs' are not in 'GenoFile' and 'GenoFileIndex'.")
  
  return(SampleIDs)
}

# https://www.well.ox.ac.uk/~gav/bgen_format/spec/v1.2.html
#' Extract sample identifiers from BGEN file (for BGEN v1.2)
#' 
#' Extract sample identifiers from BGEN file (for BGEN v1.2)
#' 
#' @param bgenFile a character of BGEN file. 
#' @examples
#' 
#' BGENFile = system.file("extdata", "example_bgen_1.2_8bits.bgen", package = "SAIGE")
#' getSampleIDsFromBGEN(BGENFile)
#' @export
getSampleIDsFromBGEN = function(bgenFile)
{
  if(!checkIfSampleIDsExist(bgenFile))
    stop("The BGEN file does not include sample identifiers. Please refer to ?checkIfSampleIDsExist and ?SAIGE.ReadGeno for more details")
  con = file(bgenFile, "rb")
  seek(con, 4)
  LH = readBin(con, n = 1, what = "integer", size = 4)
  seek(con, 4 + LH + 4)
  N = readBin(con, n = 1, what = "integer", size = 4)  # number of samples
  samplesInGeno = rep(0, N)
  
  # cycle for all samples to extract IDs
  for(i in 1:N){
    LS = readBin(con, n = 1, what = "integer", size = 2)
    sample = readChar(con, nchars = LS)
    samplesInGeno[i] = sample
  }
  
  # close connection
  close(con)
  
  return(samplesInGeno)
}


extract_genoIndex_condition = function(condition, markerInfo, genoType){
   if(condition != ""){
       	condition_original = unlist(strsplit(condition, ","))
	#if(!is.null(weight_cond)){
	#	weight_original = unlist(strsplit(weight_cond, ","))
	#}
	conditionDat = data.frame(SNP = condition_original, condIndex = seq(1,length(condition_original)))
	conditionDat = setDT(conditionDat)
		#print("head(markerInfo)")
		#print(head(markerInfo))
		#CHROM POS    ID genoIndex       ID2 genoIndex_prev

	if(genoType != "vcf"){
		markerInfo_conditionDat = merge(conditionDat, markerInfo, by.x="SNP", by.y="ID", sort = F, all.x=T)
		#print(markerInfo_conditionDat)
		#print("markerInfo_conditionDat")
		# SNP condIndex CHROM POS genoIndex  ID2 genoIndex_prev
		if(!is.null(markerInfo$ID2)){
		        #setnames(markerInfo, "genoIndex", "genoIndex2")	
			#colnames(markerInfo)[which(colnames(markerInfo) == "genoIndex")] = "genoIndex2"

			markerInfo_conditionDat = merge(markerInfo_conditionDat, markerInfo, by.x="SNP", by.y="ID2", sort = F, all.x=T)
			#print("markerInfo_conditionDat")
			#print(markerInfo_conditionDat)

			        #SNP condIndex CHROM.x POS.x genoIndex.x  ID2 genoIndex_prev.x CHROM.y
				#1: 1:33:A:C         1    <NA>    NA          NA <NA>               NA       1
				   #POS.y   ID genoIndex.y genoIndex_prev.y
				  # 1:    33 rs33       22895            23399
			setnames(markerInfo_conditionDat, "genoIndex.x", "genoIndex")
			setnames(markerInfo_conditionDat, "genoIndex.y", "genoIndex2")
			posNA = which(is.na(markerInfo_conditionDat$genoIndex) & !is.na(markerInfo_conditionDat$genoIndex2))
        		if(length(posNA) != 0){
                		markerInfo_conditionDat$genoIndex[posNA] = markerInfo_conditionDat$genoIndex2[posNA]
        		}
			markerInfo_conditionDat$genoIndex2 = NULL

		}
		markerInfo_conditionDat = markerInfo_conditionDat[which(!is.na(markerInfo_conditionDat$genoIndex)), , drop=F]

		markerInfo_conditionDat = markerInfo_conditionDat[order(markerInfo_conditionDat$condIndex), ]	
		#markerInfo_conditionDat = markerInfo_conditionDat[order(markerInfo_conditionDat$genoIndex), ]

		#markerInfo_conditionDat = markerInfo_conditionDat[which()]
       		#posInd = which(markerInfo_conditionDat %in% condition_original)
		#if(length(posInd) == length(condition_original)){
		if(nrow(markerInfo_conditionDat) == length(condition_original)){
			cond_genoIndex = markerInfo_conditionDat$genoIndex
			cond_genoIndex_prev = markerInfo_conditionDat$genoIndex_prev
       			#cond_genoIndex = genoIndex[posInd]
		}else{

			stop(length(condition_original)-nrow(markerInfo_conditionDat), " conditioning markers are not found in the geno file. Please Check.\n")	
		}	
       }else{
	        condition_group_line = paste(c("condition", condition_original), collapse = "\t")

		if(length(condition_original) == 1){
			fakem = strsplit(condition_original, split=":")[[1]]
			fakemb = paste(c(fakem[1], as.numeric(fakem[2])+2, "N", "N"), collapse=":")
			condition_group_linetemp = paste(c(condition_group_line, fakemb), collapse = "\t")
			set_iterator_inVcf(condition_group_linetemp, "1", 1, 250000000)
		}else{	

			set_iterator_inVcf(condition_group_line, "1", 1, 250000000)
		}
		cond_genoIndex = rep("0", length(condition_original))
		cond_genoIndex_prev = rep("0", length(condition_original))
       }
   }else{
	stop("condition is empty!")
   }
   return(list(cond_genoIndex = cond_genoIndex, cond_genoIndex_prev = cond_genoIndex_prev ))
}
 
