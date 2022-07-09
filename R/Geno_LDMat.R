checkGenoInput_LDmat = function(bgenFile = "",
                 bgenFileIndex = "",
                 vcfFile = "",
                 vcfFileIndex = "",
                 vcfField = "DS",
                 savFile = "",
                 savFileIndex = "",
                 sampleFile = "",
                 bedFile="",
                 bimFile="",
                 famFile=""){

    #if(is.null(sampleInModel)){
    #    stop("sampleInModel is not specified.")
    #}
    # genotype data
    if (bgenFile != "") {
        Check_File_Exist(bgenFile, "bgenFile")
        Check_File_Exist(bgenFileIndex, "bgenFileIndex")
        dosageFileType = "bgen"
    }else if (vcfFile != "") {
        Check_File_Exist(vcfFile, "vcfFile")

        if (!grepl("\\.sav$", vcfFile) && !file.exists(paste(vcfFile, ".csi", sep = ""))) {
            stop("ERROR! vcfFileIndex ", paste(vcfFile, ".csi", sep = ""), " does not exist\n")
        }
        dosageFileType = "vcf"
    }else if (savFile != "") {
        Check_File_Exist(savFile, "savFile")
        Check_File_Exist(savFileIndex, "savFileIndex")
        vcfFile = savFile
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



setGenoInput_LDmat = function(bgenFile = "",
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

  dosageFileType = checkGenoInput_LDmat(bgenFile = bgenFile,
                 bgenFileIndex = bgenFileIndex,
                 vcfFile = vcfFile,
                 vcfFileIndex = vcfFileIndex,
                 vcfField = vcfField,
                 savFile = savFile,
                 savFileIndex = savFileIndex,
                 bedFile = bedFile,
                 bimFile = bimFile,
                 famFile = famFile)
                 #sampleInModel = sampleInModel)

  #cat("length(sampleInModel) is ", length(sampleInModel), "\n")
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
    if(is.null(sampleInModel)){
      sampleInModel = data.table:::fread(famFile, header = F, colClasses = list(character = 1:4), select=c(2), data.table=F)[,1]	    
      if(any(duplicated(sampleInModel))){
        stop("ERROR: dupliated sample IDs are found in fam file.\n")
       }	      
    }
    #SampleIDs = updateSampleIDs(SampleIDs, samplesInGeno)
    #markerInfo$ID = paste0(markerInfo$CHROM,":", markerInfo$POS ,"_", markerInfo$REF, "/", markerInfo$ALT)
    setPLINKobjInCPP(bimFile, famFile, bedFile, sampleInModel, AlleleOrder)
  }

  ########## ----------  BGEN format ---------- ##########

  if(dosageFileType == "bgen"){
  #cat("length(sampleInModel) b is ", length(sampleInModel), "\n")


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


   if(is.null(sampleInModel)){
      sampleInModel = samplesInGeno	   
      if(any(duplicated(sampleInModel))){
        stop("ERROR: dupliated sample IDs are found in the bgen sample file.\n")
       }
    }
  #cat("length(sampleInModel) c is ", length(sampleInModel), "\n")

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

    if(!is.null(markerInfo$ID2)){

        posRows = which((markerInfo$ID %in% IDsToInclude) | (markerInfo$ID2 %in% IDsToInclude))
    }else{
        posRows = which((markerInfo$ID %in% IDsToInclude))

    }
    if(length(posRows) != 0)
      cat(length(posRows), " markers in idstoIncludeFile are found in the geno/dosage file.\n")
      markersInclude = c(markersInclude, markerInfo$ID[posRows])
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
   if(nrow(RangesToInclude)){
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

  if(anyInclude)
    markerInfo = subset(markerInfo, ID %in% markersInclude)

#  if(anyExclude)
#    markerInfo = subset(markerInfo, !ID %in% markersExclude)

#  anyQueue = anyInclude | anyExclude

  if(dosageFileType == "vcf"){

   if(is.null(sampleInModel)){
      sampleInModel = as.character(NULL)
    }

    setVCFobjInCPP(vcfFile, vcfFileIndex, vcfField, t_SampleInModel = sampleInModel)
    if(!is.null(IDsToInclude)){
      SNPlist = paste(c("set1", IDsToInclude), collapse = "\t")
      in_chrom="fake_chrom"
      in_beg_pd=1
      in_end_pd=250000000
      set_iterator_inVcf(SNPlist, in_chrom, in_beg_pd, in_end_pd)
    }

    if(!is.null(RangesToInclude)){
      if(length(CHROM1) > 1){
        stop("We do not support query with multiple regions for vcf file. Please only include one region in the ", rangestoIncludeFile, "\n")
      }else{
        inSNPlist=""
        in_chrom=CHROM1[1]
        in_beg_pd=START[1]
        in_end_pd=END[1]
        set_iterator_inVcf(inSNPlist, in_chrom, in_beg_pd, in_end_pd)
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

  genoList = list(markerInfo = markerInfo, genoType = dosageFileType)
  return(genoList)
}
