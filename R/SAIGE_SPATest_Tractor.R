SAIGE.Admixed = function(mu,
                        OutputFile,
                        MACCutoff_to_CollapseUltraRare,
                        groupFile,
                        genoType,
                        markerInfo,
                        traitType,
                        isImputation,
                        isCondition,
                        groups_per_chunk,
                        isOverWriteOutput,
			is_no_weight_in_groupTest,
                        chrom,
                        is_fastTest,
                        pval_cutoff_for_fastTest,
                        is_output_moreDetails) {
  OutputFileIndex = NULL
  if (is.null(OutputFileIndex))
    OutputFileIndex = paste0(OutputFile, ".index")

  outList = checkOutputFile(OutputFile, OutputFileIndex, "Region", 1, isOverWriteOutput) # Check 'Util.R'

  indexChunk = outList$indexChunk
  Start = outList$Start
  End = outList$End

  cat("Start ", Start, "\n")
  cat("End ", End, "\n")


  if (End)
  {
    message = paste0(
      "The analysis has been completed in earlier analysis. Results are saved in '",
      OutputFile,
      "'. ",
      "If you want to change parameters and restart the analysis, please use another 'OutputFile'."
    )
    return(message)
  }

  isappend = FALSE
  if (!Start) {
    isappend = TRUE
  }

  n = length(mu) #sample size

  r.corr = 0
  out.method = SKAT:::SKAT_Check_Method(method = "optimal.adj", r.corr = 0)
  method = out.method$method
  r.corr = out.method$r.corr
  regionTestType = "SKAT-O"

  region_list = checkGroupFile(groupFile)
  nRegions = region_list$nRegions
  is_weight_included = region_list$is_weight_included
  if (is_no_weight_in_groupTest & is_weight_included) {
    stop("is_no_weight_in_groupTest = TRUE but weights are found in the group file.\n")
  }

  if (is_no_weight_in_groupTest) {
    cat("No weights are used in the group test\n")
  }

  if (is_weight_included) {
    nline_per_gene = 3
  } else {
    nline_per_gene = 2
  }

  gf = file(groupFile, "r")

  cat("indexChunk is ", indexChunk, "\n")

  skipline = indexChunk * nline_per_gene
  if (indexChunk > 0 & indexChunk < nRegions) {
    for (k in 1:skipline) {
      marker_group_line_temp = readLines(gf, n = 1)
      rm(marker_group_line_temp)
    }
  }

  chrom1 = "FakeCHR"
    gc()
  num_region = 0
  mth = 0

  numberRegionsInChunk = 0
  cat("indexChunk ", indexChunk, "\n")
  cat("nRegions ", nRegions, "\n")
  cth_chunk_to_output = 1
  i = indexChunk + 1
  nEachChunk = numberRegionsInChunkR

 while (i <= nRegions) {
    #for(i in (indexChunk+1):nRegions){
    if (mth ==  numberRegionsInChunk) {
      if (i + groups_per_chunk > nRegions) {
        nregions_ro_read = nRegions - i + 1
      } else {
        nregions_ro_read = groups_per_chunk
      }
      nlinetoread = nregions_ro_read * nline_per_gene
      marker_group_line = readLines(gf, n = nlinetoread)
      RegionList = SAIGE.getRegionList_new(marker_group_line,
                                           nline_per_gene,
                                           annolist,
                                           markerInfo,
                                           chrom)
      #print(RegionList)
      cat("Read in ",
          nregions_ro_read,
          " region(s) from the group file.\n")
      mth = 0
      #numberRegionsInChunk = length(RegionList)
      numberRegionsInChunk = nregions_ro_read
    }

    #mth = mth + 1
    if (!is.null(RegionList)) {
	mainAdmixedInCPP(
	  RegionList, 
          genoType,
          OutputFile,
          traitType,
          n,
          isImputation,
          is_fastTest,
          is_output_moreDetails
        )
    }else{
        #if(!is.null(region)){
        cat(regionName, " is empty.\n")
    }
    mth = numberRegionsInChunk
    i = i + numberRegionsInChunk
 }

  message = "Analysis done!"
  message = paste0(message,
                   " The association test results for admixed samples have been saved to '",
                   OutputFile,
                   "'.")

  return(message)
}

