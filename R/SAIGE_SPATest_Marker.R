SAIGE.Marker = function(traitType,
			genoType,
			bgenFileIndex,
                        idstoIncludeFile,
                        rangestoIncludeFile,
			genoIndex_prev,
                        genoIndex,
			CHROM,
                        OutputFile,
                        OutputFileIndex = NULL,
                        nMarkersEachChunk, 
			isMoreOutput,
			isImputation,
			isFirth,
			LOCO,
			chrom,
			isCondition,
			isOverWriteOutput, 
			isAnyInclude,
			isGbyE)
{

  if(is.null(OutputFileIndex))
    OutputFileIndex = paste0(OutputFile, ".index")
  
  #genoType = objGeno$genoType

  outIndex = checkOutputFile(OutputFile, OutputFileIndex, "Marker", format(nMarkersEachChunk, scientific=F), isOverWriteOutput)    # this function is in 'Util.R'
  outIndex = outIndex$indexChunk
  isappend = FALSE
  if(outIndex != 1){
    cat("Restart the analysis from chunk:\t", outIndex, "\n")
    isappend = TRUE
  }  

  isOpenOutFile_single = openOutfile_single(traitType, isImputation, isappend, isMoreOutput, isGbyE)

  if(!isOpenOutFile_single){
    stop("Output file ", OutputFile, " can't be opened\n")
  }

  ## set up an object for genotype
  if(genoType == "plink" || genoType == "pgen"){
      #markerInfo = objGeno$markerInfo
    if(LOCO){
      genoIndex = genoIndex[which(CHROM == chrom)]
      if(!is.null(genoIndex_prev)){
	genoIndex_prev = genoIndex_prev[which(CHROM == chrom)]
      }
      CHROM = CHROM[which(CHROM == chrom)]  
      #markerInfo = markerInfo[which(markerInfo$CHROM == chrom),]  
    }
    #CHROM = markerInfo$CHROM
    #genoIndex = markerInfo$genoIndex
    ##only for one chrom
    # all markers were split into multiple chunks,
    #print(markerInfo[1:10,])
    #print(genoIndex[1:10])    

    genoIndexList = splitMarker(genoIndex, genoIndex_prev, nMarkersEachChunk, CHROM);
    nChunks = length(genoIndexList)

    if(nChunks == 0){
          stop("No markers on chrom ", chrom, " are found\n")
    }


    cat("Number of all markers to test:\t", length(genoIndex), "\n")
    cat("Number of markers in each chunk:\t", nMarkersEachChunk, "\n")
    cat("Number of chunks for all markers:\t", nChunks, "\n")
    if(outIndex > nChunks){
      cat("The analysis has been finished! Please delete ", OutputFileIndex, " if the analysis needs to be run again or set --is_overwrite_output=TRUE\n")
      is_marker_test = FALSE 
    }else{
      is_marker_test = TRUE
      i = outIndex    
    }
    
  }else if(genoType == "vcf"){
   if(!isAnyInclude){	  
    if(chrom == ""){
      stop("chrom needs to be specified for single-variant assoc tests when using VCF as input\n")
    }else{
      set_iterator_inVcf("", chrom, 1, 250000000)
    }
   }
    if(outIndex > 1){
	move_forward_iterator_Vcf(outIndex*nMarkersEachChunk)    
    }
    isVcfEnd =  check_Vcf_end()
    if(!isVcfEnd){
    	#outIndex = 1
    	genoIndex = rep("0", nMarkersEachChunk) 
    	genoIndex_prev = rep("0", nMarkersEachChunk) 
	#nChunks = outIndex + 1
	is_marker_test = TRUE
        i = outIndex
    }else{
	is_marker_test = FALSE    
	stop("No markers are left in VCF")
    }
  }else if(genoType == "bgen"){
    IDsToInclude=NULL
    RangesToInclude=NULL

    db_con <- RSQLite::dbConnect(RSQLite::SQLite(), bgenFileIndex)
    on.exit(RSQLite::dbDisconnect(db_con), add = TRUE)

   anyInclude=FALSE
   if(idstoIncludeFile != ""){
    IDsToInclude = data.table::fread(idstoIncludeFile, header = F, colClasses = c("character"), data.table=F)
    if(ncol(IDsToInclude) != 1){
      stop("'idstoIncludeFile' of ", idstoIncludeFile, " should only include one column.")
    }else{
      IDsToInclude = IDsToInclude[,1]
      cat(length(IDsToInclude), " marker IDs are found in the idstoIncludeFile file\n")
      anyInclude=TRUE  
    }
   }
  if(rangestoIncludeFile != ""){
    RangesToInclude = data.table::fread(rangestoIncludeFile, header = F, colClasses = c("character", "numeric", "numeric"), data.table=F)
    if(ncol(RangesToInclude) != 3){
      stop("rangestoIncludeFile should only include three columns.")
    }else{
      cat(nrow(RangesToInclude), " ranges are in rangestoIncludeFile\n")
    }
    anyInclude=TRUE
    colnames(RangesToInclude) = c("CHROM", "START", "END")
    #print("RangesToInclude")
    #print(RangesToInclude)
    if(nrow(RangesToInclude) > 1){cat("WARNING, only the first line of ",rangestoIncludeFile, " will be used\n")}
  } 


if(FALSE){
    if(LOCO | chrom != ""){

	if(!anyInclude){
   		query <- "SELECT chromosome, position, file_start_position, size_in_bytes
             		FROM Variant
                       		WHERE chromosome = ?"
   	}else{
        	if(!is.null(IDsToInclude) > 0){
           		id_list <- paste0("'", IDsToInclude, "'", collapse = ", ")
        		query <- paste("
  				SELECT chromosome, position, rsid, allele1, allele2, file_start_position, size_in_bytes
  				FROM Variant
  				WHERE chromosome = ?
    				AND (CONCAT(chromosome, ':', position, ':', allele1, ':', allele2)) IN (", id_list, ")", sep = "")
        	}else if(!is.null(RangesToInclude)){
			CHROM=RangesToInclude$CHROM[1]
			START=RangesToInclude$START[1]
			END=RangesToInclude$END[1]
			query <- paste("
                                SELECT chromosome, position, rsid, allele1, allele2, file_start_position, size_in_bytes
                                FROM Variant
                                WHERE chromosome = ?
                                AND  chromosome =", CHROM, 
                		"AND position >= ", START, 
                		"AND position <= ", END)
		}
   	}	
    }else{

       if(!anyInclude){
                query <- "SELECT chromosome, position, file_start_position, size_in_bytes
                        FROM Variant"
        }else{
                if(!is.null(IDsToInclude) > 0){
                        id_list <- paste0("'", IDsToInclude, "'", collapse = ", ")
                        query <- paste("
                                SELECT chromosome, position, rsid, allele1, allele2, file_start_position, size_in_bytes
                                FROM Variant
                                WHERE (CONCAT(chromosome, ':', position, ':', allele1, ':', allele2)) IN (", id_list, ")", sep = "")
                }else if(!is.null(RangesToInclude)){
                        CHROM=RangesToInclude$CHROM[1]
                        START=RangesToInclude$START[1]
                        END=RangesToInclude$END[1]
                        query <- paste("
                                SELECT chromosome, position, rsid, allele1, allele2, file_start_position, size_in_bytes
                                FROM Variant
                                WHERE  chromosome =", CHROM,
                                "AND position >= ", START,
                                "AND position <= ", END)
                }
        }
    }
   
}

   if(LOCO){

        if(!anyInclude){
                query <- "SELECT COUNT(*) AS total_variants
                        FROM Variant
                                WHERE chromosome = ?"
        }else{
                if(!is.null(IDsToInclude)){
                        id_list <- paste0("'", IDsToInclude, "'", collapse = ", ")
                        query <- paste("SELECT COUNT(*) AS total_variants
                                FROM Variant
                                WHERE chromosome = ?
                                AND (CONCAT(chromosome, ':', position, ':', allele1, ':', allele2)) IN (", id_list, ")", 
				"OR (CONCAT(chromosome, ':', position, ':', allele2, ':', allele1)) IN (", id_list, ")", 
				"OR (rsid) IN (", id_list, ")",
				sep = "")

                }else if(!is.null(RangesToInclude)){
                        CHROM=RangesToInclude$CHROM[1]
                        START=RangesToInclude$START[1]
                        END=RangesToInclude$END[1]

			query <- paste0(
  			"SELECT COUNT(*) AS total_variants
   			FROM Variant
   			WHERE chromosome = ? 
   			AND chromosome = '", CHROM, "' 
   			AND position >= ", START, "
   			AND position <= ", END
			)

                }
        }
       count_result <- RSQLite::dbGetQuery(db_con, query, params = list(chrom))
    }else{
       if(!anyInclude){
                query <- "  SELECT COUNT(*) AS total_variants
                        FROM Variant"
        }else{
                if(!is.null(IDsToInclude)){
                        id_list <- paste0("'", IDsToInclude, "'", collapse = ", ")

		   query <- paste("SELECT COUNT(*) AS total_variants
                                FROM Variant
                                WHERE (CONCAT(chromosome, ':', position, ':', allele1, ':', allele2)) IN (", id_list, ")",
                                "OR (CONCAT(chromosome, ':', position, ':', allele2, ':', allele1)) IN (", id_list, ")",
                                "OR (rsid) IN (", id_list, ")",
                                sep = "")	

                        #query <- paste("SELECT COUNT(*) AS total_variants
                        #        FROM Variant
                        #        WHERE (CONCAT(chromosome, ':', position, ':', allele1, ':', allele2)) IN (", id_list, ")", sep = "")
                }else if(!is.null(RangesToInclude)){


                        CHROM=RangesToInclude$CHROM[1]
                        START=RangesToInclude$START[1]
                        END=RangesToInclude$END[1]
		
			# Construct the SQL query
			query <- paste0("
    			SELECT COUNT(*) AS total_variants
    			FROM Variant
    			WHERE chromosome = '", CHROM, "'
			AND position >= ", START, "
    			AND position <= ", END
			)


                }
        }
	count_result <- RSQLite::dbGetQuery(db_con, query)
    } 

    total_variants <- count_result$total_variants
    cat("total_variants ", total_variants, "\n")
    nChunks = ceiling(total_variants / nMarkersEachChunk)

    if(nChunks == 0){
          stop("No markers on chrom ", chrom, " are found\n")
    }


    cat("Number of all markers to test:\t", total_variants, "\n")
    cat("Number of markers in each chunk:\t", nMarkersEachChunk, "\n")
    cat("Number of chunks for all markers:\t", nChunks, "\n")
    if(outIndex > nChunks){
      cat("The analysis has been finished! Please delete ", OutputFileIndex, " if the analysis needs to be run again or set --is_overwrite_output=TRUE\n")
      is_marker_test = FALSE
    }else{
      is_marker_test = TRUE
      i = outIndex
    }

    offset <- 0
    if(nChunks == 1){
	nMarkersEachChunk = total_variants 
    }
  }
  #chrom = "InitialChunk"
  #set_flagSparseGRM_cur_SAIGE_org()
  while(is_marker_test){
  #for(i in outIndex:nChunks)
  #{
#time_left = system.time({


    if(genoType == "plink" || genoType == "pgen"){	
      tempList = genoIndexList[[i]]
      genoIndex = as.character(format(tempList$genoIndex, scientific = FALSE))
      tempChrom = tempList$chrom
      if(!is.null(genoIndex_prev)){
      	genoIndex_prev = as.character(format(tempList$genoIndex_prev, scientific = FALSE))
      }else{
	genoIndex_prev = c("-1")
      }
    }



    if(genoType == "bgen"){
	print("chrom")
	print(chrom)
    if(LOCO | chrom != ""){
	#chrom = "01"
        if(!anyInclude){
                #query <-  paste("SELECT chromosome, position, file_start_position, size_in_bytes
                 #       FROM Variant
                 #               WHERE chromosome = ?
	#			LIMIT", nMarkersEachChunk, "OFFSET", offset)
                query <-  paste("SELECT chromosome, position, file_start_position, size_in_bytes
                        FROM Variant
			WHERE chromosome= ?
				LIMIT", nMarkersEachChunk, "OFFSET", offset)

        }else{
                if(!is.null(IDsToInclude)){
                        id_list <- paste0("'", IDsToInclude, "'", collapse = ", ")

                  
		query <- paste("SELECT chromosome, position, rsid, allele1, allele2, file_start_position, size_in_bytes
                                FROM Variant
                                WHERE chromosome = ?
                                AND (CONCAT(chromosome, ':', position, ':', allele1, ':', allele2)) IN (", id_list, ")",
                                "OR (CONCAT(chromosome, ':', position, ':', allele2, ':', allele1)) IN (", id_list, ")",
                                "OR (rsid) IN (", id_list, ")",
                                sep = "")


                        #query <- paste("SELECT chromosome, position, rsid, allele1, allele2, file_start_position, size_in_bytes
                         #       FROM Variant
                         #       WHERE chromosome = ?
                         #       AND (CONCAT(chromosome, ':', position, ':', allele1, ':', allele2)) IN (", id_list, ")
			#	LIMIT", nMarkersEachChunk, "OFFSET", offset, sep = " ")
                }else if(!is.null(RangesToInclude)){
                        CHROM=RangesToInclude$CHROM[1]
                        START=RangesToInclude$START[1]
                        END=RangesToInclude$END[1]

query <- paste0("SELECT chromosome, position, rsid, allele1, allele2, file_start_position, size_in_bytes
                 FROM Variant
                 WHERE chromosome = ?
		 AND chromosome = '", CHROM, "'
                 AND position >= ", START, "
                 AND position <= ", END, "
                 LIMIT ", nMarkersEachChunk, " 
                 OFFSET ", offset)

                }
        }
	query_result <- RSQLite::dbGetQuery(db_con, query, params = list(chrom))
	#query_result <- RSQLite::dbGetQuery(db_con, query)
    }else{

       if(!anyInclude){
                query <- paste("SELECT chromosome, position, file_start_position, size_in_bytes
                        FROM Variant
			LIMIT", nMarkersEachChunk, "OFFSET", offset)
        }else{
                if(!is.null(IDsToInclude) > 0){
                        id_list <- paste0("'", IDsToInclude, "'", collapse = ", ")
                        #query <- paste("SELECT chromosome, position, rsid, allele1, allele2, file_start_position, size_in_bytes
                         #       FROM Variant
                         #       WHERE (CONCAT(chromosome, ':', position, ':', allele1, ':', allele2)) IN (", id_list, ")
			#	LIMIT", nMarkersEachChunk, "OFFSET", offset, sep = " ")

                query <- paste("SELECT chromosome, position, rsid, allele1, allele2, file_start_position, size_in_bytes
                                FROM Variant
                                WHERE (CONCAT(chromosome, ':', position, ':', allele1, ':', allele2)) IN (", id_list, ")",
                                "OR (CONCAT(chromosome, ':', position, ':', allele2, ':', allele1)) IN (", id_list, ")",
                                "OR (rsid) IN (", id_list, ")
				LIMIT", nMarkersEachChunk, "OFFSET", offset, sep = " ")			


                }else if(!is.null(RangesToInclude)){
                        CHROM=RangesToInclude$CHROM[1]
                        START=RangesToInclude$START[1]
                        END=RangesToInclude$END[1]

query <- paste0("SELECT chromosome, position, rsid, allele1, allele2, file_start_position, size_in_bytes
                 FROM Variant
                 WHERE chromosome = '", CHROM, "' 
                 AND position >= ", START, "
                 AND position <= ", END, "
                 LIMIT ", nMarkersEachChunk, " 
                 OFFSET ", offset)

                }
        }
	query_result <- RSQLite::dbGetQuery(db_con, query)
    }

	#print(head(query_result))
	#markerInfo <- RSQLite::fetch(query_result, n = -1)
	markerInfo=query_result
	genoIndex_prev = as.character(markerInfo$file_start_position + markerInfo$size_in_bytes)
	genoIndex = as.character(markerInfo$file_start_position)

	#print(length(genoIndex))
	#print(length(genoIndex_prev))
	offset <- offset + nMarkersEachChunk
	#RSQLite::dbClearResult(query_result)
    }
    #print("tempList")
    #print(tempList)
    #print(tempList$genoIndex)
#})
#print("time_left")
#print(time_left)
    #print("genoIndex here")
    #print(genoIndex)
    # set up objects that do not change for different variants
    #if(tempChrom != chrom){
    #  setMarker("SAIGE", objNull, control, chrom, Group, ifOutGroup)
    #  chrom = tempChrom
    #}
    #ptm <- proc.time()
    #print(ptm)
    #print("gc()")
    #print(gc())
    if(genoType != "vcf"){
      cat(paste0("(",Sys.time(),") ---- Analyzing Chunk ", i, "/", nChunks, ": chrom ", chrom," ---- \n"))
    }else{
      cat(paste0("(",Sys.time(),") ---- Analyzing Chunk ", i, " :  chrom ", chrom," ---- \n"))
    }	    
    # main function to calculate summary statistics for markers in one chunk
   
   #resMarker = as.data.frame(mainMarkerInCPP(genoType, traitType, genoIndex_prev, genoIndex, isMoreOutput, isImputation)) 
   #resMarker = resMarker[which(!is.na(resMarker$BETA)), ]

  mainMarkerInCPP(genoType, traitType, genoIndex_prev, genoIndex, isMoreOutput, isImputation, isFirth)

    #timeoutput=system.time({writeOutputFile(Output = list(resMarker),
  #if(nrow(resMarker) > 0){

  if(genoType == "vcf"){
    isEnd_Output =  check_Vcf_end()
  }else{
    isEnd_Output = (i==nChunks)	
  }
  
  #}

  writeOutputFileIndex(OutputFileIndex = OutputFileIndex,
			AnalysisType = "Marker",
			nEachChunk = format(nMarkersEachChunk, scientific=F),
			indexChunk = i,
			Start = (i==1),
			End = isEnd_Output)

  #writeOutputFile(Output = list(resMarker),
  #                  OutputFile = list(OutputFile),
  #                  OutputFileIndex = OutputFileIndex,
  #                  AnalysisType = "Marker",
  #                  nEachChunk = format(nMarkersEachChunk, scientific=F),
  #                  indexChunk = i,
  #                  Start = (i==1),
  #                  End = isEnd_Output)

                    #End = (i==nChunks))})
    #print("timeoutput")
    #print(timeoutput)
    ptm <- proc.time()
    print(ptm)
    gc()
    #rm(resMarker)


    
    i = i + 1
  if(genoType == "vcf"){
    isVcfEnd =  check_Vcf_end()
    cat("isVcfEnd ", isVcfEnd, "\n")
    if(isVcfEnd){
	is_marker_test = FALSE	     
    }
  }else{
    if(i > nChunks){
      is_marker_test = FALSE
    }	    
  }
	  
  } #while(is_marker_test){

  if(genoType == "bgen"){
    RSQLite::dbDisconnect(db_con)
  }
  # information to users
  output = paste0("Analysis done! The results have been saved to '", OutputFile,"'.")

  return(output)
}


splitMarker = function(genoIndex, genoIndex_prev, nMarkersEachChunk, CHROM)
{
  genoIndexList = list()
  iTot = 1;

  uCHROM = unique(CHROM)
  for(chrom in uCHROM){
    pos = which(CHROM == chrom)
    gIdx = genoIndex[pos]
    if(!is.null(genoIndex_prev)){
	gIndprev = genoIndex_prev[pos]
    }
    M = length(gIdx)

    idxStart = seq(1, M, nMarkersEachChunk)
    idxEnd = idxStart + nMarkersEachChunk - 1

    nChunks = length(idxStart)
    idxEnd[nChunks] = M

    for(i in 1:nChunks){
      idxMarker = idxStart[i]:idxEnd[i]
      genoIndexList[[iTot]] = list(chrom = chrom,
                                   genoIndex = gIdx[idxMarker])
      if(!is.null(genoIndex_prev)){
	genoIndexList[[iTot]]$genoIndex_prev = gIndprev[idxMarker] 	     }

      iTot = iTot + 1;
    }
  }

  return(genoIndexList)
}
