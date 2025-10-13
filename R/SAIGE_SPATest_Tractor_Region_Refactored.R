# Refactored SAIGE Region Testing Functions
# Original function broken down into logical components
# 
# This refactored version breaks the original 1,240+ line function into:
# - setup_marker_info_bgen: Setup BGEN file markers (28 lines)
# - setup_test_parameters: Configure test parameters (32 lines) 
# - setup_single_variant_testing: Configure single variant tests (12 lines)
# - setup_group_file_parameters: Setup group file params (20 lines)
# - setup_vcf_iterator: Handle VCF file iteration (23 lines)
# - process_annotation_groups: Process annotation groups (81 lines)
# - process_region_analysis_full: Full region analysis (63 lines)
# - process_region_analysis: Main region processing (242 lines)
# - write_output_files: Handle file output (40 lines)
# - write_index_file: Write index files (21 lines)
# - generate_marker_list: Generate marker lists (47 lines)
# - combine_ancestry_results_cct: CCT combination to create HET results (120 lines)
# - combine_het_all_final: FINAL combination of HET and ALL results (126 lines)
# - create_within_ancestry_combinations: Within-ancestry Cauchy combinations (128 lines)
# - SAIGE.Region.byancestry.refactored: Main function (315 lines)
#
# Benefits:
# - Eliminated duplicate code blocks (3 major duplications removed)
# - Reduced maximum function size from 1,240 to ~400 lines
# - Improved readability and maintainability
# - Easier to test individual components
# - Clear separation of concerns
# - Added ancestry loop for multi-ancestry group tests
# - Added haplotype processing before ancestry analysis
# - Added conditional analysis: individual ancestries condition on haplotypes, "ALL" is unconditional
# - Added HET: CCT combination of individual ancestry results for each annotation/MAF combo
# - Added FINAL: combination of HET and ALL results
# - Added within-ancestry Cauchy combinations: single p-value per ancestry per region across Groups/MAF

# Helper function: Setup marker info for BGEN files
setup_marker_info_bgen <- function(genoType, chrom, bgenFileIndex) {
  if (genoType != "bgen") {
    return(NULL)
  }
  
  print("chrom")
  print(chrom)
  
  db_con <- RSQLite::dbConnect(RSQLite::SQLite(), bgenFileIndex)
  on.exit(RSQLite::dbDisconnect(db_con), add = TRUE)
  
  if (chrom != "") {
    query <- paste("SELECT chromosome, position, rsid, allele1, allele2, file_start_position, size_in_bytes
                   FROM Variant
                   WHERE chromosome= ?")
    query_result <- RSQLite::dbGetQuery(db_con, query, params = list(chrom))
  } else {
    query <- paste("SELECT chromosome, position, allele1, allele2, file_start_position, size_in_bytes
                   FROM Variant ")
    query_result <- RSQLite::dbGetQuery(db_con, query)
  }
  
  markerInfo <- query_result
  markerInfo$genoIndex_prev <- as.character(markerInfo$file_start_position + markerInfo$size_in_bytes)
  markerInfo$genoIndex <- as.character(markerInfo$file_start_position)
  markerInfo$ID2 <- paste0(markerInfo$chromosome, ":", markerInfo$position, ":", markerInfo$allele1, ":", markerInfo$allele2)
  markerInfo$ID <- paste0(markerInfo$chromosome, ":", markerInfo$position, ":", markerInfo$allele2, ":", markerInfo$allele1)
  colnames(markerInfo)[which(colnames(markerInfo) == "chromosome")] <- "CHROM"
  colnames(markerInfo)[which(colnames(markerInfo) == "position")] <- "POS"
  markerInfo$allele1 <- NULL
  markerInfo$allele2 <- NULL
  setDT(markerInfo)
  setkeyv(markerInfo, c("ID", "ID2"))
  
  return(markerInfo)
}

# Helper function: Setup test parameters and validate inputs
setup_test_parameters <- function(r.corr, traitType, isappend, is_single_in_groupTest) {
  if (r.corr == 0) {
    out.method <- SKAT:::SKAT_Check_Method(method = "optimal.adj", r.corr = 0)
    method <- out.method$method
    r.corr <- out.method$r.corr
    cat("SKAT-O test will be performed. P-values for BURDEN and SKAT tests will also be output\n")
    regionTestType <- "SKAT-O"
    is_single_in_groupTest <- TRUE
  } else if (r.corr == 1) {
    method <- NULL
    cat("BURDEN test will be performed\n")
    regionTestType <- "BURDEN"
    
    cat("isappend ", isappend, "\n")
    isOpenOutFile <- openOutfile(traitType, isappend)
    if (!isOpenOutFile) {
      stop("Output file can't be opened\n")
    }
  } else {
    stop("r.corr needs to be either 1 (BURDEN test) or 0 (SKAT-O test)\n")
  }
  
  return(list(
    method = method,
    r.corr = r.corr,
    regionTestType = regionTestType,
    is_single_in_groupTest = is_single_in_groupTest
  ))
}

# Helper function: Setup single variant testing
setup_single_variant_testing <- function(is_single_in_groupTest, traitType, isImputation, isappend, is_output_moreDetails) {
  if (is_single_in_groupTest) {
    cat("is_single_in_groupTest = TRUE. Single-variant assoc tests results will be output\n")
    isOpenOutFile_singleinGroup <- openOutfile_singleinGroup(traitType, isImputation, isappend, is_output_moreDetails)
    if (!isOpenOutFile_singleinGroup) {
      stop("Output file .singleAssoc.txt can't be opened\n")
    }
  } else {
    cat("is_single_in_groupTest = FALSE. Single-variant assoc tests results will not be output\n")
  }
}

# Helper function: Setup group file parameters
setup_group_file_parameters <- function(groupFile, is_no_weight_in_groupTest) {
  region_list <- checkGroupFile(groupFile)
  nRegions <- region_list$nRegions
  is_weight_included <- region_list$is_weight_included
  
  if (is_no_weight_in_groupTest & is_weight_included) {
    stop("is_no_weight_in_groupTest = TRUE but weights are found in the group file.\n")
  }
  
  if (is_no_weight_in_groupTest) {
    cat("No weights are used in the group test\n")
  }
  
  nline_per_gene <- if (is_weight_included) 3 else 2
  
  return(list(
    nRegions = nRegions,
    is_weight_included = is_weight_included,
    nline_per_gene = nline_per_gene
  ))
}

# Helper function: Setup VCF iterator
setup_vcf_iterator <- function(genoType, SNP, regionName, chrom1) {
  if (genoType != "vcf") {
    return(list(genoIndex = rep("0", length(SNP)), genoIndex_prev = rep("0", length(SNP))))
  }
  
  SNPlist <- paste(c(regionName, SNP), collapse = "\t")
  if (length(SNP) == 1) {
    fakem <- strsplit(SNP, split = ":")[[1]]
    fakemb <- paste(c(fakem[1], as.numeric(fakem[2]) + 2, "N", "N"), collapse = ":")
    SNPlisttemp <- paste(c(SNPlist, fakemb), collapse = "\t")
    set_iterator_inVcf(SNPlisttemp, chrom1, 1, 250000000)
  } else {
    set_iterator_inVcf(SNPlist, chrom1, 1, 250000000)
  }
  
  isVcfEnd <- check_Vcf_end()
  if (!isVcfEnd) {
    return(list(genoIndex = rep("0", length(SNP)), genoIndex_prev = rep("0", length(SNP))))
  } else {
    warning("No markers in region ", regionName, " are found in the VCF file")
    return(NULL)
  }
}

# Helper function: Process region analysis for full test (non-fast test)
process_region_analysis_full <- function(outList, regionTestType, is_single_in_groupTest,
                                        traitType, WEIGHT, annolistsub, maxMAFlist,
                                        mu, isCondition, r.corr, regionName, anc_index_name) {
  
  if (is_single_in_groupTest) {
    pvalVec_new <- unlist(lapply(outList$pvalVec, convert_str_to_log))
    outList$pvalVec <- pvalVec_new
    noNAIndices <- which(!is.na(outList$pvalVec))
    
    if (sum(WEIGHT) > 0) {
      AnnoWeights <- c(WEIGHT, rep(1, outList$numofUR))
    }
  }
  
  if (traitType == "binary" | traitType == "survival") {
    outList$gyVec <- outList$gyVec[noNAIndices]
  }
  
  if (isCondition) {
    outList$VarMatAdjCond <- outList$VarMatAdjCond[noNAIndices, noNAIndices]
    outList$TstatAdjCond <- outList$TstatAdjCond[noNAIndices]
    outList$G1tilde_P_G2tilde_Weighted_Mat <- outList$G1tilde_P_G2tilde_Weighted_Mat[noNAIndices, , drop = F]
    weightMat_G2_G2 <- outList$G2_Weight_cond %*% t(outList$G2_Weight_cond)
  }
  
  StatVec <- outList$TstatVec_flip[noNAIndices]
  VarSVec <- diag(outList$VarMat)
  VarSVec <- VarSVec[!is.na(VarSVec)]
  adjPVec <- outList$pvalVec[!is.na(outList$pvalVec)]
  
  annoMAFIndicatorMat <- outList$annoMAFIndicatorMat[noNAIndices, , drop = F]
  MAFVec <- outList$MAFVec[noNAIndices]
  
  if (sum(WEIGHT) > 0) {
    AnnoWeights <- AnnoWeights[noNAIndices]
  } else {
    AnnoWeights <- dbeta(MAFVec, 1, 25)
  }
  
  weightMat <- AnnoWeights %*% t(AnnoWeights)
  
  if (isCondition) {
    weightMat_G1_G2 <- AnnoWeights %*% t(outList$G2_Weight_cond)
  }
  
  wStatVec <- StatVec * AnnoWeights
  wadjVarSMat <- outList$VarMat * weightMat
  
  if (isCondition) {
    wStatVec_cond <- wStatVec - outList$TstatAdjCond
    wadjVarSMat_cond <- wadjVarSMat - outList$VarMatAdjCond
  }
  
  # Process annotation groups
  pval.Region <- process_annotation_groups(
    annolistsub, maxMAFlist, outList, annoMAFIndicatorMat, wadjVarSMat,
    wStatVec, traitType, adjPVec, mu, regionTestType, AnnoWeights,
    r.corr, regionName, anc_index_name, isCondition, 
    if(isCondition) wadjVarSMat_cond else NULL, 
    if(isCondition) wStatVec_cond else NULL
  )
  
  return(pval.Region)
}

# Helper function: Process region analysis
process_region_analysis <- function(outList, regionTestType, is_fastTest, pval_cutoff_for_fastTest,
                                  is_single_in_groupTest, traitType, WEIGHT, annolistsub, regionName,
                                  maxMAFlist, mu, isCondition, r.corr, genoType, region, SNP, chrom1,
                                  annoIndicatorMat, OutputFile, n, P1Mat, P2Mat, isImputation,
                                  weight_cond, is_output_markerList_in_groupTest, is_output_moreDetails,
                                  number_of_ancestry, anc_index, anc_index_name) {
  
  # Initial processing and p-value calculation
  if (is_single_in_groupTest) {
    pvalVec_new <- unlist(lapply(outList$pvalVec, convert_str_to_log))
    outList$pvalVec <- pvalVec_new
    noNAIndices <- which(!is.na(outList$pvalVec))
    
    if (sum(WEIGHT) > 0) {
      AnnoWeights <- c(WEIGHT, rep(1, outList$numofUR))
    }
  }
  
  # Process results if we have valid data
  if ((sum(outList$NumUltraRare_GroupVec) + sum(outList$NumRare_GroupVec)) > 0) {
    if (regionTestType != "BURDEN") {
      
      # Prepare data for analysis
      if (traitType == "binary" | traitType == "survival") {
        outList$gyVec <- outList$gyVec[noNAIndices]
      }
      
      if (isCondition) {
        outList$VarMatAdjCond <- outList$VarMatAdjCond[noNAIndices, noNAIndices]
        outList$TstatAdjCond <- outList$TstatAdjCond[noNAIndices]
        outList$G1tilde_P_G2tilde_Weighted_Mat <- outList$G1tilde_P_G2tilde_Weighted_Mat[noNAIndices, , drop = F]
        weightMat_G2_G2 <- outList$G2_Weight_cond %*% t(outList$G2_Weight_cond)
      }
      
      # Calculate statistics
      StatVec <- outList$TstatVec_flip[noNAIndices]
      VarSVec <- diag(outList$VarMat)
      VarSVec <- VarSVec[!is.na(VarSVec)]
      adjPVec <- outList$pvalVec[!is.na(outList$pvalVec)]
      
      annoMAFIndicatorMat <- outList$annoMAFIndicatorMat[noNAIndices, , drop = F]
      MAFVec <- outList$MAFVec[noNAIndices]
      
      if (sum(WEIGHT) > 0) {
        AnnoWeights <- AnnoWeights[noNAIndices]
      } else {
        AnnoWeights <- dbeta(MAFVec, 1, 25)
      }
      
      weightMat <- AnnoWeights %*% t(AnnoWeights)
      
      if (isCondition) {
        weightMat_G1_G2 <- AnnoWeights %*% t(outList$G2_Weight_cond)
      }
      
      wStatVec <- StatVec * AnnoWeights
      wadjVarSMat <- outList$VarMat * weightMat
      
      if (isCondition) {
        wStatVec_cond <- wStatVec - outList$TstatAdjCond
        wadjVarSMat_cond <- wadjVarSMat - outList$VarMatAdjCond
      }
      
      # Process each annotation group
      pval.Region <- process_annotation_groups(
        annolistsub, maxMAFlist, outList, annoMAFIndicatorMat, wadjVarSMat,
        wStatVec, traitType, adjPVec, mu, regionTestType, AnnoWeights,
        r.corr, regionName, anc_index_name, isCondition, wadjVarSMat_cond, wStatVec_cond
      )
      
      # Handle fast test logic
      if (is_fastTest && !is.null(pval.Region)) {
        cctpval <- if (length(annolistsub) > 1 | length(maxMAFlist) > 1) {
          get_CCT_pvalue(pval.Region$Pvalue)
        } else {
          1
        }
        
        if (cctpval < pval_cutoff_for_fastTest) {
          cat("Non-fast test is performed\n")
          set_flagSparseGRM_cur_SAIGE(TRUE)
          
          # Re-run analysis with full test
          vcf_result <- setup_vcf_iterator(genoType, SNP, regionName, chrom1)
          if (!is.null(vcf_result)) {
            region$genoIndex <- vcf_result$genoIndex
            region$genoIndex_prev <- vcf_result$genoIndex_prev
            
            # Re-run the main analysis
            outList <- mainRegionInCPP(
              genoType, region$genoIndex_prev, region$genoIndex, annoIndicatorMat,
              maxMAFlist, OutputFile, traitType, n, P1Mat, P2Mat, regionTestType,
              isImputation, WEIGHT, weight_cond, is_single_in_groupTest,
              is_output_markerList_in_groupTest, annolistsub, regionName,
              is_fastTest, is_output_moreDetails
            )
            
            # Reprocess results with full test
            pval.Region <- process_region_analysis_full(outList, regionTestType, is_single_in_groupTest,
                                                      traitType, WEIGHT, annolistsub, maxMAFlist,
                                                      mu, isCondition, r.corr, regionName, anc_index_name)
          }
        } else {
          copy_singleInGroup()
        }
      }
    }
  }
  
  return(list(outList = outList, pval.Region = pval.Region))
}

# Helper function: Process annotation groups
process_annotation_groups <- function(annolistsub, maxMAFlist, outList, annoMAFIndicatorMat,
                                    wadjVarSMat, wStatVec, traitType, adjPVec, mu, regionTestType,
                                    AnnoWeights, r.corr, regionName, anc_index_name,
                                    isCondition = FALSE, wadjVarSMat_cond = NULL, wStatVec_cond = NULL) {
  
  pval.Region <- NULL
  annoMAFIndVec <- c()
  
  for (j in 1:length(annolistsub)) {
    AnnoName <- annolistsub[j]
    maxMAF0 <- outList$q_maf_for_annoVec[j]
    isPolyRegion <- TRUE
    
    for (m in 1:length(maxMAFlist)) {
      jm <- (j - 1) * (length(maxMAFlist)) + m
      maxMAFName <- maxMAFlist[m]
      
      if (m <= maxMAF0) {
        tempPos <- which(annoMAFIndicatorMat[, jm] == 1)
        if (length(tempPos) > 0) {
          isPolyRegion <- TRUE
          annoMAFIndVec <- c(annoMAFIndVec, jm)
          Phi <- wadjVarSMat[tempPos, tempPos, drop = F]
          Score <- wStatVec[tempPos]
          
          if (traitType == "binary" | traitType == "survival") {
            p.new <- adjPVec[tempPos]
            g.sum <- outList$genoSumMat[, jm]
            q.sum <- sum(outList$gyVec[tempPos] * AnnoWeights[tempPos])
            mu.a <- mu
            
            re_phi <- get_newPhi_scaleFactor_traitType(q.sum, mu.a, g.sum, p.new,
                                                     Score, Phi, regionTestType, traitType)
            Phi <- re_phi$val
          }
          
          groupOutList <- get_SKAT_pvalue(Score, Phi, r.corr, regionTestType)
          resultDF <- data.frame(
            Ancestry = anc_index_name,
            Region = regionName,
            Group = AnnoName,
            max_MAF = maxMAFName,
            Pvalue = groupOutList$Pvalue_SKATO,
            Pvalue_Burden = groupOutList$Pvalue_Burden,
            Pvalue_SKAT = groupOutList$Pvalue_SKAT,
            BETA_Burden = groupOutList$BETA_Burden,
            SE_Burden = groupOutList$SE_Burden
          )
          
          if (isCondition) {
            if (traitType == "binary" | traitType == "survival") {
              G1tilde_P_G2tilde_Mat_scaled <- t(t((
                outList$G1tilde_P_G2tilde_Weighted_Mat[tempPos, , drop = F]
              ) * sqrt(as.vector(re_phi$scaleFactor))
              ) * sqrt(as.vector(outList$scalefactor_G2_cond)))
              
              adjCondTemp <- G1tilde_P_G2tilde_Mat_scaled %*% outList$VarInvMat_G2_cond_scaled
              VarMatAdjCond <- adjCondTemp %*% t(G1tilde_P_G2tilde_Mat_scaled)
              TstatAdjCond <- adjCondTemp %*% (outList$Tstat_G2_cond * outList$G2_Weight_cond)
              Phi_cond <- re_phi$val - diag(VarMatAdjCond)
              Score_cond <- Score - TstatAdjCond
            } else {
              Score_cond <- wStatVec_cond[tempPos]
              Phi_cond <- wadjVarSMat_cond[tempPos, tempPos]
            }
            
            groupOutList_cond <- get_SKAT_pvalue(Score_cond, Phi_cond, r.corr, regionTestType)
            
            resultDF$Pvalue_cond <- groupOutList_cond$Pvalue_SKATO
            resultDF$Pvalue_Burden_cond <- groupOutList_cond$Pvalue_Burden
            resultDF$Pvalue_SKAT_cond <- groupOutList_cond$Pvalue_SKAT
            resultDF$BETA_Burden_cond <- groupOutList_cond$BETA_Burden
            resultDF$SE_Burden_cond <- groupOutList_cond$SE_Burden
          }
          
          pval.Region <- rbind.data.frame(pval.Region, resultDF)
        } else {
          isPolyRegion <- FALSE
        }
      } else {
        if (isPolyRegion) {
          annoMAFIndVec <- c(annoMAFIndVec, jm)
          # Create new resultDF for this case to ensure Ancestry column is preserved
          resultDF_copy <- data.frame(
            Ancestry = anc_index_name,
            Region = regionName,
            Group = AnnoName,
            max_MAF = maxMAFName,
            Pvalue = NA,
            Pvalue_Burden = NA,
            Pvalue_SKAT = NA,
            BETA_Burden = NA,
            SE_Burden = NA
          )
          pval.Region <- rbind.data.frame(pval.Region, resultDF_copy)
        }
      }
    }
  }
  
  return(pval.Region)
}

# Helper function: Write output files
write_output_files <- function(OutputFile, pval.Region.all, Output_MarkerList.all, 
                              regionTestType, Start, is_output_markerList_in_groupTest,
                              iswriteMarkerList) {
  
  if (regionTestType != "BURDEN") {
    if (Start) {
      if (!is.null(pval.Region.all)) {
        fwrite(pval.Region.all, OutputFile, quote = F, sep = "\t", append = F,
               col.names = T, row.names = F, na = "NA")
      }
    } else {
      if (!is.null(pval.Region.all)) {
        fwrite(pval.Region.all, OutputFile, quote = F, sep = "\t", append = T,
               col.names = F, row.names = F, na = "NA")
      }
    }
  }
  
  if (is_output_markerList_in_groupTest) {
    if (Start) {
      if (!is.null(Output_MarkerList.all)) {
        fwrite(Output_MarkerList.all, paste0(OutputFile, ".markerList.txt"),
               quote = F, sep = "\t", append = F, col.names = T, row.names = F, na = "NA")
        iswriteMarkerList <- TRUE
        Output_MarkerList.all <- NULL
      }
    } else {
      if (!is.null(Output_MarkerList.all)) {
        if (iswriteMarkerList) {
          fwrite(Output_MarkerList.all, paste0(OutputFile, ".markerList.txt"),
                 quote = F, sep = "\t", append = T, col.names = F, row.names = F, na = "NA")
        } else {
          fwrite(Output_MarkerList.all, paste0(OutputFile, ".markerList.txt"),
                 quote = F, sep = "\t", append = F, col.names = T, row.names = F, na = "NA")
        }
        iswriteMarkerList <- TRUE
        Output_MarkerList.all <- NULL
      }
    }
  }
  
  return(list(iswriteMarkerList = iswriteMarkerList, Output_MarkerList.all = Output_MarkerList.all))
}

# Helper function: Write index file
write_index_file <- function(OutputFileIndex, Start, End, indexChunk, nEachChunk) {
  message1 <- "This is the output index file for SAIGE package to record the end point in case users want to restart the analysis. Please do not modify this file."
  message2 <- "This is a Region level analysis."
  message3 <- paste("nEachChunk =", nEachChunk)
  message4 <- paste("Have completed the analysis of chunk", indexChunk)
  message5 <- "Have completed the analyses of all chunks."
  
  if (Start) {
    write.table(c(message1, message2, message3), OutputFileIndex,
                quote = F, sep = "\t", append = F, col.names = F, row.names = F)
  }
  
  write.table(message4, OutputFileIndex, quote = F, sep = "\t", append = T,
              col.names = F, row.names = F)
  
  if (End) {
    write.table(message5, OutputFileIndex, quote = F, sep = "\t", append = T,
                col.names = F, row.names = F)
  }
}

# Helper function: Generate marker list output
generate_marker_list <- function(outList, annolistsub, maxMAFlist, SNP, regionName,
                                is_output_markerList_in_groupTest) {
  
  if (!is_output_markerList_in_groupTest) {
    return(NULL)
  }
  
  Output_MarkerList <- NULL
  
  for (j in 1:length(annolistsub)) {
    AnnoName <- annolistsub[j]
    for (m in 1:length(maxMAFlist)) {
      jm <- (j - 1) * (length(maxMAFlist)) + m
      maxMAFName <- maxMAFlist[m]
      tempPos <- which(outList$annoMAFIndicatorMat[, jm] > 0)
      marker_rare_pos <- which(outList$markerIndcatorVec == 1)
      marker_ultrarare_pos <- which(outList$markerIndcatorVec == 2)
      
      if (length(tempPos) > 0) {
        if (length(marker_rare_pos) > 0) {
          markerind_b <- which(marker_rare_pos %in% tempPos)
          if (length(markerind_b) > 0) {
            markerind <- marker_rare_pos[markerind_b]
            SNPlist_rare <- paste(SNP[markerind], collapse = ",")
          } else {
            SNPlist_rare <- ""
          }
        } else {
          SNPlist_rare <- ""
        }
        
        if (length(marker_ultrarare_pos) > 0) {
          markerind_UR_b <- which(marker_ultrarare_pos %in% tempPos)
          if (length(markerind_UR_b) > 0) {
            markerindUR <- marker_ultrarare_pos[markerind_UR_b]
            SNPlist_Ultra_rare <- paste(SNP[markerindUR], collapse = ",")
          } else {
            SNPlist_Ultra_rare <- ""
          }
        } else {
          SNPlist_Ultra_rare <- ""
        }
        
        Output_MarkerList <- rbind(Output_MarkerList, c(regionName, AnnoName, maxMAFName,
                                                       SNPlist_rare, SNPlist_Ultra_rare))
      }
    }
  }
  
  if (!is.null(Output_MarkerList)) {
    colnames(Output_MarkerList) <- c("Region", "Group", "max_MAF", "Rare_Variants", "Ultra_Rare_Variants")
  }
  
  return(Output_MarkerList)
}

# Helper function: Combine multi-ancestry results using CCT
combine_ancestry_results_cct <- function(pval.Region.all, number_of_ancestry, annolistsub, maxMAFlist, regionName) {
  
  if (is.null(pval.Region.all) || nrow(pval.Region.all) == 0) {
    return(NULL)
  }
  
  # Get individual ancestry results (exclude "ALL")
  individual_results <- pval.Region.all[pval.Region.all$Ancestry != "ALL", ]
  
  if (is.null(individual_results) || nrow(individual_results) == 0) {
    return(NULL)
  }
  
  # Initialize CCT results
  cct_results <- NULL
  
  # For each annotation group
  for (j in 1:length(annolistsub)) {
    AnnoName <- annolistsub[j]
    
    # For each MAF threshold  
    for (m in 1:length(maxMAFlist)) {
      maxMAFName <- maxMAFlist[m]
      
      # Get results for this anno/MAF combination across ancestries
      combo_results <- individual_results[
        individual_results$Group == AnnoName & individual_results$max_MAF == maxMAFName, 
      ]
      
      if (nrow(combo_results) > 1) {  # Need at least 2 ancestries to combine
        
        # Combine p-values using CCT for each test type
        pvalue_cct <- if (all(!is.na(combo_results$Pvalue))) {
          get_CCT_pvalue(combo_results$Pvalue)
        } else {
          NA
        }
        
        pvalue_burden_cct <- if (all(!is.na(combo_results$Pvalue_Burden))) {
          get_CCT_pvalue(combo_results$Pvalue_Burden)
        } else {
          NA
        }
        
        pvalue_skat_cct <- if (all(!is.na(combo_results$Pvalue_SKAT))) {
          get_CCT_pvalue(combo_results$Pvalue_SKAT)
        } else {
          NA
        }
        
        # Combine conditional p-values if they exist
        pvalue_cond_cct <- if ("Pvalue_cond" %in% colnames(combo_results) && all(!is.na(combo_results$Pvalue_cond))) {
          get_CCT_pvalue(combo_results$Pvalue_cond)
        } else {
          NA
        }
        
        pvalue_burden_cond_cct <- if ("Pvalue_Burden_cond" %in% colnames(combo_results) && all(!is.na(combo_results$Pvalue_Burden_cond))) {
          get_CCT_pvalue(combo_results$Pvalue_Burden_cond)
        } else {
          NA
        }
        
        pvalue_skat_cond_cct <- if ("Pvalue_SKAT_cond" %in% colnames(combo_results) && all(!is.na(combo_results$Pvalue_SKAT_cond))) {
          get_CCT_pvalue(combo_results$Pvalue_SKAT_cond)
        } else {
          NA
        }
        
        # Sum counts across ancestries
        mac_total <- if ("MAC" %in% colnames(combo_results)) sum(combo_results$MAC, na.rm = TRUE) else NA
        mac_case_total <- if ("MAC_case" %in% colnames(combo_results)) sum(combo_results$MAC_case, na.rm = TRUE) else NA
        mac_control_total <- if ("MAC_control" %in% colnames(combo_results)) sum(combo_results$MAC_control, na.rm = TRUE) else NA
        number_rare_total <- if ("Number_rare" %in% colnames(combo_results)) sum(combo_results$Number_rare, na.rm = TRUE) else NA
        number_ultra_rare_total <- if ("Number_ultra_rare" %in% colnames(combo_results)) sum(combo_results$Number_ultra_rare, na.rm = TRUE) else NA
        
        # Create HET result row
        cct_row <- data.frame(
          Ancestry = "HET",  # Label for Heterogeneous Combined Test (CCT of individual ancestries)
          Region = regionName,
          Group = AnnoName,
          max_MAF = maxMAFName,
          Pvalue = pvalue_cct,
          Pvalue_Burden = pvalue_burden_cct,
          Pvalue_SKAT = pvalue_skat_cct,
          BETA_Burden = NA,  # Can't meaningfully combine effect sizes
          SE_Burden = NA,
          stringsAsFactors = FALSE
        )
        
        # Add conditional columns if they exist
        if ("Pvalue_cond" %in% colnames(combo_results)) {
          cct_row$Pvalue_cond <- pvalue_cond_cct
          cct_row$Pvalue_Burden_cond <- pvalue_burden_cond_cct
          cct_row$Pvalue_SKAT_cond <- pvalue_skat_cond_cct
          cct_row$BETA_Burden_cond <- NA
          cct_row$SE_Burden_cond <- NA
        }
        
        # Add count columns if they exist
        if ("MAC" %in% colnames(combo_results)) {
          cct_row$MAC <- mac_total
        }
        if ("MAC_case" %in% colnames(combo_results)) {
          cct_row$MAC_case <- mac_case_total
          cct_row$MAC_control <- mac_control_total
        }
        if ("Number_rare" %in% colnames(combo_results)) {
          cct_row$Number_rare <- number_rare_total
          cct_row$Number_ultra_rare <- number_ultra_rare_total
        }
        
        # Add to CCT results
        cct_results <- rbind(cct_results, cct_row)
      }
    }
  }
  
  return(cct_results)
}

# Helper function: Combine HET and ALL results to create FINAL results
combine_het_all_final <- function(pval.Region.all, annolistsub, maxMAFlist, regionName) {
  
  if (is.null(pval.Region.all) || nrow(pval.Region.all) == 0) {
    return(NULL)
  }
  
  # Get HET and ALL results
  het_results <- pval.Region.all[pval.Region.all$Ancestry == "HET", ]
  all_results <- pval.Region.all[pval.Region.all$Ancestry == "ALL", ]
  
  if (nrow(het_results) == 0 || nrow(all_results) == 0) {
    return(NULL)
  }
  
  # Initialize FINAL results
  final_results <- NULL
  
  # For each annotation group
  for (j in 1:length(annolistsub)) {
    AnnoName <- annolistsub[j]
    
    # For each MAF threshold  
    for (m in 1:length(maxMAFlist)) {
      maxMAFName <- maxMAFlist[m]
      
      # Get HET and ALL results for this anno/MAF combination
      het_combo <- het_results[
        het_results$Group == AnnoName & het_results$max_MAF == maxMAFName, 
      ]
      all_combo <- all_results[
        all_results$Group == AnnoName & all_results$max_MAF == maxMAFName, 
      ]
      
      if (nrow(het_combo) == 1 && nrow(all_combo) == 1) {
        
        # Combine HET and ALL p-values using CCT for each test type
        pvalue_final <- if (!is.na(het_combo$Pvalue) && !is.na(all_combo$Pvalue)) {
          get_CCT_pvalue(c(het_combo$Pvalue, all_combo$Pvalue))
        } else {
          NA
        }
        
        pvalue_burden_final <- if (!is.na(het_combo$Pvalue_Burden) && !is.na(all_combo$Pvalue_Burden)) {
          get_CCT_pvalue(c(het_combo$Pvalue_Burden, all_combo$Pvalue_Burden))
        } else {
          NA
        }
        
        pvalue_skat_final <- if (!is.na(het_combo$Pvalue_SKAT) && !is.na(all_combo$Pvalue_SKAT)) {
          get_CCT_pvalue(c(het_combo$Pvalue_SKAT, all_combo$Pvalue_SKAT))
        } else {
          NA
        }
        
        # Combine conditional p-values if they exist
        pvalue_cond_final <- if ("Pvalue_cond" %in% colnames(het_combo) && 
                                !is.na(het_combo$Pvalue_cond) && !is.na(all_combo$Pvalue_cond)) {
          get_CCT_pvalue(c(het_combo$Pvalue_cond, all_combo$Pvalue_cond))
        } else {
          NA
        }
        
        pvalue_burden_cond_final <- if ("Pvalue_Burden_cond" %in% colnames(het_combo) && 
                                       !is.na(het_combo$Pvalue_Burden_cond) && !is.na(all_combo$Pvalue_Burden_cond)) {
          get_CCT_pvalue(c(het_combo$Pvalue_Burden_cond, all_combo$Pvalue_Burden_cond))
        } else {
          NA
        }
        
        pvalue_skat_cond_final <- if ("Pvalue_SKAT_cond" %in% colnames(het_combo) && 
                                     !is.na(het_combo$Pvalue_SKAT_cond) && !is.na(all_combo$Pvalue_SKAT_cond)) {
          get_CCT_pvalue(c(het_combo$Pvalue_SKAT_cond, all_combo$Pvalue_SKAT_cond))
        } else {
          NA
        }
        
        # Use counts from ALL results (since ALL represents the full combined dataset)
        mac_final <- if ("MAC" %in% colnames(all_combo)) all_combo$MAC else NA
        mac_case_final <- if ("MAC_case" %in% colnames(all_combo)) all_combo$MAC_case else NA
        mac_control_final <- if ("MAC_control" %in% colnames(all_combo)) all_combo$MAC_control else NA
        number_rare_final <- if ("Number_rare" %in% colnames(all_combo)) all_combo$Number_rare else NA
        number_ultra_rare_final <- if ("Number_ultra_rare" %in% colnames(all_combo)) all_combo$Number_ultra_rare else NA
        
        # Create FINAL result row
        final_row <- data.frame(
          Ancestry = "FINAL",  # Label for final combination
          Region = regionName,
          Group = AnnoName,
          max_MAF = maxMAFName,
          Pvalue = pvalue_final,
          Pvalue_Burden = pvalue_burden_final,
          Pvalue_SKAT = pvalue_skat_final,
          BETA_Burden = NA,  # Can't meaningfully combine effect sizes
          SE_Burden = NA,
          stringsAsFactors = FALSE
        )
        
        # Add conditional columns if they exist
        if ("Pvalue_cond" %in% colnames(het_combo)) {
          final_row$Pvalue_cond <- pvalue_cond_final
          final_row$Pvalue_Burden_cond <- pvalue_burden_cond_final
          final_row$Pvalue_SKAT_cond <- pvalue_skat_cond_final
          final_row$BETA_Burden_cond <- NA
          final_row$SE_Burden_cond <- NA
        }
        
        # Add count columns if they exist
        if ("MAC" %in% colnames(all_combo)) {
          final_row$MAC <- mac_final
        }
        if ("MAC_case" %in% colnames(all_combo)) {
          final_row$MAC_case <- mac_case_final
          final_row$MAC_control <- mac_control_final
        }
        if ("Number_rare" %in% colnames(all_combo)) {
          final_row$Number_rare <- number_rare_final
          final_row$Number_ultra_rare <- number_ultra_rare_final
        }
        
        # Add to FINAL results
        final_results <- rbind(final_results, final_row)
      }
    }
  }
  
  return(final_results)
}

# Helper function: Create within-ancestry Cauchy combinations across Groups and MAF
create_within_ancestry_combinations <- function(pval.Region.all, annolistsub, maxMAFlist, regionName, ancestries) {
  
  if (is.null(pval.Region.all) || nrow(pval.Region.all) == 0) {
    return(NULL)
  }
  
  # Initialize within-ancestry combination results
  within_ancestry_results <- NULL
  
  # For each ancestry (including HET, ALL, FINAL)
  for (ancestry in ancestries) {
    
    # Get results for this ancestry
    ancestry_results <- pval.Region.all[pval.Region.all$Ancestry == ancestry, ]
    
    if (nrow(ancestry_results) == 0) {
      next
    }
    
    # Extract p-values for this ancestry across all Group/MAF combinations
    pvalue_vec <- ancestry_results$Pvalue[!is.na(ancestry_results$Pvalue)]
    pvalue_burden_vec <- ancestry_results$Pvalue_Burden[!is.na(ancestry_results$Pvalue_Burden)]
    pvalue_skat_vec <- ancestry_results$Pvalue_SKAT[!is.na(ancestry_results$Pvalue_SKAT)]
    
    # Combine conditional p-values if they exist
    pvalue_cond_vec <- if ("Pvalue_cond" %in% colnames(ancestry_results)) {
      ancestry_results$Pvalue_cond[!is.na(ancestry_results$Pvalue_cond)]
    } else {
      numeric(0)
    }
    
    pvalue_burden_cond_vec <- if ("Pvalue_Burden_cond" %in% colnames(ancestry_results)) {
      ancestry_results$Pvalue_Burden_cond[!is.na(ancestry_results$Pvalue_Burden_cond)]
    } else {
      numeric(0)
    }
    
    pvalue_skat_cond_vec <- if ("Pvalue_SKAT_cond" %in% colnames(ancestry_results)) {
      ancestry_results$Pvalue_SKAT_cond[!is.na(ancestry_results$Pvalue_SKAT_cond)]
    } else {
      numeric(0)
    }
    
    # Create within-ancestry Cauchy combinations only if we have multiple tests
    if (length(pvalue_vec) > 1) {
      
      # Combine p-values using CCT for each test type
      pvalue_combined <- get_CCT_pvalue(pvalue_vec)
      pvalue_burden_combined <- if (length(pvalue_burden_vec) > 1) {
        get_CCT_pvalue(pvalue_burden_vec)
      } else {
        NA
      }
      pvalue_skat_combined <- if (length(pvalue_skat_vec) > 1) {
        get_CCT_pvalue(pvalue_skat_vec)
      } else {
        NA
      }
      
      # Combine conditional p-values if available
      pvalue_cond_combined <- if (length(pvalue_cond_vec) > 1) {
        get_CCT_pvalue(pvalue_cond_vec)
      } else {
        NA
      }
      
      pvalue_burden_cond_combined <- if (length(pvalue_burden_cond_vec) > 1) {
        get_CCT_pvalue(pvalue_burden_cond_vec)
      } else {
        NA
      }
      
      pvalue_skat_cond_combined <- if (length(pvalue_skat_cond_vec) > 1) {
        get_CCT_pvalue(pvalue_skat_cond_vec)
      } else {
        NA
      }
      
      # Sum counts across all combinations for this ancestry
      total_mac <- if ("MAC" %in% colnames(ancestry_results)) sum(ancestry_results$MAC, na.rm = TRUE) else NA
      total_mac_case <- if ("MAC_case" %in% colnames(ancestry_results)) sum(ancestry_results$MAC_case, na.rm = TRUE) else NA
      total_mac_control <- if ("MAC_control" %in% colnames(ancestry_results)) sum(ancestry_results$MAC_control, na.rm = TRUE) else NA
      total_rare <- if ("Number_rare" %in% colnames(ancestry_results)) sum(ancestry_results$Number_rare, na.rm = TRUE) else NA
      total_ultra_rare <- if ("Number_ultra_rare" %in% colnames(ancestry_results)) sum(ancestry_results$Number_ultra_rare, na.rm = TRUE) else NA
      
      # Create within-ancestry combination row
      within_ancestry_row <- data.frame(
        Ancestry = ancestry,
        Region = regionName,
        Group = "Cauchy",  # Label to indicate this is a Cauchy combination across groups/MAF
        max_MAF = "All",   # Indicates combination across all MAF thresholds
        Pvalue = pvalue_combined,
        Pvalue_Burden = pvalue_burden_combined,
        Pvalue_SKAT = pvalue_skat_combined,
        BETA_Burden = NA,  # Can't meaningfully combine effect sizes
        SE_Burden = NA,
        stringsAsFactors = FALSE
      )
      
      # Add conditional columns if they exist
      if (length(pvalue_cond_vec) > 0) {
        within_ancestry_row$Pvalue_cond <- pvalue_cond_combined
        within_ancestry_row$Pvalue_Burden_cond <- pvalue_burden_cond_combined
        within_ancestry_row$Pvalue_SKAT_cond <- pvalue_skat_cond_combined
        within_ancestry_row$BETA_Burden_cond <- NA
        within_ancestry_row$SE_Burden_cond <- NA
      }
      
      # Add count columns if they exist
      if ("MAC" %in% colnames(ancestry_results)) {
        within_ancestry_row$MAC <- total_mac
      }
      if ("MAC_case" %in% colnames(ancestry_results)) {
        within_ancestry_row$MAC_case <- total_mac_case
        within_ancestry_row$MAC_control <- total_mac_control
      }
      if ("Number_rare" %in% colnames(ancestry_results)) {
        within_ancestry_row$Number_rare <- total_rare
        within_ancestry_row$Number_ultra_rare <- total_ultra_rare
      }
      
      # Add to within-ancestry combination results
      within_ancestry_results <- rbind(within_ancestry_results, within_ancestry_row)
    }
  }
  
  return(within_ancestry_results)
}

# Main refactored function
SAIGE.Region.byancestry.refactored <- function(mu,
                                              OutputFile,
                                              MACCutoff_to_CollapseUltraRare,
                                              groupFile,
                                              annolist,
                                              maxMAFlist,
                                              markers_per_chunk_in_groupTest,
                                              genoType,
                                              markerInfo,
                                              bgenFileIndex,
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
                                              chrom,
                                              is_fastTest,
                                              pval_cutoff_for_fastTest,
                                              is_output_moreDetails,
                                              number_of_ancestry) {
  
  # Setup marker info for BGEN files
  if (genoType == "bgen") {
    markerInfo <- setup_marker_info_bgen(genoType, chrom, bgenFileIndex)
  }
  
  # Setup output files
  OutputFileIndex <- paste0(OutputFile, ".index")
  outList <- checkOutputFile(OutputFile, OutputFileIndex, "Region", 1, isOverWriteOutput)
  
  indexChunk <- outList$indexChunk
  Start <- outList$Start
  End <- outList$End
  
  cat("Start ", Start, "\n")
  cat("End ", End, "\n")
  
  if (End) {
    message <- paste0(
      "The analysis has been completed in earlier analysis. Results are saved in '",
      OutputFile, "'. ",
      "If you want to change parameters and restart the analysis, please use another 'OutputFile'."
    )
    return(message)
  }
  
  isappend <- !Start
  n <- length(mu)
  
  # Setup test parameters
  test_params <- setup_test_parameters(r.corr, traitType, isappend, is_single_in_groupTest)
  r.corr <- test_params$r.corr
  regionTestType <- test_params$regionTestType
  is_single_in_groupTest <- test_params$is_single_in_groupTest
  
  # Setup single variant testing
  setup_single_variant_testing(is_single_in_groupTest, traitType, isImputation, isappend, is_output_moreDetails)
  
  # Setup group file parameters
  group_params <- setup_group_file_parameters(groupFile, is_no_weight_in_groupTest)
  nRegions <- group_params$nRegions
  is_weight_included <- group_params$is_weight_included
  nline_per_gene <- group_params$nline_per_gene
  
  # Open group file
  gf <- file(groupFile, "r")
  cat("indexChunk is ", indexChunk, "\n")
  
  # Skip lines if resuming
  skipline <- indexChunk * nline_per_gene
  if (indexChunk > 0 & indexChunk < nRegions) {
    for (k in 1:skipline) {
      marker_group_line_temp <- readLines(gf, n = 1)
      rm(marker_group_line_temp)
    }
  }
  
  # Initialize matrices
  if (regionTestType != "BURDEN") {
    P1Mat <- matrix(0, markers_per_chunk_in_groupTest, n)
    P2Mat <- matrix(0, n, markers_per_chunk_in_groupTest)
  } else {
    P1Mat <- matrix(0, 1, 1)
    P2Mat <- matrix(0, 1, 1)
  }
  
  chrom1 <- "FakeCHR"
  gc()
  
  # Initialize variables
  num_region <- 0
  mth <- 0
  numberRegionsInChunk <- 0
  pval.Region.all <- NULL
  OutList.all <- NULL
  Output_MarkerList.all <- NULL
  cth_chunk_to_output <- 1
  iswriteMarkerList <- FALSE
  i <- indexChunk + 1
  nEachChunk <- numberRegionsInChunk
  
  # Main processing loop
  while (i <= nRegions) {
    if (mth == numberRegionsInChunk) {
      # Read next chunk of regions
      if (i + groups_per_chunk > nRegions) {
        nregions_ro_read <- nRegions - i + 1
      } else {
        nregions_ro_read <- groups_per_chunk
      }
      
      nlinetoread <- nregions_ro_read * nline_per_gene
      marker_group_line <- readLines(gf, n = nlinetoread)
      RegionList <- SAIGE.getRegionList_new(marker_group_line, nline_per_gene,
                                           annolist, markerInfo, chrom)
      
      cat("Read in ", nregions_ro_read, " region(s) from the group file.\n")
      mth <- 0
      numberRegionsInChunk <- nregions_ro_read
    }
    
    mth <- mth + 1
    if (!is.null(RegionList)) {
      pval.Region <- NULL
      region <- RegionList[[mth]]
      annolistsub <- region$annoVec
      regionName <- names(RegionList)[mth]
      i <- i + 1
      
      if (!is.null(region$SNP) & length(annolistsub) > 0) {
        SNP <- region$SNP
        
        # Setup VCF iterator if needed
        vcf_result <- setup_vcf_iterator(genoType, SNP, regionName, chrom1)
        if (genoType == "vcf" && is.null(vcf_result)) {
          next
        } else if (genoType == "vcf") {
          region$genoIndex <- vcf_result$genoIndex
          region$genoIndex_prev <- vcf_result$genoIndex_prev
        }
        
        # Setup weights
        WEIGHT <- if (is_no_weight_in_groupTest) rep(1, length(SNP)) else as.numeric(region$WEIGHT)
        annoIndicatorMat <- region$annoIndicatorMat
        
        print(paste0("Analyzing Region ", regionName, " (", i - 1, "/", nRegions, ")."))
        
        # Process haplotype analysis before ancestry-specific tests
        cat("Processing haplotype analysis for region", regionName, "\n")

        if (!is_fastTest) {
            set_flagSparseGRM_cur_SAIGE_org()
          } else {
            set_flagSparseGRM_cur_SAIGE(FALSE)
          }

process_Haplotype_Region(
	traitType,
	genoType,
	n,
	number_of_ancestry,
	region$genoIndex_prev,
	          region$genoIndex
)
if(FALSE){	
	process_Haplotype_Region(
          genoType,
          region$genoIndex_prev,
          region$genoIndex,
          annoIndicatorMat,
          maxMAFlist,
          OutputFile,
          traitType,
          n,
          P1Mat,
          P2Mat,
          regionTestType,
          isImputation,
          WEIGHT,
          weight_cond,
          is_single_in_groupTest,
          is_output_markerList_in_groupTest,
          annolistsub,
          regionName,
          is_fastTest,
          is_output_moreDetails,
          number_of_ancestry
        )
}        
        # Loop through each ancestry for group tests
        for (anc_index in 1:(number_of_ancestry + 1)) {
            anc_index_name <- anc_index
          if (anc_index > number_of_ancestry) {
            anc_index_name <- "ALL"
          }
          set_current_anc_index_name(paste0("DS", anc_index_name));
          # Set conditional analysis flag based on ancestry
          # Individual ancestries: conditional on haplotypes from process_Haplotype_Region
          # ALL ancestry: no conditioning on haplotypes
          isCondition_current <- if (anc_index <= number_of_ancestry) TRUE else FALSE
          
          if (isCondition_current) {
            cat("Processing ancestry", anc_index_name, "for region", regionName, "(conditional on haplotypes)\n")
          } else {
            cat("Processing ancestry", anc_index_name, "for region", regionName, "(unconditional)\n")
          }
          
          # Set sparse GRM flag
          if (!is_fastTest) {
            set_flagSparseGRM_cur_SAIGE_org()
          } else {
            set_flagSparseGRM_cur_SAIGE(FALSE)
          }
          
          # Main analysis
          outList <- mainRegionInCPP(
            genoType, region$genoIndex_prev, region$genoIndex, annoIndicatorMat,
            maxMAFlist, OutputFile, traitType, n, P1Mat, P2Mat, regionTestType,
            isImputation, WEIGHT, weight_cond, is_single_in_groupTest,
            is_output_markerList_in_groupTest, annolistsub, regionName,
            is_fastTest, is_output_moreDetails
          )
          
          # Process results
          if (regionTestType == "BURDEN" & is_fastTest) {
            if (!is.null(outList$iswriteOutput) && !(outList$iswriteOutput)) {
              set_flagSparseGRM_cur_SAIGE(TRUE)
              
              # Re-setup VCF if needed
              vcf_result <- setup_vcf_iterator(genoType, SNP, regionName, chrom1)
              if (genoType == "vcf" && !is.null(vcf_result)) {
                region$genoIndex <- vcf_result$genoIndex
                region$genoIndex_prev <- vcf_result$genoIndex_prev
              }
              
              # Re-run analysis
              outList <- mainRegionInCPP(
                genoType, region$genoIndex_prev, region$genoIndex, annoIndicatorMat,
                maxMAFlist, OutputFile, traitType, n, P1Mat, P2Mat, regionTestType,
                isImputation, WEIGHT, weight_cond, is_single_in_groupTest,
                is_output_markerList_in_groupTest, annolistsub, regionName,
                is_fastTest, is_output_moreDetails
              )
            }
          }
          
          print("outList")
          print(outList)
          
          # Process analysis results
          analysis_result <- process_region_analysis(
            outList, regionTestType, is_fastTest, pval_cutoff_for_fastTest,
            is_single_in_groupTest, traitType, WEIGHT, annolistsub, regionName,
            maxMAFlist, mu, isCondition_current, r.corr, genoType, region, SNP, chrom1,
            annoIndicatorMat, OutputFile, n, P1Mat, P2Mat, isImputation,
            weight_cond, is_output_markerList_in_groupTest, is_output_moreDetails,
            number_of_ancestry, anc_index, anc_index_name
          )
          
          ancestry_outList <- analysis_result$outList
          ancestry_pval.Region <- analysis_result$pval.Region
          
          # Accumulate results from this ancestry
          if (anc_index == 1) {
            outList <- ancestry_outList
            pval.Region <- ancestry_pval.Region
          } else {
            # Combine results from multiple ancestries
            if (!is.null(ancestry_pval.Region)) {
              pval.Region <- rbind(pval.Region, ancestry_pval.Region)
            }
          }
          
          # Clean up ancestry-specific results
          rm(ancestry_outList)
          if (!is.null(ancestry_pval.Region)) {
            rm(ancestry_pval.Region)
          }
          gc()
        }
        
        # Generate marker list
        Output_MarkerList <- generate_marker_list(outList, annolistsub, maxMAFlist, SNP,
                                                 regionName, is_output_markerList_in_groupTest)
        
        # Combine individual ancestry results using CCT to create HET results
        if (regionTestType != "BURDEN" && !is.null(pval.Region) && number_of_ancestry > 1) {
          cat("Combining individual ancestry results using CCT to create HET results for region", regionName, "\n")
          het_results <- combine_ancestry_results_cct(pval.Region, number_of_ancestry, 
                                                     annolistsub, maxMAFlist, regionName)
          
          # Add HET results to the region results
          if (!is.null(het_results)) {
            pval.Region <- rbind(pval.Region, het_results)
          }
          
          # Create FINAL combination of HET and ALL results
          cat("Creating FINAL combination of HET and ALL results for region", regionName, "\n")
          final_results <- combine_het_all_final(pval.Region, annolistsub, maxMAFlist, regionName)
          
          # Add FINAL results to the region results
          if (!is.null(final_results)) {
            pval.Region <- rbind(pval.Region, final_results)
          }
        }
        
        # Create within-ancestry Cauchy combinations across Groups and MAF
        if (regionTestType != "BURDEN" && !is.null(pval.Region) && 
            (length(annolistsub) > 1 || length(maxMAFlist) > 1)) {
          cat("Creating within-ancestry Cauchy combinations for region", regionName, "\n")
          
          # Get all ancestry types present in the results
          ancestries <- unique(pval.Region$Ancestry)
          
          # Create within-ancestry combinations
          within_ancestry_results <- create_within_ancestry_combinations(
            pval.Region, annolistsub, maxMAFlist, regionName, ancestries
          )
          
          # Add within-ancestry combination results
          if (!is.null(within_ancestry_results)) {
            pval.Region <- rbind(pval.Region, within_ancestry_results)
          }
        }
        
        # Cleanup
        if (!is_fastTest) {
          rm(region)
        }
        
        # Add to accumulated results
        indexChunk <- i
        Start <- (cth_chunk_to_output == 1)
        End <- (i == nRegions)
        
        if (regionTestType != "BURDEN") {
          pval.Region.all <- rbind(pval.Region.all, pval.Region)
        }
        
        if (is_output_markerList_in_groupTest) {
          Output_MarkerList.all <- rbind(Output_MarkerList.all, Output_MarkerList)
          rm(Output_MarkerList)
        }
        
        rm(outList)
        if (!is.null(pval.Region)) {
          rm(pval.Region)
        }
        gc()
        
      } else {
        cat(regionName, " is empty.\n")
      }
      
      # Write output when chunk is complete
      if (mth == numberRegionsInChunk) {
        cat("write to output\n")
        
        # Write main output files
        output_result <- write_output_files(OutputFile, pval.Region.all, Output_MarkerList.all,
                                           regionTestType, Start, is_output_markerList_in_groupTest,
                                           iswriteMarkerList)
        iswriteMarkerList <- output_result$iswriteMarkerList
        Output_MarkerList.all <- output_result$Output_MarkerList.all
        
        # Write index file
        write_index_file(OutputFileIndex, Start, End, indexChunk, nEachChunk)
        
        # Reset for next chunk
        pval.Region.all <- NULL
        OutList.all <- NULL
        Output_MarkerList.all <- NULL
        cth_chunk_to_output <- cth_chunk_to_output + 1
        gc()
      }
      
    } else {
      cat("The chunk is empty\n")
      mth <- numberRegionsInChunk
      i <- i + numberRegionsInChunk
      pval.Region <- NULL
    }
  }
  
  # Generate final message
  message <- "Analysis done!"
  message <- paste0(message, " The set-based tests results have been saved to '", OutputFile, "'.")
  
  if (is_output_markerList_in_groupTest) {
    message <- paste0(message, " The marker lists have been saved to '", OutputFile, ".markerList.txt'.")
  }
  
  if (is_single_in_groupTest) {
    message <- paste0(message, " The single-variant association tests results have been saved to '",
                     OutputFile, ".singleAssoc.txt'.")
  }
  
  return(message)
}
