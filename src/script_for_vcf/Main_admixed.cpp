//////// ---------- Main function for marker-level analysis in admixed samples--------- ////////////
                           //std::vector<uint32_t> & t_genoIndex,

// [[Rcpp::export]]
void mainMarkerAdmixedInCPP(
                           std::string & t_genoType,     // "PLINK", "BGEN"
                           std::string & t_traitType,
                           std::vector<std::string> & t_genoIndex_prev,
                           std::vector<std::string> & t_genoIndex,
                           bool & t_isMoreOutput,
                           bool & t_isImputation,
                           bool & t_isFirth, 
			   int t_NumberofANC)
{

  int q = t_genoIndex.size();  // number of markers
  // set up output
  std::vector<std::string> markerVec(q);  // marker IDs
  std::vector<std::string> chrVec(q);  // marker IDs
  std::vector<std::string> posVec(q);  // marker IDs
  std::vector<std::string> refVec(q);  // marker IDs
  std::vector<std::string> altVec(q);  // marker IDs
  std::vector<std::string> infoVec(q);    // marker information: CHR:POS:REF:ALT

  q = q * t_NumberofANC + 1; 
  std::vector<double> altFreqVec(q);      // allele frequencies of ALT allele, this is not always < 0.5.
  std::vector<double> altCountsVec(q);    // allele counts of ALT allele.
  std::vector<double> imputationInfoVec(q);    // imputation info of ALT allele.
  std::vector<double> missingRateVec(q);
  std::vector<double> BetaVec(q, arma::datum::nan);         // beta value for ALT allele
  std::vector<double> seBetaVec(q, arma::datum::nan);
  //std::vector<double> pvalVec(q, arma::datum::nan);
  std::vector<std::string> pvalVec(q, "NA");
  std::vector<double> TstatVec(q, arma::datum::nan);
  std::vector<double> varTVec(q, arma::datum::nan);
  //std::vector<double> pvalNAVec(q, arma::datum::nan);
  std::vector<std::string> pvalNAVec(q, "NA");
  bool isCondition = ptr_gSAIGEobj->m_isCondition;
  std::vector<double> Beta_cVec(q, arma::datum::nan);         // beta value for ALT allele
  std::vector<double> seBeta_cVec(q, arma::datum::nan);
  //std::vector<double> pval_cVec(q, arma::datum::nan);
  std::vector<std::string> pval_cVec(q, "NA");
  std::vector<double> Tstat_cVec(q, arma::datum::nan);
  std::vector<double> varT_cVec(q, arma::datum::nan);
  //std::vector<double> pvalNA_cVec(q, arma::datum::nan);
  std::vector<std::string> pvalNA_cVec(q, "NA");
  arma::rowvec G1tilde_P_G2tilde_Vec(ptr_gSAIGEobj->m_numMarker_cond);

  std::vector<bool>  isSPAConvergeVec(q);
  std::vector<double>  AF_caseVec(q);
  std::vector<double>  AF_ctrlVec(q);
  std::vector<uint32_t>  N_caseVec(q);
  std::vector<uint32_t>  N_ctrlVec(q);
  std::vector<double>  N_case_homVec(q);
  std::vector<double>  N_ctrl_hetVec(q);
  std::vector<double>  N_case_hetVec(q);
  std::vector<double>  N_ctrl_homVec(q);
  std::vector<uint32_t>  N_Vec(q);
  
   
  std::vector<uint> indexZeroVec;
  std::vector<uint> indexNonZeroVec;
  std::vector<uint> indexForMissing;


  int n = ptr_gSAIGEobj->m_n;
  arma::vec t_GVec(n);
  arma::vec t_GVecHom(n);
  arma::vec gtildeVec(n);
  arma::vec t_P2Vec;
  if(ptr_gSAIGEobj->m_isFastTest){
    ptr_gSAIGEobj->set_flagSparseGRM_cur(false);
  }else{
    ptr_gSAIGEobj->set_flagSparseGRM_cur(ptr_gSAIGEobj->m_flagSparseGRM);
  }


  bool hasVarRatio = true;;
  bool isSingleVarianceRatio = true;
  if((ptr_gSAIGEobj->m_varRatio_null).n_elem == 1){
        ptr_gSAIGEobj->assignSingleVarianceRatio(ptr_gSAIGEobj->m_flagSparseGRM_cur);
        //ptr_gSAIGEobj->assignSingleVarianceRatio(false);
  }else{
        isSingleVarianceRatio = false;
  }


  int mFirth = 0;
  int mFirthConverge = 0;

  q = t_genoIndex.size();
  
  std::vector<std::string> pvalHet_Vec(q, "NA");  
  std::vector<std::string> pvalHom_Vec(q, "NA");  
  std::vector<std::string> pvalAdmixed_Vec(q, "NA");  
  std::vector<std::string> pvalHet_cVec(q, "NA");  
  std::vector<std::string> pvalHom_cVec(q, "NA");  
  std::vector<std::string> pvalAdmixed_cVec(q, "NA");  

  arma::mat P1Mat = arma::zeros<arma::mat>(t_NumberofANC, n);
  arma::mat P2Mat = arma::zeros<arma::mat>(n, t_NumberofANC);
  arma::vec Scorevec = arma::zeros<arma::vec>(t_NumberofANC);
  arma::vec Scorevec_c = arma::zeros<arma::vec>(t_NumberofANC);
  arma::mat G1tilde_P_G2tilde_Mat(t_NumberofANC, ptr_gSAIGEobj->m_numMarker_cond);


  for(int i = 0; i < q; i++){
    if((i+1) % g_marker_chunksize == 0){
      std::cout << "Completed " << (i+1) << "/" << q << " markers in the chunk." << std::endl;
    }
    t_GVecHom.zeros();
    Scorevec.clear();
    Scorevec_c.clear();
    P1Mat.clear();
    P2Mat.clear();
    G1tilde_P_G2tilde_Mat.clear();

    for(int j = 0; j < t_NumberofANC+1; j++){    
    // information of marker
    double altFreq, altCounts, missingRate, imputeInfo, AF_case, AF_ctrl, N_case_hom, N_ctrl_het, N_case_het, N_ctrl_hom;
    std::string chr, ref, alt, marker;
    uint32_t pd, N_case, N_ctrl, N;
    //free(end);
    bool flip = false;
    std::string t_genoIndex_str = t_genoIndex.at(i);
    char* end;
    uint64_t gIndex = std::strtoull( t_genoIndex_str.c_str(), &end,10 );
    std::remove(end);

    uint64_t gIndex_prev = 0;
    if(i == 0){
        gIndex_prev = 0;
    }else{
        char* end_prev;
        std::string t_genoIndex_prev_str;
        if(t_genoType == "bgen"){
            t_genoIndex_prev_str = t_genoIndex_prev.at(i-1);
        }else if(t_genoType == "plink"){
            t_genoIndex_prev_str = t_genoIndex.at(i-1);
        }
        gIndex_prev = std::strtoull( t_genoIndex_prev_str.c_str(), &end_prev,10 );
        std::remove(end_prev);
    }




    //Main.cpp
    //PLINK or BGEN
    //uint32_t gIndex_temp = gIndex;
    bool isOutputIndexForMissing = true;
    bool isOnlyOutputNonZero = false;

   //clear vectors
   indexZeroVec.clear();
   indexNonZeroVec.clear();
   indexForMissing.clear();
   //t_GVec0.clear();
   //t_GVec.clear();
   //arma::vec timeoutput1 = getTime();
   if(j < t_NumberofANC){ 
   std::string vcfFieldtoRead = "DS" + std::to_string(j+1);
   bool isReadMarker = Unified_getOneMarker_Admixed(t_genoType, gIndex_prev, gIndex, ref, alt, marker, pd, chr, altFreq, altCounts, missingRate, imputeInfo,
                                          isOutputIndexForMissing, // bool t_isOutputIndexForMissing,
                                          indexForMissing,
                                          isOnlyOutputNonZero, // bool t_isOnlyOutputNonZero,
                                          indexNonZeroVec, t_GVec, t_isImputation, vcfFieldtoRead);
   //arma::vec timeoutput2 = getTime();
   //printTime(timeoutput1, timeoutput2, "Unified_getOneMarker");
   t_GVecHom = t_GVecHom + t_GVec;


   if(!isReadMarker){
      //std::cout << "isReadMarker " << isReadMarker << std::endl;
      g_markerTestEnd = true;
      bool isEndFile = check_Vcf_end();
      break;
    }
  }else{//test combined genotypes if(j < t_NumberofANC){
    t_GVec = t_GVecHom;
    altCounts = arma::accu(t_GVec);
    altFreq = arma::mean(altCounts)/2.0;
    missingRate = 0;
    imputeInfo = 1;
  }

   //std::cout << "t_GVec0.size()) " << t_GVec0.size() << std::endl;
   //arma::vec t_GVec(t_GVec0.size());
   //arma::vec t_GVec = arma::conv_to< arma::colvec >::from(t_GVec0);

   //arma::vec t_GVec(t_GVec0);
   //t_GVec0.clear();

   //for(uint j = 0; j < n; j++){
   //   t_GVec(j) = t_GVec0.at(j);
   //}

    //for(int indi = 0; indi < indexForNonZero.size(); indi++){
    //  std::cout << indexForNonZero[indi] << std::endl;
    //}
//   std::cout << "marker " << marker << std::endl;
//   std::cout << "indexForMissing.size() " << indexForMissing.size() << std::endl;
//   std::cout << "indexNonZeroVec.size() " << indexNonZeroVec.size() << std::endl;
    //int n = t_GVec.size();
    //arma::vec gtildeVec(n);




    std::string pds = std::to_string(pd);
    std::string info = chr+":"+pds+":"+ref+":"+alt;

    chrVec.at(i) = chr;
    posVec.at(i) = pds;
    refVec.at(i) = ref;
    altVec.at(i) = alt;
    // record basic information for the marker
    markerVec.at(i) = marker;               // marker IDs
    infoVec.at(i) = info;    // marker information: CHR:POS:REF:ALT
    int k = i*3 + j; 
    altFreqVec.at(k) = altFreq;         // allele frequencies of ALT allele, this is not always < 0.5.
    //altCountsVec.at(i) = altCounts;         // allele frequencies of ALT allele, this is not always < 0.5.
    missingRateVec.at(k) = missingRate;
    imputationInfoVec.at(k) = imputeInfo;
    // MAF and MAC are for Quality Control (QC)
    double MAF = std::min(altFreq, 1 - altFreq);
    double MAC = MAF * n * (1 - missingRate) *2;



  /*
   std::cout << "missingRate " << missingRate << std::endl;
   std::cout << "MAF " << MAF << std::endl;
   std::cout << "MAC " << MAC << std::endl;
   std::cout << "altFreq " << altFreq << std::endl;
   std::cout << "altCounts " << altCounts << std::endl;
   std::cout << "n " << n << std::endl;
   */


    // Quality Control (QC) based on missing rate, MAF, and MAC
    if((missingRate > g_missingRate_cutoff) || (MAF < g_marker_minMAF_cutoff) || (MAC < g_marker_minMAC_cutoff || imputeInfo < g_marker_minINFO_cutoff)){
      continue;
    }else{
    // Check UTIL.cpp
    //
    //
    //arma::vec timeoutput3 = getTime();
    indexZeroVec.clear();
    indexNonZeroVec.clear();


    flip = imputeGenoAndFlip_fakeflip(t_GVec, altFreq, altCounts,indexForMissing, g_impute_method, g_dosage_zerod_cutoff, g_dosage_zerod_MAC_cutoff, MAC, indexZeroVec, indexNonZeroVec);
    MAC = std::min(altCounts, 2*n-altCounts);
    MAF = std::min(altFreq, 1 - altFreq);

   if((MAF < g_marker_minMAF_cutoff) || (MAC < g_marker_minMAC_cutoff)){
        continue;
   }else{
    //arma::vec timeoutput4 = getTime();
    //printTime(timeoutput3, timeoutput4, "imputeGenoAndFlip");
    altFreqVec.at(k) = altFreq;         // allele frequencies of ALT allele, this is not always < 0.5.
    altCountsVec.at(k) = altCounts;         // allele frequencies of ALT allele, this is not always < 0.5.
    //std::cout << "MAC " << MAC << std::endl;
    //std::cout << "altFreq after flip " << altFreq << std::endl;
    //std::cout << "info " << info << std::endl;
    // analysis results for single-marker
    double Beta, seBeta, Tstat, varT, gy;
    double Beta_c, seBeta_c, Tstat_c, varT_c;
    std::string pval, pval_noSPA, pval_c, pval_noSPA_c;
    bool isSPAConverge, is_gtilde, is_Firth, is_FirthConverge;
    //arma::vec t_P2Vec;
    //arma::vec t_P2Vec;
    arma::uvec indexZeroVec_arma, indexNonZeroVec_arma;
    indexZeroVec_arma = arma::conv_to<arma::uvec>::from(indexZeroVec);
    indexNonZeroVec_arma = arma::conv_to<arma::uvec>::from(indexNonZeroVec);

    indexZeroVec.clear();
    indexNonZeroVec.clear();
    t_P2Vec.clear();
    G1tilde_P_G2tilde_Vec.clear();
    //arma::vec timeoutput5 = getTime();
    //set_varianceRatio(MAC, isSingleVarianceRatio);

    if(ptr_gSAIGEobj->m_isFastTest){
      ptr_gSAIGEobj->set_flagSparseGRM_cur(false);

      if(isSingleVarianceRatio){
        ptr_gSAIGEobj->assignSingleVarianceRatio(ptr_gSAIGEobj->m_flagSparseGRM_cur);
      }else{
        hasVarRatio = ptr_gSAIGEobj->assignVarianceRatio(MAC, ptr_gSAIGEobj->m_flagSparseGRM_cur);
      }
    }else{
      if(!isSingleVarianceRatio){
        hasVarRatio = ptr_gSAIGEobj->assignVarianceRatio(MAC, ptr_gSAIGEobj->m_flagSparseGRM_cur);
      }
    }
    //check 'Main.cpp'
    bool is_region = false;

    if(MAC > g_MACCutoffforER){
      Unified_getMarkerPval(
                    t_GVec,
                          false, // bool t_isOnlyOutputNonZero,
                          indexNonZeroVec_arma, indexZeroVec_arma, Beta, seBeta, pval, pval_noSPA,  Tstat, gy, varT,
                          altFreq, isSPAConverge, gtildeVec, is_gtilde, is_region, t_P2Vec, isCondition, Beta_c, seBeta_c, pval_c, pval_noSPA_c, Tstat_c, varT_c, G1tilde_P_G2tilde_Vec, is_Firth, is_FirthConverge, false);
    }else{
      Unified_getMarkerPval(
                    t_GVec,
                          false, // bool t_isOnlyOutputNonZero,
                          indexNonZeroVec_arma, indexZeroVec_arma, Beta, seBeta, pval, pval_noSPA, Tstat, gy, varT,
                          altFreq, isSPAConverge, gtildeVec, is_gtilde, is_region, t_P2Vec, isCondition, Beta_c, seBeta_c, pval_c, pval_noSPA_c, Tstat_c, varT_c, G1tilde_P_G2tilde_Vec, is_Firth, is_FirthConverge, true);
    }

    double pval_num;

    try {
        pval_num = std::stod(pval);
    } catch (const std::invalid_argument&) {
        std::cerr << "Argument is invalid\n";
        pval_num = 0;
    } catch (const std::out_of_range&) {
        std::cerr << "Argument is out of range for a double\n";
        pval_num = 0;
    }

    if(ptr_gSAIGEobj->m_isFastTest && pval_num < (ptr_gSAIGEobj->m_pval_cutoff_for_fastTest)){
      ptr_gSAIGEobj->set_flagSparseGRM_cur(true);

      if(!isSingleVarianceRatio){
        hasVarRatio = ptr_gSAIGEobj->assignVarianceRatio(MAC, ptr_gSAIGEobj->m_flagSparseGRM_cur);
      }else{
        ptr_gSAIGEobj->assignSingleVarianceRatio(ptr_gSAIGEobj->m_flagSparseGRM_cur);
      }


     if(MAC > g_MACCutoffforER){
      Unified_getMarkerPval(
                    t_GVec,
                          false, // bool t_isOnlyOutputNonZero,
                          indexNonZeroVec_arma, indexZeroVec_arma, Beta, seBeta, pval, pval_noSPA, Tstat, gy, varT,
                          altFreq, isSPAConverge, gtildeVec, is_gtilde, is_region, t_P2Vec, isCondition, Beta_c, seBeta_c, pval_c, pval_noSPA_c, Tstat_c, varT_c, G1tilde_P_G2tilde_Vec, is_Firth, is_FirthConverge, false);
     }else{
      Unified_getMarkerPval(
                    t_GVec,
                          false, // bool t_isOnlyOutputNonZero,
                          indexNonZeroVec_arma, indexZeroVec_arma, Beta, seBeta, pval, pval_noSPA, Tstat, gy, varT,
                          altFreq, isSPAConverge, gtildeVec, is_gtilde, is_region, t_P2Vec, isCondition, Beta_c, seBeta_c, pval_c, pval_noSPA_c, Tstat_c, varT_c, G1tilde_P_G2tilde_Vec, is_Firth, is_FirthConverge, true);


     }
     }

  if(j < t_NumberofANC){ 
   P1Mat.row(j) = sqrt(ptr_gSAIGEobj->m_varRatioVal)*gtildeVec.t();
   P2Mat.col(j) = sqrt(ptr_gSAIGEobj->m_varRatioVal)*t_P2Vec;
   Scorevec(j) = Tstat;

   if(isCondition){
	G1tilde_P_G2tilde_Mat.row(i) = G1tilde_P_G2tilde_Vec;
   	Scorevec_c(j) = Tstat_c;
   }

  } 

   if(t_traitType == "binary"){
     if(is_Firth){
       mFirth = mFirth + 1;
       if(is_FirthConverge){
                mFirthConverge = mFirthConverge + 1;
       }
     }
   }
//arma::vec timeoutput6 = getTime();
//printTime(timeoutput5, timeoutput6, "Unified_getMarkerPval");

   indexNonZeroVec_arma.clear();
   indexZeroVec_arma.clear();
   //std::cout << "isSPAConverge " << isSPAConverge << std::endl;
    BetaVec.at(k) = Beta * (1 - 2*flip);  // Beta if flip = false, -1*Beta is flip = true
    seBetaVec.at(k) = seBeta;
    pvalVec.at(k) = pval;
    pvalNAVec.at(k) = pval_noSPA;
    TstatVec.at(k) = Tstat * (1 - 2*flip);
    varTVec.at(k) = varT;

    if(isCondition){
        Beta_cVec.at(k) = Beta_c * (1 - 2*flip);  // Beta if flip = false, -1*Beta is flip = true
        seBeta_cVec.at(k) = seBeta_c;
        pval_cVec.at(k) = pval_c;
        pvalNA_cVec.at(k) = pval_noSPA_c;
        Tstat_cVec.at(k) = Tstat_c * (1 - 2*flip);
        varT_cVec.at(k) = varT_c;
    }

    if(t_traitType == "binary" || t_traitType == "survival"){
        arma::vec dosage_case = t_GVec.elem(ptr_gSAIGEobj->m_case_indices);
        arma::vec dosage_ctrl = t_GVec.elem(ptr_gSAIGEobj->m_ctrl_indices);
        AF_case = arma::mean(dosage_case) /2;
        AF_ctrl = arma::mean(dosage_ctrl) /2;
        N_case = dosage_case.n_elem;
        N_ctrl = dosage_ctrl.n_elem;
        if(flip){
           AF_case = 1-AF_case;
           AF_ctrl = 1-AF_ctrl;
        }
        isSPAConvergeVec.at(k) = isSPAConverge;
        AF_caseVec.at(k) = AF_case;
        AF_ctrlVec.at(k) = AF_ctrl;
        N_caseVec.at(k) = N_case;
        N_ctrlVec.at(k) = N_ctrl;

        arma::uvec N_case_ctrl_het_hom0;
        if(t_isMoreOutput){
          N_case_ctrl_het_hom0 = arma::find(dosage_case <= 2 && dosage_case >=1.5);
          N_case_homVec.at(k)  = N_case_ctrl_het_hom0.n_elem;
          N_case_ctrl_het_hom0 = arma::find(dosage_case < 1.5 && dosage_case >= 0.5);
          N_case_hetVec.at(k) = N_case_ctrl_het_hom0.n_elem;
          N_case_ctrl_het_hom0 = arma::find(dosage_ctrl <= 2 && dosage_ctrl >=1.5);
          N_ctrl_homVec.at(k) = N_case_ctrl_het_hom0.n_elem;
          N_case_ctrl_het_hom0 = arma::find(dosage_ctrl < 1.5 && dosage_ctrl >= 0.5);
          N_ctrl_hetVec.at(k) = N_case_ctrl_het_hom0.n_elem;
          if(flip){
                N_case_homVec.at(k) = N_case - N_case_hetVec.at(k) -  N_case_homVec.at(k);
                N_ctrl_homVec.at(k) = N_ctrl - N_ctrl_hetVec.at(k) - N_ctrl_homVec.at(k);
          }
        }
      }else if(t_traitType == "quantitative"){
        N_Vec.at(k) = n;
      }

     }// if((MAF < g_marker_minMAF_cutoff) || (MAC < g_marker_minMAC_cutoff)){

     } //    if((missingRate > g_missingRate_cutoff) || (MAF < g_marker_minMAF_cutoff) || (MAC < g_marker_minMAC_cutoff || imputeInfo < g_marker_minINFO_cutoff)){
    //arma::vec timeoutput3 = getTime();
    //printTime(timeoutput2, timeoutput3, "after Unified_getOneMarker");
    //t_GVec.clear();
   
    if(j == t_NumberofANC){
   	arma::mat VarMat = P1Mat * P2Mat;
	double P_het_admixed = get_jointScore_pvalue(Scorevec, VarMat);
	double P_hom_admixed = pval;
	arma::vec pvecforcct = {P_het_admixed, P_hom_admixed};
    	double P_cct_admixed = CCT_cpp(pvecforcct);
	std::string P_het_admixed_str = convertDoubletoStringPval(P_het_admixed);
	std::string P_hom_admixed_str = convertDoubletoStringPval(P_hom_admixed);
	std::string P_cct_admixed_str = convertDoubletoStringPval(P_cct_admixed);
	pvalHet_Vec.at(i) = P_het_admixed_str;
	pvalHom_Vec.at(i) = P_hom_admixed_str;
	pvalAdmixed_Vec.at(i) = P_cct_admixed_str;

	if(isCondition){
		

      	    arma::mat AdjCondMat = G1tilde_P_G2tilde_Mat * (ptr_gSAIGEobj->m_VarInvMat_cond);
            arma::mat VarMatAdjCond = AdjCondMat * (G1tilde_P_G2tilde_Mat.t());
            arma::vec TstatAdjCond = AdjCondMat * (ptr_gSAIGEobj->m_Tstat_cond);
	    arma::vec Scorevec_cond = Scorevec - TstatAdjCond;
	    arma::mat VarMat_cond = VarMat - VarMatAdjCond;
	    double P_het_admixed_cond = get_jointScore_pvalue(Scorevec_cond, VarMat_cond);	
	    double P_hom_admixed_cond = pval_c;	
	    arma::vec pvecforcct_cond = {P_het_admixed_cond, P_hom_admixed_cond};
            double P_cct_admixed_cond = CCT_cpp(pvecforcct_cond);
	    std::string P_het_admixed_cond_str = convertDoubletoStringPval(P_het_admixed_cond);
            std::string P_hom_admixed_cond_str = convertDoubletoStringPval(P_hom_admixed_cond);
            std::string P_cct_admixed_cond_str = convertDoubletoStringPval(P_cct_admixed_cond);
            pvalHet_cVec.at(i) = P_het_admixed_cond_str;
            pvalHom_cVec.at(i) = P_hom_admixed_cond_str;
            pvalAdmixed_cVec.at(i) = P_cct_admixed_cond_str; 
	
	}	
    }

    }//for(int j = 0; j < t_NumberofANC; j++){
  }//i

  //output
  writeOutfile_single_admixed_new(t_isMoreOutput,
      t_isImputation,
      isCondition,
      t_isFirth,
  mFirth,
  mFirthConverge,
  t_traitType,
  chrVec,
  posVec,
  markerVec,
  refVec,
  altVec,
  altCountsVec,
  altFreqVec,
  imputationInfoVec,
  missingRateVec,
  BetaVec,
  seBetaVec,
  TstatVec,
  varTVec,
  pvalVec,
  pvalNAVec,
  isSPAConvergeVec,
  Beta_cVec,
  seBeta_cVec,
  Tstat_cVec,
  varT_cVec,
  pval_cVec,
  pvalNA_cVec,
  AF_caseVec,
  AF_ctrlVec,
  N_caseVec,
  N_ctrlVec,
  N_case_homVec,
  N_ctrl_hetVec,
  N_case_hetVec,
  N_ctrl_homVec,
  N_Vec,
  pvalHet_Vec,
  pvalHom_Vec,
  pvalAdmixed_Vec,
  pvalHet_cVec,
  pvalHom_cVec,
  pvalAdmixed_cVec,
  t_NumberofANC);

}



void writeOutfile_single_admixed_new(bool t_isMoreOutput,
                        bool t_isImputation,
                        bool t_isCondition,
                        bool t_isFirth,
                         int mFirth,
                         int mFirthConverge,
                        std::string t_traitType,
                        std::vector<std::string> & chrVec,
                        std::vector<std::string> & posVec,
                        std::vector<std::string> & markerVec,
                        std::vector<std::string> & refVec,
                        std::vector<std::string> & altVec,
                        std::vector<double> & altCountsVec,
                        std::vector<double> & altFreqVec,
                        std::vector<double> & imputationInfoVec,
                        std::vector<double> & missingRateVec,
                        std::vector<double> & BetaVec,
                        std::vector<double> & seBetaVec,
                        std::vector<double> & TstatVec,
                        std::vector<double> & varTVec,
                        std::vector<std::string> & pvalVec,
                        std::vector<std::string> & pvalNAVec,
                        std::vector<bool>  & isSPAConvergeVec,
                        std::vector<double> & Beta_cVec,
                        std::vector<double> & seBeta_cVec,
                        std::vector<double> & Tstat_cVec,
                        std::vector<double> & varT_cVec,
                        std::vector<std::string> & pval_cVec,
                        std::vector<std::string> & pvalNA_cVec,
                        std::vector<double> & AF_caseVec,
                        std::vector<double> & AF_ctrlVec,
                        std::vector<uint32_t> & N_caseVec,
                        std::vector<uint32_t> & N_ctrlVec,
                        std::vector<double>  & N_case_homVec,
                        std::vector<double>  & N_ctrl_hetVec,
                        std::vector<double>  & N_case_hetVec,
                        std::vector<double>  & N_ctrl_homVec,
                        std::vector<uint32_t> & N_Vec,
  std::vector<std::string> & pvalHet_Vec,
  std::vector<std::string> & pvalHom_Vec,
  std::vector<std::string> & pvalAdmixed_Vec,
  std::vector<std::string> & pvalHet_cVec,
  std::vector<std::string> & pvalHom_cVec,
  std::vector<std::string> & pvalAdmixed_cVec,
  int t_NumberofANC
){
  int numtest = 0;
  for(unsigned int j = 0; j < pvalHom_Vec.size(); j++){
        if(pvalHomVec.at(j) != "NA"){
                numtest = numtest + 1;
                OutFile_single << chrVec.at(j);
                OutFile_single << "\t";
                OutFile_single << posVec.at(j);
                OutFile_single << "\t";
                OutFile_single << markerVec.at(j);
                OutFile_single << "\t";
                OutFile_single << refVec.at(j);
                OutFile_single << "\t";
                OutFile_single << altVec.at(j);

	    for(unsigned int i = 0; i < t_NumberofANC+1; i++){
		unsigned int k = j*(t_NumberofANC+1) + i;


                OutFile_single << "\t";
                OutFile_single << altCountsVec.at(k);
                OutFile_single << "\t";
                OutFile_single << altFreqVec.at(k);
                OutFile_single << "\t";

                if(t_isImputation){
                        OutFile_single << imputationInfoVec.at(k);
                        OutFile_single << "\t";
                }else{
                        OutFile_single << missingRateVec.at(k);
                        OutFile_single << "\t";

                }
                OutFile_single << BetaVec.at(k);
                OutFile_single << "\t";
                OutFile_single << seBetaVec.at(k);
                OutFile_single << "\t";
                OutFile_single << TstatVec.at(k);
                OutFile_single << "\t";
                OutFile_single << varTVec.at(k);
                OutFile_single << "\t";
                OutFile_single << pvalVec.at(k);
                OutFile_single << "\t";

                if(t_traitType == "binary" || t_traitType == "survival"){
                        OutFile_single << pvalNAVec.at(k);
                        OutFile_single << "\t";
                        OutFile_single << std::boolalpha << isSPAConvergeVec.at(k);
                        OutFile_single << "\t";
                }
                if(t_isCondition){
                        OutFile_single << Beta_cVec.at(k);
                        OutFile_single << "\t";
                        OutFile_single << seBeta_cVec.at(k);
                        OutFile_single << "\t";
                        OutFile_single << Tstat_cVec.at(k);
                        OutFile_single << "\t";
                        OutFile_single << varT_cVec.at(k);
                        OutFile_single << "\t";
                        OutFile_single << pval_cVec.at(k);
                        OutFile_single << "\t";
                        if(t_traitType == "binary" || t_traitType == "survival"){
                                OutFile_single << pvalNA_cVec.at(k);
                                OutFile_single << "\t";
                        }
                }
                if(t_traitType == "binary" || t_traitType == "survival"){
                        OutFile_single << AF_caseVec.at(k);
                        OutFile_single << "\t";
                        OutFile_single << AF_ctrlVec.at(k);
                        OutFile_single << "\t";
                        OutFile_single << N_caseVec.at(k);
                        OutFile_single << "\t";
                        OutFile_single << N_ctrlVec.at(k);

                        if(t_isMoreOutput){
                                OutFile_single << "\t";
                                OutFile_single << N_case_homVec.at(k);
                                OutFile_single << "\t";
                                OutFile_single << N_case_hetVec.at(k);
                                OutFile_single << "\t";
                                OutFile_single << N_ctrl_homVec.at(k);
                                OutFile_single << "\t";
                                OutFile_single << N_ctrl_hetVec.at(k);
                        }
                        //OutFile_single << "\n";
                }else if(t_traitType == "quantitative"){
                        OutFile_single << N_Vec.at(k);
                        //OutFile_single << "\n";

                }
        }//for(unsigned int i = 0; i < t_NumberofANC+1; i++){
 
	 OutFile_single << "\t";
         OutFile_single << pvalHet_Vec.at(j);
	 OutFile_single << "\t";
         OutFile_single << pvalHom_Vec.at(j);
	 OutFile_single << "\t";
         OutFile_single << pvalAdmixed_Vec.at(j);
	 OutFile_single << "\t";
         OutFile_single << pvalHet_cVec.at(j);
	 OutFile_single << "\t";
         OutFile_single << pvalHom_cVec.at(j);
	 OutFile_single << "\t";
         OutFile_single << pvalAdmixed_cVec.at(j);
	 OutFile_single << "\n";

   }

}

  std::cout << numtest << " markers were tested." << std::endl;
  if(t_traitType == "binary"){
      if(t_isFirth){
        std::cout << "Firth approx was applied to " << mFirth << " markers. " << mFirthConverge << " sucessfully converged." <<std::endl;
      }
   }
}

