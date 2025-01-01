    // marker-level information
if(g_isadmixed){
    i = q0;
    double altFreq, altCounts, missingRate, imputeInfo;
    std::vector<uint32_t> indexForMissing;
    std::string chr, ref, alt, pds, marker;
    uint32_t pd;
    bool flip = false;
    bool isOutputIndexForMissing = true;
    bool isOnlyOutputNonZero = false;

    //GVec.resize(t_n);
    //GVec.zeros();
    GVec = GVec_sumdosage;
    chr = chrVec.at(i-1);
    pds = posVec.at(i-1);
    ref = refVec.at(i-1);
    alt = "SUM";
    std::string info = chr+":"+pds+":"+ref+":"+alt; 
    marker = info+"_SUM";
    altCounts = arma::accu(GVec);
    altFreq = altCounts/(2 * t_n);
    double MAF = std::min(altFreq, 1 - altFreq);
    double w0;
    double MAC = MAF * 2 * t_n;   

    arma::uvec indexZeroVec_arma, indexNonZeroVec_arma;
    MAF = std::min(altFreq, 1 - altFreq);
    MAC = std::min(altCounts, t_n *2 - altCounts);
    chrVec.at(i) = chr;
    posVec.at(i) = pds;
    refVec.at(i) = ref;
    altVec.at(i) = alt;  	
    markerVec.at(i) = marker;             // marker IDs
    infoVec.at(i) = 1;                 // marker information: CHR:POS:REF:ALT
    altFreqVec.at(i) = altFreq;           // allele frequencies of ALT allele, this is not always < 0.5.
    missingRateVec.at(i) = 0;
    altCountsVec.at(i) = altCounts;
    MACVec.at(i) = MAC;
    MAFVec.at(i) = MAF;
    imputationInfoVec.at(i) = 1;
    arma::uvec indexZeroVec_arma = arma::find(GVec == 0);
    arma::uvec indexNonZeroVec_arma = arma::find(GVec != 0);
      

    if(isWeightCustomized){
        w0 = t_weight(i);
    }else{
        w0 = boost::math::pdf(beta_dist, MAF);
    }


    uint nNonZero = indexNonZeroVec_arma.n_elem;

    if(MAC > g_region_minMAC_cutoff){  // not Ultra-Rare Variants

      indicatorVec.at(i) = 1;
      if(!isSingleVarianceRatio){
        hasVarRatio = ptr_gSAIGEobj->assignVarianceRatio(MAC, ptr_gSAIGEobj->m_flagSparseGRM_cur);
      }else{
        ptr_gSAIGEobj->assignSingleVarianceRatio(ptr_gSAIGEobj->m_flagSparseGRM_cur);
      }


      if(t_regionTestType != "BURDEN" || t_isSingleinGroupTest){ //perform single-variant assoc tests

        //indexZeroVec_arma = arma::conv_to<arma::uvec>::from(indexZeroVec);

        //set_varianceRatio(MAC, isSingleVarianceRatio);
        if(MAC > g_MACCutoffforER){
          Unified_getMarkerPval(
                    GVec,
                    false, // bool t_isOnlyOutputNonZero,
          indexNonZeroVec_arma, indexZeroVec_arma, Beta, seBeta, pval, pval_noSPA, Tstat, gy, varT, altFreq, isSPAConverge, gtildeVec, is_gtilde, true, P2Vec, isCondition, Beta_c, seBeta_c, pval_c, pval_noSPA_c, Tstat_c, varT_c, G1tilde_P_G2tilde_Vec, is_Firth, is_FirthConverge, false);

        }else{
          Unified_getMarkerPval(
                    GVec,
                    false, // bool t_isOnlyOutputNonZero,
          indexNonZeroVec_arma, indexZeroVec_arma, Beta, seBeta, pval, pval_noSPA, Tstat, gy, varT, altFreq, isSPAConverge, gtildeVec, is_gtilde, true, P2Vec, isCondition, Beta_c, seBeta_c, pval_c, pval_noSPA_c, Tstat_c, varT_c, G1tilde_P_G2tilde_Vec, is_Firth, is_FirthConverge, true);
        }



        BetaVec.at(i) = Beta * (1 - 2*flip);  // Beta if flip = false, -1 * Beta is flip = true
        seBetaVec.at(i) = seBeta;
        pvalVec.at(i) = pval;
        //pvalVec_val.at(i) = std::stod(pval);
        pvalNAVec.at(i) = pval_noSPA;
        TstatVec.at(i) = Tstat * (1 - 2*flip);
        TstatVec_flip.at(i) = Tstat;
        gyVec.at(i) = gy;
        varTVec.at(i) = varT;
        isSPAConvergeVec.at(i) = isSPAConverge;

        if(isCondition){
          Beta_cVec.at(i) = Beta_c * (1 - 2*flip);  // Beta if flip = false, -1 * Beta is flip = true
          seBeta_cVec.at(i) = seBeta_c;
          pval_cVec.at(i) = pval_c;
          pvalNA_cVec.at(i) = pval_noSPA_c;
          Tstat_cVec.at(i) = Tstat_c * (1 - 2*flip);
          varT_cVec.at(i) = varT_c;
          G1tilde_P_G2tilde_Weighted_Mat.row(i) = G1tilde_P_G2tilde_Vec % w0G2Vec_cond.t() * w0;
        }


     }//if(t_regionTestType != "BURDEN" || t_isSingleinGroupTest){

     arma::vec dosage_case, dosage_ctrl;
     if(t_traitType == "binary" || t_traitType == "survival"){
        dosage_case = GVec.elem(ptr_gSAIGEobj->m_case_indices);
        dosage_ctrl = GVec.elem(ptr_gSAIGEobj->m_ctrl_indices);
        MACcasegroup = arma::accu(dosage_case);
        MACcontrolgroup = arma::accu(dosage_ctrl);
    }
     if(t_regionTestType != "BURDEN" || t_isSingleinGroupTest){
      if(t_traitType == "binary" || t_traitType == "survival"){
        AF_case = arma::mean(dosage_case) /2;
        AF_ctrl = arma::mean(dosage_ctrl) /2;
        if(flip){
          AF_case = 1-AF_case;
           AF_ctrl = 1-AF_ctrl;
        }
        AF_caseVec.at(i) = AF_case;
        AF_ctrlVec.at(i) = AF_ctrl;
        N_case = dosage_case.n_elem;
        N_ctrl = dosage_ctrl.n_elem;
        N_caseVec.at(i) = N_case;
        N_ctrlVec.at(i) = N_ctrl;


        arma::uvec N_case_ctrl_het_hom0;
        if(t_isMoreOutput){
          N_case_ctrl_het_hom0 = arma::find(dosage_case <= 2 && dosage_case >=1.5);
          N_case_homVec.at(i)  = N_case_ctrl_het_hom0.n_elem;
          N_case_ctrl_het_hom0 = arma::find(dosage_case < 1.5 && dosage_case >= 0.5);
          N_case_hetVec.at(i) = N_case_ctrl_het_hom0.n_elem;
          N_case_ctrl_het_hom0 = arma::find(dosage_ctrl <= 2 && dosage_ctrl >=1.5);
          N_ctrl_homVec.at(i) = N_case_ctrl_het_hom0.n_elem;
          N_case_ctrl_het_hom0 = arma::find(dosage_ctrl < 1.5 && dosage_ctrl >= 0.5);
          N_ctrl_hetVec.at(i) = N_case_ctrl_het_hom0.n_elem;
          if(flip){
                N_case_homVec.at(i) = N_case - N_case_hetVec.at(i) -  N_case_homVec.at(i);
                N_ctrl_homVec.at(i) = N_ctrl - N_ctrl_hetVec.at(i) - N_ctrl_homVec.at(i);
          }
        }
      }else if(t_traitType == "quantitative"){
        N_Vec.at(i) = t_n;
      }
     } //if(t_regionTestType != "BURDEN" || t_isSingleinGroupTest){

      //arma::vec timeoutput4a = getTime();
      //printTime(timeoutput3ab, timeoutput4a, "Unified_getOneMarker 3");

    }//    if(MAC > g_region_minMAC_cutoff){  // not Ultra-Rare Variants
//
}//if(g_isadmixed){	
