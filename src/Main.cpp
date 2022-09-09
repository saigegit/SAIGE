// This includes the main codes to connect C++ and R

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include <vector>
#include <thread>         // std::this_thread::sleep_for
#include <chrono>         // std::chrono::seconds
// std::this_thread::sleep_for (std::chrono::seconds(1));
#include <cstdio>         // std::remove
#include <fstream>
#include <string.h>
// Currently, omp does not work well, will check it later
// error: SET_VECTOR_ELT() can only be applied to a 'list', not a 'character'
// remove all Rcpp::List to check if it works
// #include <omp.h>
// // [[Rcpp::plugins(openmp)]]]

#include "Main.hpp"
#include "PLINK.hpp"
#include "BGEN.hpp"
#include "VCF.hpp"
#include "SAIGE_test.hpp"
#include "UTIL.hpp"
#include "CCT.hpp"

#include <Rcpp.h>
#include "getMem.hpp"

#include <boost/math/distributions/beta.hpp>

// global objects for different genotype formats

static PLINK::PlinkClass* ptr_gPLINKobj = NULL;
static BGEN::BgenClass* ptr_gBGENobj = NULL;
static VCF::VcfClass* ptr_gVCFobj = NULL;
// global objects for different analysis methods
static SAIGE::SAIGEClass* ptr_gSAIGEobj = NULL;
//single, SAIGE
//Region, SAIGE-GENE+


// global variables for analysis
std::string g_impute_method;      // "mean", "minor", or "drop", //drop is not allowed
double g_missingRate_cutoff;
//unsigned int g_omp_num_threads;
double g_marker_minMAF_cutoff;
double g_marker_minMAC_cutoff;
double g_region_minMAC_cutoff;    // for Rare Variants (RVs) whose MAC < this value, we aggregate these variants like SAIGE-GENE+ 
double g_marker_minINFO_cutoff;
arma::vec g_region_maxMAF_cutoff;
double g_min_gourpmac_for_burdenonly;
double g_maxMAFLimit;
unsigned int g_region_maxMarkers_cutoff;   // maximal number of markers in one chunk, only used for region-based analysis to reduce memory usage
bool g_isOutputMoreDetails;
int g_marker_chunksize;

std::string g_method_to_CollapseUltraRare;
double g_DosageCutoff_for_UltraRarePresence;

double g_dosage_zerod_MAC_cutoff;
double g_dosage_zerod_cutoff;
bool g_markerTestEnd = false;
arma::vec g_weights_beta(2);

bool  g_is_Firth_beta;
double g_pCutoffforFirth;
double g_MACCutoffforER;


std::ofstream OutFile;
std::ofstream OutFile_singleInGroup;
std::ofstream OutFile_single;
std::ofstream OutFile_singleInGroup_temp;


std::string g_outputFilePrefixGroup;
std::string g_outputFilePrefixSingleInGroup;
std::string g_outputFilePrefixSingleInGroup_temp;
std::string g_outputFilePrefixSingle;





// [[Rcpp::export]]
void setAssocTest_GlobalVarsInCPP(std::string t_impute_method,
		double t_missing_cutoff,
                               double t_min_maf_marker,
                               double t_min_mac_marker,
                               double t_min_info_marker,
			       double t_dosage_zerod_cutoff,
                               double t_dosage_zerod_MAC_cutoff,
			       arma::vec & t_weights_beta, 
			       std::string t_outputFilePrefix,
			       double t_MACCutoffforER)
{
  g_impute_method = t_impute_method;
  g_missingRate_cutoff = t_missing_cutoff;
  g_marker_minMAF_cutoff = t_min_maf_marker;
  g_marker_minMAC_cutoff = t_min_mac_marker;
  g_marker_minINFO_cutoff = t_min_info_marker;
  g_dosage_zerod_cutoff = t_dosage_zerod_cutoff;
  g_dosage_zerod_MAC_cutoff = t_dosage_zerod_MAC_cutoff;
  g_weights_beta = t_weights_beta;
  g_outputFilePrefixGroup = t_outputFilePrefix;
  g_outputFilePrefixSingleInGroup = t_outputFilePrefix + ".singleAssoc.txt";
  g_outputFilePrefixSingleInGroup_temp = t_outputFilePrefix + ".singleAssoc.txt_temp";
  g_outputFilePrefixSingle = t_outputFilePrefix;
  g_MACCutoffforER = t_MACCutoffforER;
}
// [[Rcpp::export]]
void setMarker_GlobalVarsInCPP(
			       bool t_isOutputMoreDetails,
			       int t_marker_chunksize
			       )

{
  g_isOutputMoreDetails = t_isOutputMoreDetails;
  g_marker_chunksize = t_marker_chunksize;
}


//double t_DosageCutoff_for_UltraRarePresence,
			       //std::string t_method_to_CollapseUltraRare,

  //g_method_to_CollapseUltraRare = t_method_to_CollapseUltraRare;
  //g_DosageCutoff_for_UltraRarePresence = t_DosageCutoff_for_UltraRarePresence;
// [[Rcpp::export]]
void setRegion_GlobalVarsInCPP(
                               arma::vec t_max_maf_region,
                               unsigned int t_max_markers_region,
			       double t_MACCutoff_to_CollapseUltraRare, 
			       double t_min_gourpmac_for_burdenonly)
{
  g_region_maxMAF_cutoff = t_max_maf_region;
  g_maxMAFLimit = g_region_maxMAF_cutoff.max();
  g_region_maxMarkers_cutoff = t_max_markers_region;
  g_region_minMAC_cutoff = t_MACCutoff_to_CollapseUltraRare;
  g_min_gourpmac_for_burdenonly = t_min_gourpmac_for_burdenonly;
}



//////// ---------- Main function for marker-level analysis --------- ////////////
                           //std::vector<uint32_t> & t_genoIndex,

// [[Rcpp::export]]
void mainMarkerInCPP(
                           std::string & t_genoType,     // "PLINK", "BGEN"
			   std::string & t_traitType,
			   std::vector<std::string> & t_genoIndex_prev,
			   std::vector<std::string> & t_genoIndex,
			   bool & t_isMoreOutput,
			   bool & t_isImputation,
			   bool & t_isFirth)
{

  int q = t_genoIndex.size();  // number of markers
  // set up output
  std::vector<std::string> markerVec(q);  // marker IDs
  std::vector<std::string> chrVec(q);  // marker IDs
  std::vector<std::string> posVec(q);  // marker IDs
  std::vector<std::string> refVec(q);  // marker IDs
  std::vector<std::string> altVec(q);  // marker IDs

  std::vector<std::string> infoVec(q);    // marker information: CHR:POS:REF:ALT
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
  //if(isCondition){
  std::vector<double> Beta_cVec(q, arma::datum::nan);         // beta value for ALT allele
  std::vector<double> seBeta_cVec(q, arma::datum::nan);
  //std::vector<double> pval_cVec(q, arma::datum::nan);
  std::vector<std::string> pval_cVec(q, "NA");
  std::vector<double> Tstat_cVec(q, arma::datum::nan);
  std::vector<double> varT_cVec(q, arma::datum::nan);
  //std::vector<double> pvalNA_cVec(q, arma::datum::nan);
  std::vector<std::string> pvalNA_cVec(q, "NA");
  //}
  arma::rowvec G1tilde_P_G2tilde_Vec(ptr_gSAIGEobj->m_numMarker_cond);

  std::vector<bool>  isSPAConvergeVec(q);
  std::vector<double>  AF_caseVec(q);
  std::vector<double>  AF_ctrlVec(q);
  std::vector<uint32_t>  N_caseVec(q);
  std::vector<uint32_t>  N_ctrlVec(q);
    //if(t_isMoreOutput){
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

  for(int i = 0; i < q; i++){
    if((i+1) % g_marker_chunksize == 0){
      std::cout << "Completed " << (i+1) << "/" << q << " markers in the chunk." << std::endl;
    }
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
   //
   bool isReadMarker = Unified_getOneMarker(t_genoType, gIndex_prev, gIndex, ref, alt, marker, pd, chr, altFreq, altCounts, missingRate, imputeInfo,
                                          isOutputIndexForMissing, // bool t_isOutputIndexForMissing,
                                          indexForMissing,
                                          isOnlyOutputNonZero, // bool t_isOnlyOutputNonZero,
                                          indexNonZeroVec, t_GVec, t_isImputation);
   //arma::vec timeoutput2 = getTime();   
   //printTime(timeoutput1, timeoutput2, "Unified_getOneMarker"); 
//
   if(!isReadMarker){
      //std::cout << "isReadMarker " << isReadMarker << std::endl;
      g_markerTestEnd = true;
      bool isEndFile = check_Vcf_end();
      break;
    }


   //std::cout << "t_GVec0.size()) " << t_GVec0.size() << std::endl;
   //arma::vec t_GVec(t_GVec0.size());
   //arma::vec t_GVec = arma::conv_to< arma::colvec >::from(t_GVec0);

   //arma::vec t_GVec(t_GVec0);
   //t_GVec0.clear(); 

   //for(uint j = 0; j < n; j++){
   //	t_GVec(j) = t_GVec0.at(j);	
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
    altFreqVec.at(i) = altFreq;         // allele frequencies of ALT allele, this is not always < 0.5.
    //altCountsVec.at(i) = altCounts;         // allele frequencies of ALT allele, this is not always < 0.5.
    missingRateVec.at(i) = missingRate;
    imputationInfoVec.at(i) = imputeInfo;



    // MAF and MAC are for Quality Control (QC)
    double MAF = std::min(altFreq, 1 - altFreq);
    double MAC = MAF * n * (1 - missingRate) *2;
    
    
   /*
   std::cout << "missingRate " << missingRate << std::endl;
   std::cout << "MAF " << MAF << std::endl;
   std::cout << "MAC " << MAC << std::endl;
   std::cout << "altFreq " << altFreq << std::endl;
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


    flip = imputeGenoAndFlip(t_GVec, altFreq, altCounts,indexForMissing, g_impute_method, g_dosage_zerod_cutoff, g_dosage_zerod_MAC_cutoff, MAC, indexZeroVec, indexNonZeroVec);
   
//arma::vec timeoutput4 = getTime();
//printTime(timeoutput3, timeoutput4, "imputeGenoAndFlip");
    altFreqVec.at(i) = altFreq;         // allele frequencies of ALT allele, this is not always < 0.5.
    altCountsVec.at(i) = altCounts;         // allele frequencies of ALT allele, this is not always < 0.5.
    MAC = std::min(altCounts, 2*n-altCounts);
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
    BetaVec.at(i) = Beta * (1 - 2*flip);  // Beta if flip = false, -1*Beta is flip = true       
    seBetaVec.at(i) = seBeta;       
    pvalVec.at(i) = pval;
    pvalNAVec.at(i) = pval_noSPA;
    TstatVec.at(i) = Tstat * (1 - 2*flip);
    varTVec.at(i) = varT;

    if(isCondition){
    	Beta_cVec.at(i) = Beta_c * (1 - 2*flip);  // Beta if flip = false, -1*Beta is flip = true
    	seBeta_cVec.at(i) = seBeta_c;
    	pval_cVec.at(i) = pval_c;
    	pvalNA_cVec.at(i) = pval_noSPA_c;
    	Tstat_cVec.at(i) = Tstat_c * (1 - 2*flip);
    	varT_cVec.at(i) = varT_c;
    }
	
    if(t_traitType == "binary"){ 
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
      isSPAConvergeVec.at(i) = isSPAConverge;
      AF_caseVec.at(i) = AF_case;
      AF_ctrlVec.at(i) = AF_ctrl;

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
      N_Vec.at(i) = n;

    }

    
   } //    if((missingRate > g_missingRate_cutoff) || (MAF < g_marker_minMAF_cutoff) || (MAC < g_marker_minMAC_cutoff || imputeInfo < g_marker_minINFO_cutoff)){
 
  
  
    //t_GVec.clear();
  }

  //output
  writeOutfile_single(t_isMoreOutput,
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
  N_Vec);

}




// a unified function to get single marker from genotype file
bool Unified_getOneMarker(std::string & t_genoType,   // "PLINK", "BGEN", "Vcf"
                               uint64_t & t_gIndex_prev,        // different meanings for different genoType
                               uint64_t & t_gIndex,        // different meanings for different genoType
                               std::string& t_ref,       // REF allele
                               std::string& t_alt,       // ALT allele (should probably be minor allele, otherwise, computation time will increase)
                               std::string& t_marker,    // marker ID extracted from genotype file
                               uint32_t& t_pd,           // base position
                               std::string& t_chr,       // chromosome
                               double& t_altFreq,        // frequency of ALT allele
                               double& t_altCounts,      // counts of ALT allele
                               double& t_missingRate,    // missing rate
                               double& t_imputeInfo,     // imputation information score, i.e., R2 (all 1 for PLINK)
                               bool & t_isOutputIndexForMissing,               // if true, output index of missing genotype data
                               std::vector<uint>& t_indexForMissing,     // index of missing genotype data
                               bool & t_isOnlyOutputNonZero,                   // if true, only output a vector of non-zero genotype. (NOTE: if ALT allele is not minor allele, this might take much computation time)
                               std::vector<uint>& t_indexForNonZero, //
			       arma::vec & t_GVec,
			       bool t_isImputation 
			       )     // the index of non-zero genotype in the all subjects. Only valid if t_isOnlyOutputNonZero == true.
{
  //arma::vec GVec(ptr_gSAIGEobj->m_n);
  bool isBoolRead = true;
  if(t_genoType == "plink"){
   bool isTrueGenotype = true;
   //t_gIndex_prev is after reading the last marker

   //arma::vec timeoutput1 = getTime();
   ptr_gPLINKobj->getOneMarker(t_gIndex_prev, t_gIndex, t_ref, t_alt, t_marker, t_pd, t_chr, t_altFreq, t_altCounts, t_missingRate, t_imputeInfo, 
                                       t_isOutputIndexForMissing, t_indexForMissing, t_isOnlyOutputNonZero, t_indexForNonZero,
                                       isTrueGenotype, t_GVec);   // t_isTrueGenotype, only used for PLINK format.
   //arma::vec timeoutput2 = getTime();
     //   printTime(timeoutput1, timeoutput2, "Unified_getOneMarker a");
  }
  
  if(t_genoType == "bgen"){
    //bool isBoolRead = true;
    ptr_gBGENobj->getOneMarker(t_gIndex_prev, t_gIndex, t_ref, t_alt, t_marker, t_pd, t_chr, t_altFreq, t_altCounts, t_missingRate, t_imputeInfo, 
                                      t_isOutputIndexForMissing, t_indexForMissing, t_isOnlyOutputNonZero, t_indexForNonZero,
                                      isBoolRead, t_GVec, t_isImputation);
  }

  if(t_genoType == "vcf"){
    ptr_gVCFobj->getOneMarker(t_ref, t_alt, t_marker, t_pd, t_chr, t_altFreq, t_altCounts, t_missingRate, t_imputeInfo,
                                      t_isOutputIndexForMissing, t_indexForMissing, t_isOnlyOutputNonZero, t_indexForNonZero, isBoolRead, t_GVec, t_isImputation);
    ptr_gVCFobj->move_forward_iterator(1);
  }	  
  
  return isBoolRead;
}


// [[Rcpp::export]]
uint32_t Unified_getSampleSizeinGeno(std::string & t_genoType){
    uint32_t N0;
    if(t_genoType == "plink"){
	N0 = ptr_gPLINKobj->getN0();
    }
    if(t_genoType == "bgen"){
	N0 = ptr_gBGENobj->getN0();
    }	    

    if(t_genoType == "vcf"){
       N0 = ptr_gVCFobj->getN0();
    }
    return(N0);
}	

// [[Rcpp::export]]
uint32_t Unified_getSampleSizeinAnalysis(std::string & t_genoType){
    uint32_t N;
    if(t_genoType == "plink"){
        N = ptr_gPLINKobj->getN();
    }
    if(t_genoType == "bgen"){
        N = ptr_gBGENobj->getN();
    }

    if(t_genoType == "vcf"){
       N = ptr_gVCFobj->getN();
    }
    return(N);
}




// a unified function to get marker-level p-value
void Unified_getMarkerPval(
                           arma::vec & t_GVec,
                           bool t_isOnlyOutputNonZero,
			   arma::uvec & t_indexForNonZero_vec,
			   arma::uvec & t_indexForZero_vec,
                           double& t_Beta, 
                           double& t_seBeta, 
                           std::string& t_pval,
			   std::string& t_pval_noSPA,
                           double& t_Tstat,
			   double& t_gy,
                           double& t_varT,
                           double t_altFreq, 
			   bool & t_isSPAConverge,
			   arma::vec & t_gtilde,
			   bool & is_gtilde,
			   bool  is_region,
                           arma::vec & t_P2Vec,
			   bool t_isCondition,
                           double& t_Beta_c, 
                           double& t_seBeta_c, 
                           std::string& t_pval_c,
			   std::string& t_pval_noSPA_c,
                           double& t_Tstat_c,
                           double& t_varT_c,
			   arma::rowvec & t_G1tilde_P_G2tilde_Vec, 
			   bool & t_isFirth,
			   bool & t_isFirthConverge, 
			   bool t_isER) 
{
    if(t_isOnlyOutputNonZero == true)
      Rcpp::stop("When using SAIGE method to calculate marker-level p-values, 't_isOnlyOutputNonZero' should be false.");   


    ptr_gSAIGEobj->getMarkerPval(t_GVec, t_indexForNonZero_vec, t_indexForZero_vec, t_Beta, t_seBeta, t_pval, t_pval_noSPA, t_altFreq, t_Tstat, t_gy, t_varT, t_isSPAConverge, t_gtilde, is_gtilde, is_region, t_P2Vec, t_isCondition, t_Beta_c, t_seBeta_c, t_pval_c, t_pval_noSPA_c, t_Tstat_c, t_varT_c, t_G1tilde_P_G2tilde_Vec, t_isFirth, t_isFirthConverge, t_isER); //SAIGE_new.cpp
    
    //t_indexForNonZero_vec.clear();
  
}

//////// ---------- Main functions to set objects for different genotype format --------- ////////////

// [[Rcpp::export]]
void setPLINKobjInCPP(std::string t_bimFile,
                      std::string t_famFile,
                      std::string t_bedFile,
                      std::vector<std::string> & t_SampleInModel,
                      std::string t_AlleleOrder)
{
  ptr_gPLINKobj = new PLINK::PlinkClass(t_bimFile,
                                        t_famFile,
                                        t_bedFile,
                                        //t_SampleInModel,
                                        t_AlleleOrder);
  ptr_gPLINKobj->setPosSampleInPlink(t_SampleInModel);
  
}

// [[Rcpp::export]]
void setBGENobjInCPP(std::string t_bgenFileName,
                     std::string t_bgenFileIndex,
                     std::vector<std::string> & t_SampleInBgen,
                     std::vector<std::string> & t_SampleInModel,
                     std::string t_AlleleOrder)
{
  std::cout << "t_SampleInBgen " << t_SampleInBgen.size() << std::endl;
  ptr_gBGENobj = new BGEN::BgenClass(t_bgenFileName,
                                     t_bgenFileIndex,
                                     t_SampleInBgen,
                                     t_SampleInModel,
				     false,
				     false,
                                     t_AlleleOrder);
  //int n = ptr_gBGENobj->getN();
}


// [[Rcpp::export]]
void setVCFobjInCPP(std::string t_vcfFileName,
            std::string t_vcfFileIndex,
            std::string t_vcfField,
            std::vector<std::string> & t_SampleInModel)
{
  ptr_gVCFobj = new VCF::VcfClass(t_vcfFileName,
		  		t_vcfFileIndex,
				t_vcfField,
				false,
				t_SampleInModel);

bool isEnd = ptr_gVCFobj->check_iterator_end();

}



//////// ---------- Main functions to set objects for different analysis methods --------- ////////////

// [[Rcpp::export]]
void setSAIGEobjInCPP(arma::mat & t_XVX,
        arma::mat & t_XXVX_inv,
        arma::mat & t_XV,
        arma::mat & t_XVX_inv_XV,
	arma::mat & t_Sigma_iXXSigma_iX,
        arma::mat & t_X,
        arma::vec &  t_S_a,
        arma::vec & t_res,
        arma::vec & t_mu2,
        arma::vec & t_mu,
        arma::vec & t_varRatio_sparse,
        arma::vec & t_varRatio_null,
	arma::vec & t_cateVarRatioMinMACVecExclude,
        arma::vec & t_cateVarRatioMaxMACVecInclude,
        double t_SPA_Cutoff,
        arma::vec & t_tauvec,
        std::string t_traitType,
        arma::vec & t_y,
	std::string t_impute_method,
	bool t_flagSparseGRM,
	bool t_isFastTest,
	double t_pval_cutoff_for_fastTest,
        arma::umat & t_locationMat,
        arma::vec & t_valueVec,
        int t_dimNum,
	bool t_isCondition,
	std::vector<uint32_t> & t_condition_genoIndex,
	bool t_is_Firth_beta,
	double t_pCutoffforFirth,
	arma::vec & t_offset,
	arma::vec & t_resout)
{
  // check SAIGE.cpp
  ptr_gSAIGEobj = new SAIGE::SAIGEClass(
	t_XVX,
        t_XXVX_inv,
        t_XV,
        t_XVX_inv_XV,
	t_Sigma_iXXSigma_iX,
        t_X,
        t_S_a,
        t_res,
        t_mu2,
        t_mu,
	t_varRatio_sparse,
        t_varRatio_null,
	t_cateVarRatioMinMACVecExclude,
	t_cateVarRatioMaxMACVecInclude,
        t_SPA_Cutoff,
        t_tauvec,
        t_traitType,
        t_y,
	t_impute_method,
	t_flagSparseGRM,
	t_isFastTest,
	t_pval_cutoff_for_fastTest,
	t_locationMat,
	t_valueVec,
	t_dimNum, 
	t_isCondition,
	t_condition_genoIndex,
	t_is_Firth_beta,
        t_pCutoffforFirth,
	t_offset, 
	t_resout);
  //ptr_gSAIGEobj->m_flagSparseGRM = false;
}



// [[Rcpp::export]]
void setSparseSigmaInCPP(int r, arma::umat & t_locationMatinR, arma::vec & t_valueVecinR)
{
  ptr_gSAIGEobj->setupSparseMat(r, t_locationMatinR, t_valueVecinR);
  ptr_gSAIGEobj->m_flagSparseGRM = true;
}


// [[Rcpp::export]]
Rcpp::List RegionSetUpConditional_binary_InCPP(arma::vec & t_weight_cond){

	unsigned int q_cond = (ptr_gSAIGEobj->m_VarInvMat_cond).n_rows;
  	boost::math::beta_distribution<> beta_dist(g_weights_beta[0], g_weights_beta[1]);
  	arma::vec w0G2Vec_cond(q_cond);
  	double w0G2_cond, MAFG2_cond;
        for(unsigned int ci = 0; ci < q_cond; ci++){
		if(!(t_weight_cond.is_zero())){
			w0G2_cond = t_weight_cond(ci);
		}else{
                	MAFG2_cond = (ptr_gSAIGEobj->m_MAF_cond)[ci];
                	w0G2_cond = boost::math::pdf(beta_dist, MAFG2_cond);
		}
                w0G2Vec_cond.at(ci) = w0G2_cond;
        }
	arma::mat m_VarMat_weighted_cond = (w0G2Vec_cond * (w0G2Vec_cond.t())) % (ptr_gSAIGEobj->m_VarMat_cond);

	 Rcpp::List OutList = Rcpp::List::create(Rcpp::Named("VarMat_G2_cond") = m_VarMat_weighted_cond,
                                          Rcpp::Named("Score_G2_cond") = ptr_gSAIGEobj->m_Tstat_cond,
                                          Rcpp::Named("pval_G2_cond") = ptr_gSAIGEobj->m_p_cond,
                                          Rcpp::Named("gsum_G2_cond") = ptr_gSAIGEobj->m_gsum_cond,
                                          Rcpp::Named("qsum_G2_cond") = ptr_gSAIGEobj->m_qsum_cond
					  );
	 return(OutList);

}


//////// ---------- Main function for region-level analysis --------- ////////////
// [[Rcpp::export]]
Rcpp::List mainRegionInCPP(
                           std::string t_genoType,     // "PLINK", "BGEN"
                           std::vector<std::string> & t_genoIndex_prev,
                           std::vector<std::string> & t_genoIndex,
			   arma::mat & annoIndicatorMat,
			   arma::vec & maxMAFVec, 
                           std::string t_outputFile,
			   std::string t_traitType,
                           unsigned int t_n,           // sample size  
                           arma::mat & P1Mat,            // edited on 2021-08-19: to avoid repeated memory allocation of P1Mat and P2Mat
                           arma::mat & P2Mat, 
			   std::string t_regionTestType, 
			   bool t_isImputation,
			   arma::vec & t_weight,
			   arma::vec & t_weight_cond,
			   bool t_isSingleinGroupTest,
			   bool t_isOutputMarkerList, 
			   std::vector<std::string> & annoStringVec,
			   std::string regionName, 
			   bool t_isFastTest, 
			   bool t_isMoreOutput) 
{


  //create the output list
  Rcpp::List OutList = Rcpp::List::create();
  //arma::vec timeoutput1 = getTime();
  bool isWeightCustomized = false;
  unsigned int q0 = t_genoIndex.size();                 // number of markers (before QC) in one region
  if(!(t_weight.is_zero()) && t_weight.n_elem == q0){
     isWeightCustomized = true;	
  } 	  
  unsigned int q_anno = annoIndicatorMat.n_cols;
  unsigned int q_maf = maxMAFVec.n_elem;
  unsigned int q_anno_maf = q_anno*q_maf;
  arma::mat genoURMat(t_n, q_anno_maf, arma::fill::zeros);
  unsigned int q = q0 + q_anno_maf;
  arma::imat annoMAFIndicatorMat(q, q_anno_maf, arma::fill::zeros);
  arma::ivec annoMAFIndicatorVec(q_anno_maf);
  annoMAFIndicatorVec.zeros();
  arma::vec maxMAFperAnno(q_anno, arma::fill::zeros);
  arma::vec MAFIndicatorVec(maxMAFVec.n_elem);
  MAFIndicatorVec.zeros();


  //setting up conditional markers  
  unsigned int q_cond = (ptr_gSAIGEobj->m_VarInvMat_cond).n_rows;	
  arma::rowvec G1tilde_P_G2tilde_Vec(q_cond);
  boost::math::beta_distribution<> beta_dist(g_weights_beta[0], g_weights_beta[1]);

  bool isCondition = ptr_gSAIGEobj->m_isCondition;
  arma::vec w0G2Vec_cond(q_cond);
  double w0G2_cond, MAFG2_cond;
  if(isCondition){
	for(unsigned int ci = 0; ci < q_cond; ci++){
		if(!(t_weight_cond.is_zero())){
			w0G2_cond = t_weight_cond(ci);
		}else{	
			MAFG2_cond = (ptr_gSAIGEobj->m_MAF_cond)[ci];	
  			w0G2_cond = boost::math::pdf(beta_dist, MAFG2_cond);
		}	
		w0G2Vec_cond.at(ci) = w0G2_cond;
	}
  }
  arma::mat w0G2Mat_cond(q_cond, q_cond);
  w0G2Mat_cond = w0G2Vec_cond * (w0G2Vec_cond.t());
  arma::mat genoSumMat(t_n, q_anno_maf, arma::fill::zeros); //for Phi_cc for binary traits and BURDEN test
  arma::vec genoSumcount_noweight(q_anno_maf, arma::fill::zeros);
  //arma::sp_mat genoSumMat_sp(t_n, q_anno_maf); //for Phi_cc for binary traits and BURDEN test
  std::vector<double> Beta_cVec(q, arma::datum::nan);         // beta value for ALT allele
  std::vector<double> seBeta_cVec(q, arma::datum::nan);
  //std::vector<double> pval_cVec(q, arma::datum::nan);
  std::vector<std::string> pval_cVec(q, "NA");
  std::vector<double> Tstat_cVec(q, arma::datum::nan);
  std::vector<double> varT_cVec(q, arma::datum::nan);
  //std::vector<double> pvalNA_cVec(q, arma::datum::nan);
  std::vector<std::string> pvalNA_cVec(q, "NA");
  arma::mat G1tilde_P_G2tilde_Weighted_Mat(q, q_cond);
  //group test output
  //arma::vec MAC_GroupVec = arma::zeros<vec>(q_anno_maf);
  arma::vec MAC_GroupVec(q_anno_maf);
  MAC_GroupVec.zeros();
  arma::vec MACCase_GroupVec(q_anno_maf);
  MACCase_GroupVec.zeros();
  arma::vec MACControl_GroupVec(q_anno_maf);
  MACControl_GroupVec.zeros();
  arma::vec NumRare_GroupVec(q_anno_maf);
  NumRare_GroupVec.zeros();
  arma::vec NumUltraRare_GroupVec(q_anno_maf);
  NumUltraRare_GroupVec.zeros();
  arma::vec gtildeVec;
  double MACgroup, MACcasegroup, MACcontrolgroup, AF_case, AF_ctrl, N_case_hom, N_ctrl_het, N_case_het, N_ctrl_hom;
  uint32_t N_case, N_ctrl, N; 

  //single-variant assoc output
  arma::uvec indicatorVec(q, arma::fill::zeros);       // 0: does not pass QC, 1: non-URV, 2: URV
  std::vector<std::string> markerVec(q);
  std::vector<std::string> chrVec(q);  // marker IDs
  std::vector<std::string> posVec(q);  // marker IDs
  std::vector<std::string> refVec(q);  // marker IDs
  std::vector<std::string> altVec(q);  // marker IDs
  std::vector<std::string> infoVec(q);    // marker information: CHR:POS:REF:ALT
  std::vector<double> altFreqVec(q, arma::datum::nan);      // allele frequencies of ALT allele, this is not always < 0.5.
  std::vector<double> MACVec(q, arma::datum::nan);      // allele frequencies of ALT allele, this is not always < 0.5.
  std::vector<double> MAFVec(q, arma::datum::nan);      // allele frequencies of ALT allele, this is not always < 0.5.
  std::vector<double> AF_caseVec(q, arma::datum::nan);      // allele frequencies of ALT allele, this is not always < 0.5.
  std::vector<double> AF_ctrlVec(q, arma::datum::nan);      // allele frequencies of ALT allele, this is not always < 0.5.
  std::vector<uint32_t> N_caseVec(q, arma::datum::nan);      // allele frequencies of ALT allele, this is not always < 0.5.
  std::vector<uint32_t> N_ctrlVec(q, arma::datum::nan);      // allele frequencies of ALT allele, this is not always < 0.5.
  std::vector<double>  N_case_homVec(q, arma::datum::nan);
  std::vector<double>  N_ctrl_hetVec(q, arma::datum::nan);
  std::vector<double>  N_case_hetVec(q, arma::datum::nan);
  std::vector<double>  N_ctrl_homVec(q, arma::datum::nan);
  std::vector<uint32_t> N_Vec(q, arma::datum::nan);      // allele frequencies of ALT allele, this is not always < 0.5.
  std::vector<double> altCountsVec(q, arma::datum::nan);    // allele counts of ALT allele.
  std::vector<double> imputationInfoVec(q, arma::datum::nan);    // imputation info of ALT allele.
  std::vector<double> missingRateVec(q, arma::datum::nan);
  std::vector<double> BetaVec(q, arma::datum::nan);         // beta value for ALT allele
  std::vector<double> seBetaVec(q, arma::datum::nan);
  std::vector<std::string> pvalVec(q, "NA");
  //std::vector<double> pvalVec_val(q, arma::datum::nan);
  std::vector<double> TstatVec(q, arma::datum::nan);
  std::vector<double> TstatVec_flip(q, arma::datum::nan);
  std::vector<double> gyVec(q, arma::datum::nan);
  std::vector<double> varTVec(q, arma::datum::nan);
  std::vector<std::string> pvalNAVec(q, "NA");  
  std::vector<bool>  isSPAConvergeVec(q);


  unsigned int m1 = g_region_maxMarkers_cutoff;     // number of markers in all chunks expect for the last chunk
  // P2Mat should be of dimension: t_n * m1 
  // Suppose that 
  // n is the sample size in analysis 
  // m (<q) is the number of markers that pass the marker-level QC (e.g., g_missingRate_cutoff and g_region_maxMAF_cutoff)
  // VarMat (m x m) is the variance matrix of these m markers
  // VarMat = P1Mat %*% P2Mat, where P1Mat is of (m x n) and P2Mat is of (n x m)
  // Added on 09-17-2021: we collapse all ultra-rare variants (URV) to get one "fake" marker. 
  // That part has been moved to function mainRegionURVInCPP()
  
  std::vector<unsigned int> mPassCVVec;
    
  // conduct marker-level analysis
  std::string pval, pval_noSPA,  pval_c, pval_noSPA_c;
  double Beta, seBeta, Tstat, varT, gy;
  double Beta_c, seBeta_c, Tstat_c, varT_c;
  bool isSPAConverge, is_gtilde, is_Firth, is_FirthConverge;
  arma::vec P1Vec(t_n), P2Vec(t_n);
  arma::vec GVec(t_n);
  arma::vec GZeroVec(t_n);
  std::vector<uint> indexZeroVec;
  std::vector<uint> indexNonZeroVec;

  //variance ratio
  bool hasVarRatio = true;;
  bool isSingleVarianceRatio = true;

  if((ptr_gSAIGEobj->m_varRatio_null).n_elem > 1){
    isSingleVarianceRatio = false;
  }else{
    ptr_gSAIGEobj->assignSingleVarianceRatio(ptr_gSAIGEobj->m_flagSparseGRM_cur);
  }  
  // initiate chunk information
  unsigned int nchunks = 0; //number of chunks
  unsigned int ichunk = 0; //ith chunk
  unsigned int i1InChunk = 0; //i1th marker in ith chunk
  unsigned int i1 = 0;    // index of Markers (non-URV)
  unsigned int i2 = 0;    // index of Markers (Ultra-Rare Variants, URV)
  unsigned int jm; 

  double cctpval;
  double cctpval_cond;
  // cycle for q0 markers
  for(unsigned int i = 0; i < q0; i++)
  {
    // marker-level information
    double altFreq, altCounts, missingRate, imputeInfo;
    std::vector<uint32_t> indexForMissing;
    std::string chr, ref, alt, marker;
    uint32_t pd;
    bool flip = false;
    bool isOutputIndexForMissing = true;
    bool isOnlyOutputNonZero = false;

    GVec.resize(t_n);
    GVec.zeros();

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


    bool isReadMarker = Unified_getOneMarker(t_genoType, gIndex_prev, gIndex, ref, alt, marker, pd, chr, altFreq, altCounts, missingRate, imputeInfo,
                                          isOutputIndexForMissing, // bool t_isOutputIndexForMissing,
                                          indexForMissing,
                                          isOnlyOutputNonZero, // bool t_isOnlyOutputNonZero,
                                          indexNonZeroVec,
					  GVec,
					  t_isImputation);

   //arma::vec timeoutput2a = getTime();
   //printTime(timeoutput1a, timeoutput2a, "Unified_getOneMarker");

   if(!isReadMarker){
      std::cout << "ERROR: Reading " <<  i << "th marker failed." << std::endl;
      break;
    }	    
    std::string pds = std::to_string(pd);
    std::string info = chr+":"+std::to_string(pd)+":"+ref+":"+alt;

    double MAF = std::min(altFreq, 1 - altFreq);
    double w0;
    double MAC = MAF * 2 * t_n * (1 - missingRate);   // checked on 08-10-2021
    flip = imputeGenoAndFlip(GVec, altFreq, altCounts, indexForMissing, g_impute_method, g_dosage_zerod_cutoff, g_dosage_zerod_MAC_cutoff, MAC, indexZeroVec, indexNonZeroVec);
    arma::uvec indexZeroVec_arma, indexNonZeroVec_arma;
    MAF = std::min(altFreq, 1 - altFreq);
    MAC = std::min(altCounts, t_n *2 - altCounts);
    chrVec.at(i) = chr;
    posVec.at(i) = pds;
    refVec.at(i) = ref;
    altVec.at(i) = alt;
    markerVec.at(i) = marker;             // marker IDs
    infoVec.at(i) = info;                 // marker information: CHR:POS:REF:ALT
    altFreqVec.at(i) = altFreq;           // allele frequencies of ALT allele, this is not always < 0.5.
    missingRateVec.at(i) = missingRate;
    altCountsVec.at(i) = altCounts;
    MACVec.at(i) = MAC;
    MAFVec.at(i) = MAF;
    imputationInfoVec.at(i) = imputeInfo;

//arma::vec timeoutput3a = getTime();
    //printTime(timeoutput2a, timeoutput3a, "Unified_getOneMarker 2");
   if((missingRate > g_missingRate_cutoff) || (MAF > g_maxMAFLimit) || (MAF < g_marker_minMAF_cutoff) || (MAC < g_marker_minMAC_cutoff) || (imputeInfo < g_marker_minINFO_cutoff)){
      continue;
   }else{ 

    if(isWeightCustomized){
        w0 = t_weight(i);
    }else{
        w0 = boost::math::pdf(beta_dist, MAF);
    }


    indexNonZeroVec_arma = arma::conv_to<arma::uvec>::from(indexNonZeroVec);
    uint nNonZero = indexNonZeroVec_arma.n_elem;

    if(MAC > g_region_minMAC_cutoff){  // not Ultra-Rare Variants

      indicatorVec.at(i) = 1;      
      if(i1InChunk == 0){
        std::cout << "Start analyzing chunk " << ichunk << "....." << std::endl;
      }

      if(!isSingleVarianceRatio){
        hasVarRatio = ptr_gSAIGEobj->assignVarianceRatio(MAC, ptr_gSAIGEobj->m_flagSparseGRM_cur);
      }else{
        ptr_gSAIGEobj->assignSingleVarianceRatio(ptr_gSAIGEobj->m_flagSparseGRM_cur);
      }


      if(t_regionTestType != "BURDEN" || t_isSingleinGroupTest){ //perform single-variant assoc tests 
 
        indexZeroVec_arma = arma::conv_to<arma::uvec>::from(indexZeroVec);

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


        if(t_regionTestType != "BURDEN"){
          P1Mat.row(i1InChunk) = sqrt(ptr_gSAIGEobj->m_varRatioVal)*gtildeVec.t();
          P2Mat.col(i1InChunk) = sqrt(ptr_gSAIGEobj->m_varRatioVal)*P2Vec;
	}
     }//if(t_regionTestType != "BURDEN" || t_isSingleinGroupTest){ 

     i1 += 1;
     i1InChunk += 1;
     //arma::vec timeoutput3aa = getTime();
     //printTime(timeoutput3a, timeoutput3aa, "Unified_getOneMarker 3a");
     arma::vec dosage_case, dosage_ctrl;
     if(t_traitType == "binary"){
        dosage_case = GVec.elem(ptr_gSAIGEobj->m_case_indices);
        dosage_ctrl = GVec.elem(ptr_gSAIGEobj->m_ctrl_indices);
        MACcasegroup = arma::accu(dosage_case);
        MACcontrolgroup = arma::accu(dosage_ctrl);	
    }
      //arma::vec timeoutput3ab = getTime();
      //printTime(timeoutput3aa, timeoutput3ab, "Unified_getOneMarker 3b");
      //arma::vec MAFIndicatorVec(maxMAFVec.n_elem);
      MAFIndicatorVec.zeros();
      MAFIndicatorVec.elem( find(maxMAFVec >= MAF) ).ones();	
      annoMAFIndicatorVec.zeros();
      for(unsigned int j = 0; j < q_anno; j++){
        if(annoIndicatorMat(i,j) == 1){
		maxMAFperAnno(j) = std::max(maxMAFperAnno(j), MAF);
		for(unsigned int m = 0; m < q_maf; m++){
			if(MAFIndicatorVec(m) == 1){
  				//arma::vec timeoutput3ab0 = getTime();
				jm = j*q_maf + m;	
				annoMAFIndicatorVec(jm) = 1;
				MAC_GroupVec(jm) = MAC_GroupVec(jm) + MAC;

				if(t_traitType == "binary"){
					MACCase_GroupVec(jm) = MACCase_GroupVec(jm) + MACcasegroup;
					MACControl_GroupVec(jm) = MACControl_GroupVec(jm) + MACcontrolgroup;
					//genoSumMat.col(jm) = genoSumMat.col(jm) + w0*GVec;
				}
				for(unsigned int k = 0; k < nNonZero; k++){	
					genoSumMat(indexNonZeroVec_arma(k), jm) = genoSumMat(indexNonZeroVec_arma(k), jm) + w0*GVec(indexNonZeroVec_arma(k));

					genoSumcount_noweight(jm) = genoSumcount_noweight(jm) + GVec(indexNonZeroVec_arma(k));
				}
	
  //arma::vec timeoutput3ab2 = getTime();
   //      printTime(timeoutput3ab1, timeoutput3ab2, "Unified_getOneMarker 3b2");
				NumRare_GroupVec(jm) = NumRare_GroupVec(jm) + 1;
			}

		}	
        }
      }
  //arma::vec timeoutput3ac = getTime();
   //    printTime(timeoutput3ab, timeoutput3ac, "Unified_getOneMarker 3c");
     annoMAFIndicatorMat.row(i) = annoMAFIndicatorVec.t();
     if(t_regionTestType != "BURDEN" || t_isSingleinGroupTest){
      if(t_traitType == "binary"){
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

    }else{   // Ultra-Rare Variants (URV)
      indicatorVec.at(i) = 2;
      arma::vec MAFIndicatorVec(maxMAFVec.n_elem);
      MAFIndicatorVec.zeros();
      MAFIndicatorVec.elem( find(maxMAFVec >= MAF) ).ones();
      annoMAFIndicatorVec.zeros();
      for(unsigned int j = 0; j < q_anno; j++){
        if(annoIndicatorMat(i,j) == 1){
		maxMAFperAnno(j) = std::max(maxMAFperAnno(j), MAF);
		for(unsigned int m = 0; m < q_maf; m++){
                        if(MAFIndicatorVec(m) == 1){
                        	jm = j*q_maf + m;
				annoMAFIndicatorVec(jm) = 2;
				if(!isWeightCustomized){
				  for(unsigned int k = 0; k < nNonZero; k++){
					genoURMat(indexNonZeroVec_arma(k), jm) = std::max(genoURMat(indexNonZeroVec_arma(k), jm), GVec(indexNonZeroVec_arma(k)));
				  }
					//genoURMat.col(jm) = arma::max(genoURMat.col(jm), GVec);
				}else{

                                  for(unsigned int k = 0; k < nNonZero; k++){
					genoURMat(indexNonZeroVec_arma(k), jm) = std::max(genoURMat(indexNonZeroVec_arma(k), jm) , t_weight(i) * (GVec(indexNonZeroVec_arma(k))));
					//weightURMat_cnt(indexNonZeroVec_arma(k), jm) = weightURMat_cnt(indexNonZeroVec_arma(k), jm) + 1;
				  }
				}	
				NumUltraRare_GroupVec(jm) = NumUltraRare_GroupVec(jm) + 1;
			}
		}
	}
      }
      annoMAFIndicatorMat.row(i) = annoMAFIndicatorVec.t();

      i2 += 1;
    }
 }//else if((missingRate > g_missingRate_cutoff) || (MAF > g_maxMAFLimit) || (MAF < g_marker_minMAF_cutoff) || (MAC < g_marker_minMAC_cutoff) || (imputeInfo < g_marker_minINFO_cutoff)){
//
    
    if(i1InChunk == m1 ){
      std::cout << "In chunks 0-" << ichunk << ", " << i2 << " markers are ultra-rare and " << i1 << " markers are not ultra-rare." << std::endl;
      if(t_regionTestType != "BURDEN"){
        P1Mat.save(t_outputFile + "_P1Mat_Chunk_" + std::to_string(ichunk) + ".bin");
        P2Mat.save(t_outputFile + "_P2Mat_Chunk_" + std::to_string(ichunk) + ".bin");
      }

      mPassCVVec.push_back(m1);
      ichunk += 1;
      i1InChunk = 0;
      nchunks = nchunks + 1;
    }
    Rcpp::checkUserInterrupt();


} //  for(unsigned int i = 0; i < q0; i++)



//the second last chunk
  if(i1InChunk != 0){
    std::cout << "In chunks 0-" << ichunk << ", " << i2 << " markers are ultra-rare and " << i1 << " markers are not ultra-rare." << std::endl;
    if(t_regionTestType != "BURDEN"){	  
      P1Mat = P1Mat.rows(0, i1InChunk - 1);
      P2Mat = P2Mat.cols(0, i1InChunk - 1);
      P1Mat.save(t_outputFile + "_P1Mat_Chunk_" + std::to_string(ichunk) + ".bin");
      P2Mat.save(t_outputFile + "_P2Mat_Chunk_" + std::to_string(ichunk) + ".bin");
    }
    ichunk = ichunk + 1;
    //}
    mPassCVVec.push_back(i1InChunk);
    nchunks = nchunks + 1; 
    i1InChunk = 0;
  }


//for all UR variants
if(i2 > 0){
  int m1new = std::max(m1, q_anno_maf);
  //if(isWeightCustomized){
    //weightURMat_cnt.replace(0, 1); 	  
    //genoURMat = genoURMat / weightURMat_cnt;
    //weightURMat.ones();    
  //}

  if(t_regionTestType != "BURDEN"){
    P1Mat.resize(m1new, P1Mat.n_cols);
    P2Mat.resize(P2Mat.n_rows,m1new);
  }

  arma::mat XV, XXVX_inv;
  ptr_gSAIGEobj->extract_XV_XXVX_inv(XV, XXVX_inv);
  //the last chunk for UR
  unsigned int i;
  for(unsigned int j = 0; j < q_anno; j++){
     for(unsigned int m = 0; m < q_maf; m++){
	jm = j*q_maf+m;
	arma::vec genoURVec = genoURMat.col(jm);
	int n = genoURVec.size();
	arma::uvec indexForNonZero = arma::find(genoURVec != 0);
	i = q0 + jm;
	markerVec.at(i) = "UR";             // marker IDs
	if(indexForNonZero.n_elem > 0){
	//URindVec.push_back(jm+1);
	  double altFreq = arma::mean(genoURVec)/2;
	  double altCounts = arma::accu(genoURVec);
	  double missingRate = 0;
	  double imputeInfo = 1;
    	  std::string chr, ref, alt, marker;
    	  bool flip = false;
	  std::string info = "UR";	
    	  double MAF = std::min(altFreq, 1 - altFreq);
	  double w0;
	  double MAC = MAF * 2 * t_n * (1 - missingRate);
	  std::vector<uint32_t> indexForMissing;
    	  flip = imputeGenoAndFlip(genoURVec, altFreq, altCounts, indexForMissing, g_impute_method, g_dosage_zerod_cutoff, g_dosage_zerod_MAC_cutoff, MAC, indexZeroVec, indexNonZeroVec);
	  if(isWeightCustomized){
	    //genoSumMat.col(jm) = genoSumMat.col(jm) + genoURVec % (weightURMat.col(jm));
	    //genoSumMat.col(jm) = genoSumMat.col(jm) + genoURVec;
	      for(unsigned int k = 0; k < indexForNonZero.n_elem; k++){
                genoSumMat(indexForNonZero(k), jm) = genoSumMat(indexForNonZero(k), jm) + genoURVec(indexForNonZero(k));
	 	genoSumcount_noweight(jm) = genoSumcount_noweight(jm) + genoURVec(indexForNonZero(k));
	      }
    	  }else{
               w0 = boost::math::pdf(beta_dist, MAF);
	       for(unsigned int k = 0; k < indexForNonZero.n_elem; k++){
                genoSumMat(indexForNonZero(k), jm) = genoSumMat(indexForNonZero(k), jm) + genoURVec(indexForNonZero(k)) * w0;
	 	genoSumcount_noweight(jm) = genoSumcount_noweight(jm) + genoURVec(indexForNonZero(k));
	      }
    	  }

	 //genoSumMat.col(1).print("genoSumMat.col(1)");


  	 if(t_regionTestType != "BURDEN"){
	    arma::vec genoSumMatvec1 = genoSumMat.col(jm);
	    arma::vec genoSumMatvec2 = XV * genoSumMatvec1;
	    arma::vec genoSumMatvec3 = genoSumMatvec1 - XXVX_inv * genoSumMatvec2;
	    genoSumMat.col(jm) = genoSumMatvec3;
          }

	  MAC = MAF * 2 * t_n * (1 - missingRate);   // checked on 08-10-2021


      if(!isSingleVarianceRatio){
        hasVarRatio = ptr_gSAIGEobj->assignVarianceRatio(MAC, ptr_gSAIGEobj->m_flagSparseGRM_cur);
      }else{
        ptr_gSAIGEobj->assignSingleVarianceRatio(ptr_gSAIGEobj->m_flagSparseGRM_cur);
      }


	  if(t_regionTestType != "BURDEN" || t_isSingleinGroupTest){


	    annoMAFIndicatorVec.zeros();
	    annoMAFIndicatorVec(jm) = 1;
	    annoMAFIndicatorMat.row(i) = annoMAFIndicatorVec.t();
    	    infoVec.at(i) = info;                 // marker information: CHR:POS:REF:ALT
    	    altFreqVec.at(i) = altFreq;	// allele frequencies of ALT allele, this is not always < 0.5.
	    altCountsVec.at(i) = altCounts;
	    missingRateVec.at(i) = missingRate;
    	    MACVec.at(i) = MAC;
    	    MAFVec.at(i) = MAF;

            arma::uvec indexZeroVec_arma, indexNonZeroVec_arma;
            indexZeroVec_arma = arma::conv_to<arma::uvec>::from(indexZeroVec);
            indexNonZeroVec_arma = arma::conv_to<arma::uvec>::from(indexNonZeroVec);

	    if(MAC <= g_MACCutoffforER && t_traitType == "binary"){	

              ptr_gSAIGEobj->getMarkerPval(genoURVec, indexNonZeroVec_arma, indexZeroVec_arma, Beta, seBeta, pval, pval_noSPA, altFreq, Tstat, gy, varT, isSPAConverge, gtildeVec, is_gtilde, true, P2Vec, isCondition, Beta_c, seBeta_c, pval_c, pval_noSPA_c, Tstat_c, varT_c, G1tilde_P_G2tilde_Vec, is_Firth, is_FirthConverge, true);
	    }else{
              ptr_gSAIGEobj->getMarkerPval(genoURVec, indexNonZeroVec_arma, indexZeroVec_arma, Beta, seBeta, pval, pval_noSPA, altFreq, Tstat, gy, varT, isSPAConverge, gtildeVec, is_gtilde, true, P2Vec, isCondition, Beta_c, seBeta_c, pval_c, pval_noSPA_c, Tstat_c, varT_c, G1tilde_P_G2tilde_Vec, is_Firth, is_FirthConverge, false);

	    }

            BetaVec.at(i) = Beta* (1 - 2*flip);
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

    	    chrVec.at(i) = "UR";
    	    posVec.at(i) = "UR";
    	    refVec.at(i) = "UR";
    	    altVec.at(i) = "UR";


	    std::string str = std::to_string(maxMAFVec.at(m));
	    str.erase ( str.find_last_not_of('0') + 1, std::string::npos );
    	    markerVec.at(i) = regionName + ":" + annoStringVec.at(j) + ":" + str ;
            arma::vec dosage_case, dosage_ctrl;
            MAC_GroupVec(jm) = MAC_GroupVec(jm) + MAC;
            if(t_traitType == "binary"){
                        dosage_case = genoURVec.elem(ptr_gSAIGEobj->m_case_indices);
                        dosage_ctrl = genoURVec.elem(ptr_gSAIGEobj->m_ctrl_indices);
                        MACcasegroup = arma::accu(dosage_case);
                        MACcontrolgroup = arma::accu(dosage_ctrl);
                        MACCase_GroupVec(jm) = MACCase_GroupVec(jm) + MACcasegroup;
                        MACControl_GroupVec(jm) = MACControl_GroupVec(jm) + MACcontrolgroup;
            }

          if(t_traitType == "binary"){
            AF_case = arma::mean(dosage_case) /2;
            AF_ctrl = arma::mean(dosage_ctrl) /2;
            if(flip){
              AF_case = 1-AF_case;
              AF_ctrl = 1-AF_ctrl;
            }
          AF_caseVec.at(i) = AF_case;
          AF_ctrlVec.at(i) = AF_ctrl;
          N_caseVec.at(i) = dosage_case.n_elem;
          N_ctrlVec.at(i) = dosage_ctrl.n_elem;

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
          N_Vec.at(i) = n;
        }
      if(t_regionTestType != "BURDEN"){	
        P1Mat.row(i1InChunk) = sqrt(ptr_gSAIGEobj->m_varRatioVal)*gtildeVec.t();
	//gtildeVec.print("gtildeVec");
        P2Mat.col(i1InChunk) = sqrt(ptr_gSAIGEobj->m_varRatioVal)*P2Vec;
	//P2Vec.print("P2Vec");

      } //if(t_regionTestType != "BURDEN"){
    }else{//if(t_regionTestType != "BURDEN" || t_isSingleinGroupTest){	
            MAC_GroupVec(jm) = MAC_GroupVec(jm) + MAC;
	                arma::vec dosage_case, dosage_ctrl;
            if(t_traitType == "binary"){
                        dosage_case = genoURVec.elem(ptr_gSAIGEobj->m_case_indices);
                        dosage_ctrl = genoURVec.elem(ptr_gSAIGEobj->m_ctrl_indices);
                        MACcasegroup = arma::accu(dosage_case);
                        MACcontrolgroup = arma::accu(dosage_ctrl);
                        MACCase_GroupVec(jm) = MACCase_GroupVec(jm) + MACcasegroup;
                        MACControl_GroupVec(jm) = MACControl_GroupVec(jm) + MACcontrolgroup;
            }
    }
  
    i1InChunk = i1InChunk + 1;
    i1 = i1 + 1;
    }
   }  
  }


  if(i1InChunk != 0){
    nchunks = nchunks + 1;
    if(t_regionTestType != "BURDEN"){
      P1Mat = P1Mat.rows(0, i1InChunk - 1);
      P2Mat = P2Mat.cols(0, i1InChunk - 1);
//if(nchunks != 1){
      P1Mat.save(t_outputFile + "_P1Mat_Chunk_" + std::to_string(ichunk) + ".bin");
      P2Mat.save(t_outputFile + "_P2Mat_Chunk_" + std::to_string(ichunk) + ".bin");
    }
      ichunk = ichunk + 1;
//    }
    mPassCVVec.push_back(i1InChunk);
  }
   //std::cout << "P1Mat.n_rows ok2 " << P1Mat.n_rows << std::endl; 

  }// if(i2 > 0)    
  int mPassCVVecsize = mPassCVVec.size();
  nchunks = mPassCVVecsize;


arma::mat VarMat;
//(i1, i1);

if(t_regionTestType != "BURDEN"){
  VarMat.resize(i1, i1);	
  if(nchunks == 1){
    VarMat = P1Mat * P2Mat;
  }

  // the region includes more markers than limitation, so P1Mat and P2Mat have been put in hard drive
  if(nchunks > 1)
  {
    int first_row = 0, first_col = 0, last_row = 0, last_col = 0;
    
    for(unsigned int index1 = 0; index1 < nchunks; index1++)
    {
      last_row = first_row + mPassCVVec.at(index1) - 1;
      
      std::string P1MatFile = t_outputFile + "_P1Mat_Chunk_" + std::to_string(index1) + ".bin";
      
        P1Mat.load(P1MatFile);
   
       //P1Mat.print("P1Mat");

      if(P1Mat.n_cols == 0) continue;
      
      // off-diagonal sub-matrix
      for(unsigned int index2 = 0; index2 < index1; index2++)
      {
        std::cout << "Analyzing chunks (" << index1 << "/" << nchunks - 1 << ", " << index2 << "/" << nchunks - 1 << ")........" << std::endl;
       
        P2Mat.load(t_outputFile + "_P2Mat_Chunk_" + std::to_string(index2) + ".bin");
       //P2Mat.print("P2Mat");
        
        if(P2Mat.n_cols == 0) continue;
        arma::mat offVarMat = P1Mat * P2Mat;
        last_col = first_col + mPassCVVec.at(index2) - 1;
        
        VarMat.submat(first_row, first_col, last_row, last_col) = offVarMat;
        VarMat.submat(first_col, first_row, last_col, last_row) = offVarMat.t();
        first_col = last_col + 1;
      }
      
      // diagonal sub-matrix
      last_col = first_col + mPassCVVec.at(index1) - 1;
      std::cout << "Analyzing chunks (" << index1 << "/" << nchunks - 1 << ", " << index1 << "/" << nchunks - 1 << ")........" << std::endl;
      //std::cout << "P2Mat.n_cols " << P2Mat.n_cols << std::endl;
      P2Mat.load(t_outputFile + "_P2Mat_Chunk_" + std::to_string(index1) + ".bin");
      
      arma::mat diagVarMat = P1Mat * P2Mat;
      
      VarMat.submat(first_row, first_col, last_row, last_col) = diagVarMat;
      
      first_row = last_row + 1;
      first_col = 0;
      Rcpp::checkUserInterrupt();
    }
     
    for(unsigned int index1 = 0; index1 < nchunks; index1++)
    {
      std::string P1MatFile = t_outputFile + "_P1Mat_Chunk_" + std::to_string(index1) + ".bin";
      std::string P2MatFile = t_outputFile + "_P2Mat_Chunk_" + std::to_string(index1) + ".bin";
      const char* File1 = P1MatFile.c_str();
      const char* File2 = P2MatFile.c_str();
      std::remove(File1);
      std::remove(File2);
    }
    
  }
}  //if(t_regionTestType != "BURDEN"){


//read and test single markers done
//arma::vec timeoutput2 = getTime();
//printTime(timeoutput1, timeoutput2, "read and test single markers done");
//check the max MAF of all markers
//if there are sets containing the same markers
arma::uvec q_maf_for_anno(q_anno);
for(unsigned int j = 0; j < q_anno; j++){
	arma::uvec jtemp = find(maxMAFVec >= maxMAFperAnno(j));
	q_maf_for_anno(j) = jtemp.min();
}

//If only conduct Burden test
//arma::vec BURDEN_pval_Vec(q_anno_maf);
std::vector<std::string> BURDEN_pval_Vec(q_anno_maf, "NA");
//BURDEN_pval_Vec.fill(-1.0);
//arma::vec BURDEN_pval_cVec(q_anno_maf);
std::vector<std::string> BURDEN_pval_cVec(q_anno_maf, "NA");
//BURDEN_pval_cVec.fill(-1.0);
std::vector<std::string> BURDEN_AnnoName_Vec(q_anno_maf);
std::vector<std::string> BURDEN_maxMAFName_Vec(q_anno_maf);
std::vector<double> BURDEN_Beta_Vec(q_anno_maf);
std::vector<double> BURDEN_seBeta_Vec(q_anno_maf);
//std::vector<double> BURDEN_pval_cVec(q_anno_maf);
std::vector<double> BURDEN_Beta_cVec(q_anno_maf);
std::vector<double> BURDEN_seBeta_cVec(q_anno_maf);



bool iswriteOutput = false;
bool isregion = true;
if(!ptr_gSAIGEobj->m_flagSparseGRM){
	isregion = false;
}
//Rcpp::DataFrame OUT_BURDEN = Rcpp::DataFrame::create();
unsigned int i= 0;
unsigned int q_maf_m;
bool isPolyMarker = true;
std::string AnnoName;
double maxMAFName;
if(t_regionTestType == "BURDEN"){
     for(unsigned int j = 0; j < q_anno; j++){
       q_maf_m = q_maf_for_anno(j);
       AnnoName = annoStringVec[j];
       isPolyMarker = true;	
       for(unsigned int m = 0; m < q_maf; m++){
	maxMAFName = maxMAFVec(m); 
	jm = j*q_maf+m;
	i = jm;
      if(m <= q_maf_m){
        arma::vec genoSumVec = genoSumMat.col(jm);
        int n = genoSumVec.size();
        arma::uvec indexNonZeroVec_arma = arma::find(genoSumVec != 0);
	arma::uvec indexZeroVec_arma = arma::find(genoSumVec == 0);
        //double altFreq = arma::mean(genoSumVec)/2;
        //double altCounts = arma::accu(genoSumVec);
	double altCounts =  genoSumcount_noweight(i);
	double altFreq = altCounts/(2*t_n);
	double MAC, MAF;
	if(altFreq > 1){
		MAF = 1;
		MAC = t_n;
	}else{
		MAF = std::min(altFreq, 1 - altFreq);
		MAC = MAF*2*t_n;
	}
	//std::cout << "altCounts " << altCounts << std::endl;
	//std::cout << "altFreq " << altCounts << std::endl;
        double missingRate = 0;
        double imputeInfo = 1;
        std::string chr, ref, alt, marker;
        bool flip = false;
        std::string info = "UR";
        //double MAF = std::min(altFreq, 1 - altFreq);
        double w0;
        //double MAC = MAF * 2 * t_n * (1 - missingRate);
        if(indexNonZeroVec_arma.n_elem > 0 && MAC >= g_min_gourpmac_for_burdenonly){
	  //if(MAC >= g_min_gourpmac_for_burdenonly){
          isPolyMarker = true;   
	  std::vector<uint32_t> indexForMissing;

          if(!isSingleVarianceRatio){
            hasVarRatio = ptr_gSAIGEobj->assignVarianceRatio(MAC, ptr_gSAIGEobj->m_flagSparseGRM_cur);
          }
	  //arma::vec timeoutput_getp = getTime();
	  if(MAC <= g_MACCutoffforER && t_traitType == "binary"){
          ptr_gSAIGEobj->getMarkerPval(genoSumVec, indexNonZeroVec_arma, indexZeroVec_arma, Beta, seBeta, pval, pval_noSPA, altFreq, Tstat, gy, varT, isSPAConverge, gtildeVec, is_gtilde, isregion, P2Vec, isCondition, Beta_c, seBeta_c, pval_c, pval_noSPA_c, Tstat_c, varT_c, G1tilde_P_G2tilde_Vec, is_Firth, is_FirthConverge, true);
	  }else{
          ptr_gSAIGEobj->getMarkerPval(genoSumVec, indexNonZeroVec_arma, indexZeroVec_arma, Beta, seBeta, pval, pval_noSPA, altFreq, Tstat, gy, varT, isSPAConverge, gtildeVec, is_gtilde, isregion, P2Vec, isCondition, Beta_c, seBeta_c, pval_c, pval_noSPA_c, Tstat_c, varT_c, G1tilde_P_G2tilde_Vec, is_Firth, is_FirthConverge, false);
	  }	  
	  //arma::vec timeoutput_getp2 = getTime();
	  //printTime(timeoutput_getp, timeoutput_getp2, "get p  done");
	  
	  if(isCondition){
	    BURDEN_pval_cVec.at(i) = pval_c;
	    BURDEN_Beta_cVec.at(i) = Beta_c;
	    BURDEN_seBeta_cVec.at(i) = seBeta_c;
          }

	   BURDEN_AnnoName_Vec.at(i) = AnnoName;
	   std::string str = std::to_string(maxMAFName);
	   str.erase ( str.find_last_not_of('0') + 1, std::string::npos );
	   BURDEN_maxMAFName_Vec.at(i) = str;
	   BURDEN_pval_Vec.at(i) = pval;
	   BURDEN_Beta_Vec.at(i) = Beta;
	   BURDEN_seBeta_Vec.at(i) = seBeta;

         }else{
	   isPolyMarker = false;	 
	 }
     }else{
	if(isPolyMarker){
	  if(isCondition){
            BURDEN_pval_cVec.at(i) = pval_c;
            BURDEN_Beta_cVec.at(i) = Beta_c;
            BURDEN_seBeta_cVec.at(i) = seBeta_c;
          }

           BURDEN_AnnoName_Vec.at(i) = AnnoName;
           BURDEN_maxMAFName_Vec.at(i) = std::to_string(maxMAFName);
           BURDEN_pval_Vec.at(i) = pval;
           BURDEN_Beta_Vec.at(i) = Beta;
           BURDEN_seBeta_Vec.at(i) = seBeta; 
	 }
     } 	     

   }
 }
	  std::vector<double> nonMissingPvalVec_std, nonMissingPvalVec_c_std;
	  double burden_p, burden_p_cond;
	  for(unsigned int i = 0; i < BURDEN_pval_Vec.size(); i++){
		if(BURDEN_pval_Vec.at(i) != "NA"){
	  		burden_p = std::stod(BURDEN_pval_Vec.at(i));
			nonMissingPvalVec_std.push_back(burden_p);
			if(isCondition){
				burden_p_cond = std::stod(BURDEN_pval_cVec.at(i));
				nonMissingPvalVec_c_std.push_back(burden_p_cond);
			}
		}	
	  } 

	  arma::vec nonMissingPvalVec = arma::conv_to< arma::vec >::from(nonMissingPvalVec_std);
	  //arma::uvec nonMissingPvalVecInd = arma::find(BURDEN_pval_Vec >= 0);
	  //arma::vec nonMissingPvalVec = BURDEN_pval_Vec.elem(nonMissingPvalVecInd);
	   cctpval = CCT_cpp(nonMissingPvalVec);
           if(isCondition){
	        //arma::vec nonMissingPvalVec_cond = BURDEN_pval_cVec.elem(nonMissingPvalVecInd);
		arma::vec nonMissingPvalVec_cond  = arma::conv_to< arma::vec >::from(nonMissingPvalVec_c_std);
	   	cctpval_cond = CCT_cpp(nonMissingPvalVec_cond);	   
           }

	
	   if(!t_isFastTest){
		iswriteOutput = true;
	   }else{
		if(ptr_gSAIGEobj->m_flagSparseGRM_cur){
			iswriteOutput = true;
		}else{
		  if(cctpval >= 0.1){
			iswriteOutput = true;
		  } 	
		}
	   }
	
	   if(iswriteOutput){
		writeOutfile_BURDEN(regionName,
			BURDEN_AnnoName_Vec,
			BURDEN_maxMAFName_Vec,
			BURDEN_pval_Vec,
			BURDEN_Beta_Vec,
			BURDEN_seBeta_Vec,
			BURDEN_pval_cVec,
			BURDEN_Beta_cVec,
			BURDEN_seBeta_cVec,
			MAC_GroupVec,
			MACCase_GroupVec,
			MACControl_GroupVec,
			NumRare_GroupVec,
			NumUltraRare_GroupVec,
			cctpval,
			cctpval_cond,
			q_anno,
			q_maf,
			isCondition,
			t_traitType);
	   }else{
	     OutList.push_back(iswriteOutput, "iswriteOutput");
	   }

 //arma::vec timeoutput3 = getTime();
 //printTime(timeoutput2, timeoutput3, "burden test done");
 }else{
  q_maf_for_anno = q_maf_for_anno + 1;
  OutList.push_back(MAC_GroupVec, "MAC_GroupVec");
  OutList.push_back(q_maf_for_anno, "q_maf_for_annoVec");
  if(t_traitType == "binary"){
    OutList.push_back(MACCase_GroupVec, "MACCase_GroupVec");
    OutList.push_back(MACControl_GroupVec, "MACCtrl_GroupVec");
    OutList.push_back(genoSumMat, "genoSumMat");
    OutList.push_back(gyVec, "gyVec");
  }

    OutList.push_back(VarMat, "VarMat");	
    OutList.push_back(MAFVec, "MAFVec");	
    OutList.push_back(TstatVec_flip, "TstatVec_flip");	
  //arma::mat scaled_m_VarInvMat_cond;
    if(isCondition){
  //std::cout << "okk5" << std::endl;
      arma::mat AdjCondMat = G1tilde_P_G2tilde_Weighted_Mat * (ptr_gSAIGEobj->m_VarInvMat_cond / (w0G2Mat_cond));
      arma::mat VarMatAdjCond = AdjCondMat * (G1tilde_P_G2tilde_Weighted_Mat.t());
      arma::vec TstatAdjCond = AdjCondMat * (ptr_gSAIGEobj->m_Tstat_cond % w0G2Vec_cond ); 
      OutList.push_back(G1tilde_P_G2tilde_Weighted_Mat, "G1tilde_P_G2tilde_Weighted_Mat"); 
      OutList.push_back(ptr_gSAIGEobj->m_scalefactor_G2_cond, "scalefactor_G2_cond");
      OutList.push_back(ptr_gSAIGEobj->m_VarInvMat_cond_scaled_weighted, "VarInvMat_G2_cond_scaled"); 
      OutList.push_back(ptr_gSAIGEobj->m_Tstat_cond, "Tstat_G2_cond"); //m_Tstat_cond is weighted
      OutList.push_back(ptr_gSAIGEobj->m_G2_Weight_cond, "G2_Weight_cond");
      OutList.push_back(TstatAdjCond, "TstatAdjCond");
      OutList.push_back(VarMatAdjCond, "VarMatAdjCond"); 
    }

    if(!t_isFastTest){
            iswriteOutput = true;
     }else{
        if(ptr_gSAIGEobj->m_flagSparseGRM_cur){
            iswriteOutput = true;
        }
     }
}


 int numofUR = q_anno_maf;
 int numofUR0;
 int mFirth = 0;
if(t_isSingleinGroupTest){
 // OutList.push_back(pvalVec_val, "pvalVec");
 OutList.push_back(pvalVec, "pvalVec");
if(iswriteOutput){
  numofUR0 = writeOutfile_singleInGroup(t_isMoreOutput,
      t_isImputation,
      isCondition,
      is_Firth,
      mFirth,
      is_FirthConverge,
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
      OutFile_singleInGroup);
 }else{
  OutFile_singleInGroup_temp.open(g_outputFilePrefixSingleInGroup_temp.c_str(), std::ofstream::out);
  numofUR0 = writeOutfile_singleInGroup(t_isMoreOutput,
      t_isImputation,
      isCondition,
      is_Firth,
      mFirth,
      is_FirthConverge,
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
      OutFile_singleInGroup_temp);
  OutFile_singleInGroup_temp.close();

 }//iswriteOutput
 OutList.push_back(numofUR, "numofUR");
}

OutList.push_back(iswriteOutput, "iswriteOutput");

 if(t_isOutputMarkerList){
	OutList.push_back(indicatorVec, "markerIndcatorVec");
 }

  

  OutList.push_back(NumRare_GroupVec, "NumRare_GroupVec");
  OutList.push_back(NumUltraRare_GroupVec, "NumUltraRare_GroupVec");


 if(t_regionTestType != "BURDEN" || t_isOutputMarkerList){
  OutList.push_back(annoMAFIndicatorMat, "annoMAFIndicatorMat");
 }
  
  return OutList;
}








// [[Rcpp::export]]
void assign_conditionMarkers_factors(
                           std::string t_genoType,     // "plink", "bgen", "vcf"
			   std::vector<std::string> & t_genoIndex_prev,
                           std::vector<std::string> & t_genoIndex,
                           unsigned int t_n, 
			   arma::vec & t_weight_cond
			   )           // sample size
{
  ptr_gSAIGEobj->set_flagSparseGRM_cur(ptr_gSAIGEobj->m_flagSparseGRM);
  bool isImpute = false;	
  unsigned int q = t_genoIndex.size();
  arma::mat P1Mat(q, t_n);
  arma::mat P2Mat(t_n, q);
  arma::mat VarInvMat(q, q);
  arma::vec TstatVec(q);
  //arma::vec pVec(q);
    std::vector<std::string> pVec(q, "NA");
  arma::vec MAFVec(q);
  arma::vec gyVec(q);
  arma::vec w0G2_cond_Vec(q);
  arma::vec gsumVec(t_n, arma::fill::zeros);
  //double beta1 = g_weights_beta[0];
  //double beta2 = g_weights_beta[1];  
  boost::math::beta_distribution<> beta_dist(g_weights_beta[0], g_weights_beta[1]);
  //boost::math::beta_distribution<> beta_dist(beta1, beta2);
  //g_weights_beta.print();
  //boost::math::beta_distribution<> beta_dist(1, 25);
  //std::vector<double> GVec0(t_n);
  arma::vec GVec(t_n);
  std::string pval, pval_noSPA;
  double Beta, seBeta, Tstat, varT, gy, w0G2_cond;
  bool isSPAConverge, is_gtilde, is_Firth, is_FirthConverge;
  arma::vec P2Vec(t_n);

  //std::vector<uint> indexZeroVec;
  //std::vector<uint> indexNonZeroVec;
  std::string pval_c, pval_noSPA_c;
  double Beta_c, seBeta_c, Tstat_c, varT_c;
  arma::rowvec G1tilde_P_G2tilde_Vec;
  bool isCondition = false;
  for(unsigned int i = 0; i < q; i++)
  {
    // marker-level information
    double altFreq, altCounts, missingRate, imputeInfo;
    std::vector<uint32_t> indexForMissing;
    std::vector<uint> indexZeroVec;
    std::vector<uint> indexNonZeroVec;
    std::string chr, ref, alt, marker;
    uint32_t pd;
    bool flip = false;

    bool isOutputIndexForMissing = true;
    bool isOnlyOutputNonZero = false; 




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
            gIndex_prev = std::strtoull( t_genoIndex_prev_str.c_str(), &end_prev,10 );
        }else if(t_genoType == "plink"){
            t_genoIndex_prev_str = t_genoIndex.at(i-1);
            gIndex_prev = std::strtoull( t_genoIndex_prev_str.c_str(), &end_prev,10 );
        }
	//else if(t_genoType == "vcf"){
         //   gIndex_prev = 0;	
	//}
	std::remove(end_prev);
    }

    bool isReadMarker = Unified_getOneMarker(t_genoType, gIndex_prev, gIndex, ref, alt, marker, pd, chr, altFreq, altCounts, missingRate, imputeInfo,
                                          isOutputIndexForMissing, // bool t_isOutputIndexForMissing,
                                          indexForMissing,
                                          isOnlyOutputNonZero, // bool t_isOnlyOutputNonZero,
                                          indexNonZeroVec, GVec, isImpute);
     //arma::vec GVec(GVec0);
     //GVec0.clear();
    if(!isReadMarker){
      break;
    }

    std::string info = chr+":"+std::to_string(pd)+"_"+ref+"/"+alt;

  double MAF = std::min(altFreq, 1 - altFreq);
  double MAC = MAF * 2 * t_n * (1 - missingRate);

  bool hasVarRatio;
  if((ptr_gSAIGEobj->m_varRatio_null).n_elem == 1){
	ptr_gSAIGEobj->assignSingleVarianceRatio(ptr_gSAIGEobj->m_flagSparseGRM);	
  }else{
	hasVarRatio = ptr_gSAIGEobj->assignVarianceRatio(MAC, ptr_gSAIGEobj->m_flagSparseGRM);
	//if(!hasVarRatio){
	//	std::cout << "Error! Conditioning marker " << info << " has MAC " << MAC << " and does not have variance ratio estimated." << std::endl;
	//	exit(EXIT_FAILURE);
	//}	
  }	  
 
  flip = imputeGenoAndFlip(GVec, altFreq, altCounts, indexForMissing, g_impute_method, g_dosage_zerod_cutoff, g_dosage_zerod_MAC_cutoff, MAC, indexZeroVec, indexNonZeroVec);


 arma::uvec indexZeroVec_arma, indexNonZeroVec_arma;
       indexZeroVec_arma = arma::conv_to<arma::uvec>::from(indexZeroVec);
       indexNonZeroVec_arma = arma::conv_to<arma::uvec>::from(indexNonZeroVec);


  MAF = std::min(altFreq, 1 - altFreq);
  MAC = std::min(altCounts, 2*t_n-altCounts);

  arma::vec gtildeVec;

  if(MAC > g_MACCutoffforER){
     Unified_getMarkerPval(
                    GVec,
                    false, // bool t_isOnlyOutputNonZero,
                    indexNonZeroVec_arma, indexZeroVec_arma, Beta, seBeta, pval, pval_noSPA,  Tstat, gy, varT, altFreq, isSPAConverge, gtildeVec, is_gtilde, true, P2Vec, isCondition, Beta_c, seBeta_c, pval_c, pval_noSPA_c, Tstat_c, varT_c, G1tilde_P_G2tilde_Vec, is_Firth, is_FirthConverge, false);
  }else{
     Unified_getMarkerPval(
                    GVec,
                    false, // bool t_isOnlyOutputNonZero,
                    indexNonZeroVec_arma, indexZeroVec_arma, Beta, seBeta, pval, pval_noSPA,  Tstat, gy, varT, altFreq, isSPAConverge, gtildeVec, is_gtilde, true, P2Vec, isCondition, Beta_c, seBeta_c, pval_c, pval_noSPA_c, Tstat_c, varT_c, G1tilde_P_G2tilde_Vec, is_Firth, is_FirthConverge, true);

  }
      P1Mat.row(i) = sqrt(ptr_gSAIGEobj->m_varRatioVal)*gtildeVec.t();
      P2Mat.col(i) = sqrt(ptr_gSAIGEobj->m_varRatioVal)*P2Vec;
      //P1Mat.row(i) = gtildeVec.t();
      //P2Mat.col(i) = P2Vec;
      MAFVec(i) = MAF;
      //w0G2_cond = boost::math::pdf(beta_dist, MAF);
      //t_weight_cond.print();
      if(MAF == 0.0){
	std::cerr << "ERROR: Conditioning marker is monomorphic\n";
      }	      


     if(!t_weight_cond.is_zero()){
	 w0G2_cond = t_weight_cond(i);
    }else{
	 w0G2_cond = boost::math::pdf(beta_dist, MAF);
    }
     w0G2_cond_Vec(i) = w0G2_cond;
     gyVec(i) = gy * w0G2_cond;
     gsumVec = gsumVec + GVec * w0G2_cond;
     TstatVec(i) = Tstat;
     //pVec(i) = std::stod(pval);
     pVec.at(i) = pval;
  }
  arma::mat VarMat = P1Mat * P2Mat;

  VarInvMat =  arma::pinv(VarMat);
  //VarInvMat = VarMat.i();
  double qsum = arma::accu(gyVec);
  arma::vec gsumtildeVec;


  ptr_gSAIGEobj->getadjG(gsumVec, gsumtildeVec);
  ptr_gSAIGEobj->assignConditionFactors(
		   			P2Mat,
					VarInvMat,
					VarMat,
					TstatVec,
				        w0G2_cond_Vec,	
					MAFVec,
					qsum,
					gsumtildeVec,
					pVec);

}

// [[Rcpp::export]]
void assign_conditionMarkers_factors_binary_region(
			   arma::vec & scalefactor_G2_cond){
	//std::cout << "assign_conditionMarkers_factors_binary_region" << std::endl;
	ptr_gSAIGEobj->assignConditionFactors_scalefactor(scalefactor_G2_cond);
}

// [[Rcpp::export]]
void set_iterator_inVcf(std::string & variantList, std::string & chrom, int & beg_pd, int & end_pd){
   if(!variantList.empty()){
	ptr_gVCFobj->set_iterator(variantList);	
   }else{
	ptr_gVCFobj->set_iterator(chrom, beg_pd, end_pd);
   }	   
}	

// [[Rcpp::export]]
bool check_Vcf_end(){
	bool isEnd = false;
	isEnd = ptr_gVCFobj->check_iterator_end();
	return(isEnd);
}


// [[Rcpp::export]]
void move_forward_iterator_Vcf(int i){
	ptr_gVCFobj->move_forward_iterator(i);
}



// [[Rcpp::export]]
arma::vec fast_logistf_fit(arma::mat & x,
		arma::vec & y,
		arma::vec & weight,
		arma::vec & offset,
		bool firth,
		arma::uvec & col_fit,
    	arma::vec init, 
	int maxit, 
	int maxstep, 
	int maxhs, 
	double lconv, 
	double gconv, 
	double xconv, 
	bool & isfirthconverge){
  int n = x.n_rows;
  int k = x.n_cols;
  arma::vec beta = init;
  int iter = 0;
  arma::vec pi_0 = -x * beta - offset; 
  pi_0 = arma::exp(pi_0) + 1;
  arma::vec pi = 1/pi_0;
  int evals = 1;
  arma::vec beta_old;
  arma::mat oneVec(k, 1 , arma::fill::ones);
  arma::mat XX_covs(k, k, arma::fill::zeros);
  isfirthconverge = false;
  while(iter <= maxit){
	beta_old = beta;
	arma::vec wpi = weight % pi % (1 - pi);
	arma::vec wpi_sqrt = arma::sqrt(wpi);
	arma::vec W2 = weight % wpi_sqrt;
	arma::mat XW2(n, k, arma::fill::zeros);
	for(int j = 0; j < k; j++){
                XW2.col(j) = x.col(j) % W2;
        }       

	arma::mat Q;
	arma::mat R;
	arma::qr_econ(Q, R, XW2);
	arma::vec h = Q % Q * oneVec;
	arma::vec U_star(2, arma::fill::zeros);
	arma::vec ypih;
	if(firth){
		ypih = (weight % (y - pi)) + (h % (0.5 - pi));
	}else{
		ypih = (weight % (y - pi));
	}
	//ypih.print();
	arma::vec xcol(n, arma::fill::zeros);
	U_star = x.t() * ypih;
	
	arma::mat XX_XW2(n, k, arma::fill::zeros);
	for(int j = 0; j < k; j++){
		xcol = x.col(j);
		XX_XW2.col(j) = xcol % wpi_sqrt; 	
	}
	arma::mat XX_Fisher = XX_XW2.t() * (XX_XW2);
	bool isinv = arma::inv_sympd (XX_covs, XX_Fisher); 
	if(!isinv){
		break;
	}	
	//}
	arma::vec delta = XX_covs * U_star;
	delta.replace(arma::datum::nan, 0);	

	double mx = arma::max(arma::abs(delta))/maxstep;
	if(mx > 1){
		delta = delta/mx;
	}
	evals = evals + 1;
	iter = iter + 1;
	beta = beta + delta;
	pi_0 = -x * beta - offset;
  	pi_0 = arma::exp(pi_0) + 1;
  	pi = 1/pi_0;
	if((iter == maxit) || ( (arma::max(arma::abs(delta)) <= xconv) & (abs(U_star).is_zero(gconv)))){
		isfirthconverge = true;
		break;
	}
  }
	arma::mat var;
	if(XX_covs.has_nan()){
		var = XX_covs;
		beta = arma::datum::nan;
	}
	return beta;
}

// [[Rcpp::export]]
void closeGenoFile(std::string & t_genoType)
{
  if(t_genoType == "bgen"){
    ptr_gBGENobj->closegenofile();
  }else if(t_genoType == "vcf"){
    ptr_gVCFobj->closegenofile();
  }else if(t_genoType == "plink"){
    ptr_gPLINKobj->closegenofile();
  }	  
}

// [[Rcpp::export]]
bool openOutfile(std::string t_traitType, bool isappend){
	bool isopen;
	if(!isappend){
	OutFile.open(g_outputFilePrefixGroup.c_str());
	isopen = OutFile.is_open();
	if(isopen){
		OutFile << "Region\tGroup\tmax_MAF\tPvalue_Burden\tBETA_Burden\tSE_Burden\t";
		if(ptr_gSAIGEobj->m_isCondition){
			OutFile << "Pvalue_Burden_c\tBeta_Burden_c\tseBeta_Burden_c\t";
		}
		OutFile << "MAC\t";
		if(t_traitType == "binary"){	
			OutFile << "MAC_case\tMAC_control\t";
		}
		OutFile << "Number_rare\tNumber_ultra_rare\n";
	}
     }else{
	OutFile.open(g_outputFilePrefixGroup.c_str(), std::ofstream::out | std::ofstream::app) ;				
	isopen = OutFile.is_open();
     }
     return(isopen);
}

// [[Rcpp::export]]
bool openOutfile_singleinGroup(std::string t_traitType, bool t_isImputation, bool isappend, bool t_isMoreOutput){
        bool isopen;
     if(!isappend){
        OutFile_singleInGroup.open(g_outputFilePrefixSingleInGroup.c_str());
	isopen = OutFile_singleInGroup.is_open();
        if(isopen){
                OutFile_singleInGroup << "CHR\tPOS\tMarkerID\tAllele1\tAllele2\tAC_Allele2\tAF_Allele2\t";
		if(t_isImputation){
			OutFile_singleInGroup << "imputationInfo\t";
		}else{
			
			OutFile_singleInGroup << "MissingRate\t";
		}
		OutFile_singleInGroup << "BETA\tSE\tTstat\tvar\tp.value\t";
	        if(t_traitType == "binary"){
                        OutFile_singleInGroup << "p.value.NA\tIs.SPA\t";
                }	

                if(ptr_gSAIGEobj->m_isCondition){
                        OutFile_singleInGroup << "BETA_c\tSE_c\tTstat_c\tvar_c\tp.value_c\t";
			if(t_traitType == "binary"){
				OutFile_singleInGroup << "p.value.NA_c\t";		
			}
                }
		
	        if(t_traitType == "binary"){
                        OutFile_singleInGroup << "AF_case\tAF_ctrl\tN_case\tN_ctrl";
 			if(t_isMoreOutput){
                                OutFile_singleInGroup << "\tN_case_hom\tN_case_het\tN_ctrl_hom\tN_ctrl_het";
                        }
			OutFile_singleInGroup << "\n";
                }else if(t_traitType == "quantitative"){
			OutFile_singleInGroup << "N\n";	
			
		}	

        }
    }else{
       OutFile_singleInGroup.open(g_outputFilePrefixSingleInGroup.c_str(), std::ofstream::out | std::ofstream::app);
       isopen = OutFile_singleInGroup.is_open();
    }

    //if(t_testType != "BURDEN" && t_isfastTest){
    //  OutFile_singleInGroup_temp.open(g_outputFilePrefixSingleInGroup_temp.c_str(), std::ofstream::out | std::ofstream::app);
    //}
      return(isopen);
}


// [[Rcpp::export]]
bool openOutfile_single(std::string t_traitType, bool t_isImputation, bool isappend, bool t_isMoreOutput){
      bool isopen;
      if(!isappend){
        OutFile_single.open(g_outputFilePrefixSingle.c_str());
        isopen = OutFile_single.is_open();
        if(isopen){
                OutFile_single << "CHR\tPOS\tMarkerID\tAllele1\tAllele2\tAC_Allele2\tAF_Allele2\t";
                if(t_isImputation){
                        OutFile_single << "imputationInfo\t";
                }else{

                        OutFile_single << "MissingRate\t";
                }
                OutFile_single << "BETA\tSE\tTstat\tvar\tp.value\t";
                if(t_traitType == "binary"){
                        OutFile_single << "p.value.NA\tIs.SPA\t";
                }

                if(ptr_gSAIGEobj->m_isCondition){
                        OutFile_single << "BETA_c\tSE_c\tTstat_c\tvar_c\tp.value_c\t";
			if(t_traitType == "binary"){
				OutFile_single << "p.value.NA_c\t";
			}
                }
		

                if(t_traitType == "binary"){
                        OutFile_single << "AF_case\tAF_ctrl\tN_case\tN_ctrl";

			if(t_isMoreOutput){	
				OutFile_single << "\tN_case_hom\tN_case_het\tN_ctrl_hom\tN_ctrl_het";
			}
			OutFile_single << "\n";

                }else if(t_traitType == "quantitative"){
                        OutFile_single << "N\n";

                }

        }
     }else{
        OutFile_single.open(g_outputFilePrefixSingle.c_str(), std::ofstream::out | std::ofstream::app);
        isopen = OutFile_single.is_open();
     }

     return(isopen);
}


void writeOutfile_single(bool t_isMoreOutput,
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
                        std::vector<uint32_t> & N_Vec

){
  int numtest = 0;
  for(unsigned int k = 0; k < pvalVec.size(); k++){
        if(pvalVec.at(k) != "NA"){
                numtest = numtest + 1;
                OutFile_single << chrVec.at(k);
                OutFile_single << "\t";
                OutFile_single << posVec.at(k);
                OutFile_single << "\t";
                OutFile_single << markerVec.at(k);
                OutFile_single << "\t";
                OutFile_single << refVec.at(k);
                OutFile_single << "\t";
                OutFile_single << altVec.at(k);
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

                if(t_traitType == "binary"){
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
			if(t_traitType == "binary"){	
                        	OutFile_single << pvalNA_cVec.at(k);
                        	OutFile_single << "\t";
			}	
                }
                if(t_traitType == "binary"){
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
                        OutFile_single << "\n";
                }else if(t_traitType == "quantitative"){
                        OutFile_single << N_Vec.at(k);
                        OutFile_single << "\n";

                }
        }
  }
  std::cout << numtest << " markers were tested." << std::endl;
  if(t_traitType == "binary"){
      if(t_isFirth){
        std::cout << "Firth approx was applied to " << mFirth << " markers. " << mFirthConverge << " sucessfully converged." <<std::endl;
      }
   }
}


// [[Rcpp::export]]
void set_flagSparseGRM_cur_SAIGE(bool t_flagSparseGRM_cur){
	ptr_gSAIGEobj->set_flagSparseGRM_cur(t_flagSparseGRM_cur);
}

// [[Rcpp::export]]
void set_flagSparseGRM_cur_SAIGE_org(){
	ptr_gSAIGEobj->set_flagSparseGRM_cur(ptr_gSAIGEobj->m_flagSparseGRM);
}



void set_varianceRatio(double MAC, bool isSingleVarianceRatio){
    bool hasVarRatio;
    if(!ptr_gSAIGEobj->m_isFastTest){
       ptr_gSAIGEobj->set_flagSparseGRM_cur(ptr_gSAIGEobj->m_flagSparseGRM);
       if(!isSingleVarianceRatio){
         hasVarRatio = ptr_gSAIGEobj->assignVarianceRatio(MAC, ptr_gSAIGEobj->m_flagSparseGRM);
       }else{
         ptr_gSAIGEobj->assignSingleVarianceRatio(ptr_gSAIGEobj->m_flagSparseGRM);
       }
     }else{
       ptr_gSAIGEobj->set_flagSparseGRM_cur(false);
       if(!isSingleVarianceRatio){
         hasVarRatio = ptr_gSAIGEobj->assignVarianceRatio(MAC, false);
       }else{
         ptr_gSAIGEobj->assignSingleVarianceRatio(false);
       }
     }
}




void writeOutfile_BURDEN(std::string regionName,
			std::vector<std::string>  & BURDEN_AnnoName_Vec,
			std::vector<std::string> & BURDEN_maxMAFName_Vec,
			std::vector<std::string> & BURDEN_pval_Vec,
			std::vector<double> & BURDEN_Beta_Vec,
			std::vector<double> & BURDEN_seBeta_Vec,
			std::vector<std::string> & BURDEN_pval_cVec,
			std::vector<double> & BURDEN_Beta_cVec,
			std::vector<double> & BURDEN_seBeta_cVec,
			arma::vec & MAC_GroupVec,
			arma::vec & MACCase_GroupVec,
			arma::vec & MACControl_GroupVec,
			arma::vec & NumRare_GroupVec,
			arma::vec & NumUltraRare_GroupVec,
			double cctpval,
			double cctpval_cond,
			unsigned int q_anno,
			unsigned int q_maf,
			bool isCondition,
			std::string t_traitType){
     unsigned int i;
     for(unsigned int j = 0; j < q_anno; j++){
       for(unsigned int m = 0; m < q_maf; m++){
           i = j*q_maf+m;
	   if(BURDEN_pval_Vec.at(i) != "NA"){
           OutFile << regionName;
           OutFile << "\t";
           OutFile << BURDEN_AnnoName_Vec.at(i);
           OutFile << "\t";
           OutFile << BURDEN_maxMAFName_Vec.at(i);
           OutFile << "\t";
           OutFile << BURDEN_pval_Vec.at(i);
           OutFile << "\t";
           OutFile << BURDEN_Beta_Vec.at(i);
           OutFile << "\t";
           OutFile << BURDEN_seBeta_Vec.at(i);
           OutFile << "\t";
           if(isCondition){
               OutFile << BURDEN_pval_cVec.at(i);
               OutFile << "\t";
               OutFile << BURDEN_Beta_cVec.at(i);
               OutFile << "\t";
               OutFile << BURDEN_seBeta_cVec.at(i);
               OutFile << "\t";
	   }
	   OutFile << MAC_GroupVec(i);
           OutFile << "\t";
           if(t_traitType == "binary"){
               OutFile << MACCase_GroupVec(i);
               OutFile << "\t";
               OutFile << MACControl_GroupVec(i);
               OutFile << "\t";
           }
           OutFile << NumRare_GroupVec(i);
           OutFile << "\t";
           OutFile << NumUltraRare_GroupVec(i);
           OutFile << "\n";
	}
     }
   }  
     OutFile << regionName;
     OutFile << "\tCauchy\tNA\t";
     OutFile << cctpval;
     OutFile << "\tNA\tNA\t";	
     if(isCondition){
	OutFile << cctpval_cond;
	OutFile << "\tNA\tNA\t";
     }
     OutFile << "NA\t";
     if(t_traitType == "binary"){
        OutFile << "NA\t";
        OutFile << "NA\t";
     }
     OutFile << "NA\t";
     OutFile << "NA\n";	
}


int writeOutfile_singleInGroup(bool t_isMoreOutput,
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
                        std::ofstream & t_OutFile_singleInGroup){
  int numofUR = 0;
  for(unsigned int k = 0; k < pvalVec.size(); k++){
        //if(std::isfinite(pvalVec.at(k))){
        //if(!std::isnan(pvalVec.at(k))){
        if(pvalVec.at(k) != "NA"){
		if(chrVec.at(k) == "UR"){
                        numofUR = numofUR + 1;
                }
                t_OutFile_singleInGroup << chrVec.at(k);
                t_OutFile_singleInGroup << "\t";
                t_OutFile_singleInGroup << posVec.at(k);
                t_OutFile_singleInGroup << "\t";
                t_OutFile_singleInGroup << markerVec.at(k);
                t_OutFile_singleInGroup << "\t";
                t_OutFile_singleInGroup << refVec.at(k);
                t_OutFile_singleInGroup << "\t";
                t_OutFile_singleInGroup << altVec.at(k);
                t_OutFile_singleInGroup << "\t";
                t_OutFile_singleInGroup << altCountsVec.at(k);
                t_OutFile_singleInGroup << "\t";
                t_OutFile_singleInGroup << altFreqVec.at(k);
                t_OutFile_singleInGroup << "\t";

                if(t_isImputation){
                        t_OutFile_singleInGroup << imputationInfoVec.at(k);
                        t_OutFile_singleInGroup << "\t";
                }else{
                        t_OutFile_singleInGroup << missingRateVec.at(k);
                        t_OutFile_singleInGroup << "\t";

                }
                t_OutFile_singleInGroup << BetaVec.at(k);
                t_OutFile_singleInGroup << "\t";
                t_OutFile_singleInGroup << seBetaVec.at(k);
                t_OutFile_singleInGroup << "\t";
                t_OutFile_singleInGroup << TstatVec.at(k);
                t_OutFile_singleInGroup << "\t";
                t_OutFile_singleInGroup << varTVec.at(k);
                t_OutFile_singleInGroup << "\t";
                t_OutFile_singleInGroup << pvalVec.at(k);
                t_OutFile_singleInGroup << "\t";

                if(t_traitType == "binary"){
                        t_OutFile_singleInGroup << pvalNAVec.at(k);
                        t_OutFile_singleInGroup << "\t";
                        t_OutFile_singleInGroup << std::boolalpha << isSPAConvergeVec.at(k);
                        t_OutFile_singleInGroup << "\t";
                }
                if(t_isCondition){
                        t_OutFile_singleInGroup << Beta_cVec.at(k);
                        t_OutFile_singleInGroup << "\t";
                        t_OutFile_singleInGroup << seBeta_cVec.at(k);
                        t_OutFile_singleInGroup << "\t";
                        t_OutFile_singleInGroup << Tstat_cVec.at(k);
                        t_OutFile_singleInGroup << "\t";
                        t_OutFile_singleInGroup << varT_cVec.at(k);
                        t_OutFile_singleInGroup << "\t";
                        t_OutFile_singleInGroup << pval_cVec.at(k);
                        t_OutFile_singleInGroup << "\t";
			if(t_traitType == "binary"){
                        	t_OutFile_singleInGroup << pvalNA_cVec.at(k);
                        	t_OutFile_singleInGroup << "\t";
			}
                }
                if(t_traitType == "binary"){
                        t_OutFile_singleInGroup << AF_caseVec.at(k);
                        t_OutFile_singleInGroup << "\t";
                        t_OutFile_singleInGroup << AF_ctrlVec.at(k);
                        t_OutFile_singleInGroup << "\t";
                        t_OutFile_singleInGroup << N_caseVec.at(k);
                        t_OutFile_singleInGroup << "\t";
                        t_OutFile_singleInGroup << N_ctrlVec.at(k);

                        if(t_isMoreOutput){
                                t_OutFile_singleInGroup << "\t";
                                t_OutFile_singleInGroup << N_case_homVec.at(k);
                                t_OutFile_singleInGroup << "\t";
                                t_OutFile_singleInGroup << N_case_hetVec.at(k);
                                t_OutFile_singleInGroup << "\t";
                                t_OutFile_singleInGroup << N_ctrl_homVec.at(k);
                                t_OutFile_singleInGroup << "\t";
                                t_OutFile_singleInGroup << N_ctrl_hetVec.at(k);
                        }
                        t_OutFile_singleInGroup << "\n";
                }else if(t_traitType == "quantitative"){
                        t_OutFile_singleInGroup << N_Vec.at(k);
                        t_OutFile_singleInGroup << "\n";

                }
        }
  }
  return(numofUR);

}

// [[Rcpp::export]]
void copy_singleInGroup(){
  std::ifstream ini_file;
  ini_file.open(g_outputFilePrefixSingleInGroup_temp.c_str());
  if (!ini_file)
  {
    std::cout << "Error in Opening the temp file!" << std::endl;
  }
  std::string str;
  while (getline(ini_file, str))
  {
    OutFile_singleInGroup << str << "\n"; 	
  }
  ini_file.close();
}


