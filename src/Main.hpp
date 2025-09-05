#ifndef MAIN_HPP
#define MAIN_HPP

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

void setAssocTest_GlobalVarsInCPP(std::string t_impute_method,
                double t_missing_cutoff,
                               double t_min_maf_marker,
                               double t_min_mac_marker,
                               double t_min_info_marker,
                               double t_dosage_zerod_cutoff,
                               double t_dosage_zerod_MAC_cutoff,
                               arma::vec & t_weights_beta,
                               std::string t_outputFilePrefix,
                               double t_MACCutoffforER,
                               bool t_isadmixed);


void setAssocTest_GlobalVarsInCPP_indexInModel_male(arma::uvec & t_indexInModel_male);

void setAssocTest_GlobalVarsInCPP_X_PARregion_mat(arma::umat & t_X_PARregion_mat);

void processMale_XnonPAR(arma::vec & t_GVec,  uint32_t& t_pd , arma::umat & t_XPARregion);


void setMarker_GlobalVarsInCPP(
                               bool t_isOutputMoreDetails,
                               int t_marker_chunksize
                               );


void setRegion_GlobalVarsInCPP(
                               arma::vec t_max_maf_region,
                               unsigned int t_max_markers_region,
                               double t_MACCutoff_to_CollapseUltraRare,
			       double t_min_gourpmac_for_burdenonly,
			        arma::vec t_r_corr);



void mainMarkerInCPP(
                           std::string & t_genoType,     // "PLINK", "BGEN"
                           std::string & t_traitType,
			   std::vector<std::string> & t_genoIndex_prev,
                           std::vector<std::string>  & t_genoIndex,
                           bool & t_isMoreOutput,
                           bool & t_isImputation,
			   bool & t_isFirth);

bool Unified_getOneMarker(std::string & t_genoType,   // "PLINK", "BGEN"
				uint64_t & t_gIndex_prev,
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
                               std::vector<uint>& t_indexForNonZero,
			       arma::vec & t_GVec,
			       bool t_isImputation);

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
			   bool t_isER);


Rcpp::List mainRegionInCPP(
                           std::string t_genoType,     // "PLINK", "BGEN"
			   std::vector<std::string> & t_genoIndex_prev,
                           std::vector<std::string> & t_genoIndex,
                           arma::mat & annoIndicatorMat,
                           arma::vec & maxMAFVec,
                           std::string t_outputFile,
                           std::string t_traitType,
                           unsigned int t_n,           // sample size
                           arma::mat P1Mat,            // edited on 2021-08-19: to avoid repeated memory allocation of P1Mat and P2Mat
                           arma::mat P2Mat,
                           std::string t_regionTestType,
			   bool t_isImputation,
			   arma::vec & t_weight,
			   arma::vec & t_weight_cond,
			   bool t_isSingleinGroupTest,
			   bool t_isOutputMarkerList,
			   std::vector<std::string> & annoStringVec,
                           std::string regionName,
			   bool t_isFastTest,
			   bool t_isMoreOutput);



void setPLINKobjInCPP(std::string t_bimFile,
                      std::string t_famFile,
                      std::string t_bedFile,
                      std::vector<std::string> & t_SampleInModel,
                      std::string t_AlleleOrder);



void setBGENobjInCPP(std::string t_bgenFileName,
                     std::string t_bgenFileIndex,
                     std::vector<std::string> & t_SampleInBgen,
                     std::vector<std::string> & t_SampleInModel,
                     std::string t_AlleleOrder);

void setVCFobjInCPP(std::string t_vcfFileName,
            std::string t_vcfFileIndex,
            std::string t_vcfField,
            std::vector<std::string> & t_SampleInModel);



void setSAIGEobjInCPP(arma::mat & t_XVX,
        arma::mat & t_XXVX_inv,
        arma::mat & t_XV,
        arma::mat & t_XVX_inv_XV,
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
	arma::vec & t_offset);

void assign_conditionMarkers_factors(
                           std::string t_genoType,     // "plink", "bgen", "vcf"
                           std::vector<std::string> & t_genoIndex_prev,
                           std::vector<std::string> & t_genoIndex,
                           unsigned int t_n,
                           arma::vec & t_weight_cond
                           );


void assign_conditionMarkers_factors_binary_region(
                           arma::vec & scalefactor_G2_cond);

void set_iterator_inVcf(std::string & variantList);

void set_iterator_inVcf(std::string & variantList, std::string & chrom, int & beg_pd, int & end_pd);

bool check_Vcf_end();

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
        double xconv);

Rcpp::List RegionSetUpConditional_binary_InCPP(arma::vec & t_weight_cond);

void closeGenoFile(std::string & t_genoType);

bool openOutfile(std::string t_traitType, bool isappend);


bool openOutfile_singleinGroup(std::string t_traitType, bool t_isImputation, bool isappend, bool t_isMoreOutput);

bool openOutfile_single(std::string t_traitType, bool t_isImputation, bool isappend, bool t_isMoreOutput);

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
                        std::vector<uint32_t> & N_Vec);




int writeOutfile_singleinGroup(bool t_isMoreOutput,
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
                        std::vector<uint32_t> & N_Vec);


void set_flagSparseGRM_cur_SAIGE(bool t_flagSparseGRM_cur);

void set_flagSparseGRM_cur_SAIGE_org();


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
                        std::string cctpval,
                        std::string cctpval_cond,
                        unsigned int q_anno,
                        unsigned int q_maf,
                        bool isCondition,
                        std::string t_traitType);

void copy_singleInGroup();

void set_varianceRatio(double MAC, bool isSingleVarianceRatio);

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
                        std::ofstream & t_OutFile_singleInGroup);

uint32_t Unified_getSampleSizeinGeno(std::string & t_genoType);
uint32_t Unified_getSampleSizeinAnalysis(std::string & t_genoType);

int writeOutfile_singleInadmixed(bool t_isMoreOutput,
                        bool t_isImputation,
                        bool t_isCondition,
                        bool t_isFirth,
                        int mFirth,
                        int mFirthConverge,
			int q0,
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
			std::vector<uint32_t> & nonNAIndexVec,
                        std::ofstream & t_OutFile_singleInGroup,
			std::string markerName);

bool openOutfile_single_admixed(std::string t_traitType, bool t_isCondition, bool t_isMoreOutput, std::vector<std::string> & pvalVec, std::ofstream & t_OutFile_singleInGroup);

void mainAdmixedInCPP_inner(
                           std::string t_genoType,     // "PLINK", "BGEN"
                           std::vector<std::string> & t_genoIndex_prev,
                           std::vector<std::string> & t_genoIndex,
                           std::string t_outputFile,
                           std::string t_traitType,
                           unsigned int t_n,           // sample size
                           std::string t_regionTestType,
                           bool t_isImputation,
                           arma::vec & t_weight,
                           arma::vec & t_weight_cond,
                           std::string regionName,
                           bool t_isFastTest,
                           bool t_isMoreOutput,
			   bool t_isWriteHeader
);

void mainAdmixedInCPP(
        Rcpp::List & RegionList,
        std::string t_genoType,
        std::string t_outputFile,
        std::string t_traitType,
        unsigned int t_n,           // sample size
        std::string t_regionTestType,
        arma::vec & t_weight_cond,
        bool t_isImputation,
        bool t_isFastTest,
        bool t_isMoreOutput,
	bool t_isWriteHeader);


void mainMarkerAdmixedInCPP(
                           std::string & t_genoType,     // "PLINK", "BGEN"
                           std::string & t_traitType,
                           std::vector<std::string> & t_genoIndex_prev,
                           std::vector<std::string> & t_genoIndex,
                           bool & t_isMoreOutput,
                           bool & t_isImputation,
                           bool & t_isFirth,
                           int t_NumberofANC,
			                              double t_pvalcutoff_of_haplotype);

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
  int t_NumberofANC,
    std::vector<std::string> & pvalHap_Vec
);

bool openOutfile_single_admixed_new(std::string t_traitType, bool t_isImputation, bool isappend, bool t_isMoreOutput, int t_NumberofANC);


void assign_conditionHaplotypes(
			   std::string t_traitType,
                           std::string t_genoType,     //"vcf"
                           unsigned int t_n,
                           int t_NumberofANC,
                           uint64_t gIndex,
                           uint64_t gIndex_prev,
                           arma::vec & nanc_case_vec,
                           arma::vec & nanc_ctrl_vec,
                           arma::vec & nanc_vec,
                           arma::uvec & not_nan_anc_indices_vec,
			                              double t_pvalcutoff_of_haplotype,
						      bool & isconditiononHaplo
                           );



#endif
