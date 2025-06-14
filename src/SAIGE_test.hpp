
#ifndef SAIGE_HPP
#define SAIGE_HPP

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>


namespace SAIGE{

class SAIGEClass
{
    private:
      arma::mat m_XVX;
      arma::mat m_XVX_inv_XV;
      arma::mat m_X;
      arma::mat m_Sigma_iXXSigma_iX;
      arma::vec m_res;
      arma::vec m_resout;
      arma::vec m_mu;
      arma::vec m_mu2;
      arma::vec m_tauvec;
      arma::vec  m_S_a;
      std::string m_traitType; 
      std::string m_impute_method;
      std::vector<arma::vec> m_res_mt_vec;	
      //std::vector<arma::mat> m_XV_mt_vec;
    public:
        // Default constructor
       SAIGEClass() {
            // Initialize m_res_mt_vec as an empty vector
            m_res_mt_vec = std::vector<arma::vec>();
	    //m_XV_mt_vec = std::vector<arma::mat>();
        }	


      std::vector<uint32_t> m_condition_genoIndex;
      arma::mat m_XXVX_inv;
      arma::mat m_XV;
      int m_n, m_p; //MAIN Dimensions: sample size, number of covariates
      double m_varRatioVal;
      arma::vec m_varRatio_sparse;
      arma::vec m_varRatio_null;
      arma::vec m_y;

      std::vector<std::string> m_traitType_vec;

      unsigned int m_itrait, m_startip, m_startin, m_startic, m_endip, m_endin, m_endic;
      arma::uvec m_ip;


      bool m_isOutputAFinCaseCtrl;
      bool m_isOutputNinCaseCtrl;
      bool m_isOutputHetHomCountsinCaseCtrl;
      arma::uvec m_case_indices;
      arma::uvec m_ctrl_indices;
      arma::uvec m_case_hom_indices;
      arma::uvec m_case_het_indices;
      arma::uvec m_ctrl_hom_indices;
      arma::uvec m_ctrl_het_indices;
      int m_n_case;
      int m_n_ctrl;
      arma::sp_mat m_SigmaMat_sp;
      bool m_flagSparseGRM;
      bool m_flagSparseGRM_cur;
      bool m_isFastTest;
      double m_pval_cutoff_for_fastTest; 
      double m_SPA_Cutoff;
      arma::umat m_locationMat;
      arma::vec m_valueVec;
      int m_dimNum;	
      arma::vec m_cateVarRatioMinMACVecExclude; 
      arma::vec m_cateVarRatioMaxMACVecInclude;
      arma::mat m_P2Mat_cond;
      int m_numMarker_cond;
      arma::mat m_VarInvMat_cond;
      arma::mat m_VarMat_cond;
      arma::vec m_Tstat_cond;
      //arma::vec m_G2_Weight_cond;
      arma::mat m_G2_Weight_cond;
      arma::vec m_MAF_cond;
      //arma::mat m_qsum_cond_Mat;
      arma::mat  m_qsum_cond;
      arma::vec m_gsum_cond;
      arma::mat m_gsum_cond_Mat;
      std::vector<std::string> m_p_cond;
      arma::vec m_scalefactor_G2_cond;
      arma::mat m_scalefactor_G2_cond_Mat;
      arma::mat m_VarInvMat_cond_scaled_weighted;
      //arma::mat m_VarInvMat_cond_region_binary;
      bool m_isCondition;
      bool m_is_Firth_beta;
      double m_pCutoffforFirth;
     arma::vec  m_offset;	
      bool m_isVarPsadj;
      bool m_islog10p;




arma::mat m_XVX_mt;
arma::mat m_XV_mt;
arma::mat m_XXVX_inv_mt;
arma::mat m_XVX_inv_XV_mt;
arma::mat m_Sigma_iXXSigma_iX_mt;
arma::mat m_X_mt;
arma::mat m_S_a_mt;
arma::mat m_res_mt;
arma::mat m_resout_mt;
arma::mat m_mu2_mt;
arma::mat m_mu_mt;
arma::mat m_varRatio_sparse_mt;
arma::mat m_varRatio_null_mt;
arma::mat    m_varRatio_null_sample_mt;
arma::mat    m_varRatio_null_noXadj_mt;


arma::mat m_tauvec_mt;
arma::mat m_y_mt;
arma::umat m_locationMat_mt;
arma::mat m_valueVec_mt;
arma::mat m_offset_mt;
arma::umat m_sampleindices_mt;

arma::uvec m_sampleindices_vec;
arma::uvec m_sampleIndexLenVec;
arma::uvec m_colXvec;
  ////////////////////// -------------------- functions ---------------------------------- //////////////////////
SAIGEClass(
        arma::mat & t_XVX_mt,
        arma::mat & t_XXVX_inv_mt,
        arma::mat & t_XV_mt,
        arma::mat & t_XVX_inv_XV_mt,
        arma::mat & t_Sigma_iXXSigma_iX_mt,
        arma::mat & t_X_mt,
        arma::mat & t_S_a_mt,
        arma::mat & t_res_mt,
        arma::mat & t_mu2_mt,
        arma::mat & t_mu_mt,
        arma::mat & t_varRatio_sparse_mt,
        arma::mat & t_varRatio_null_mt,

        arma::vec & t_cateVarRatioMinMACVecExclude,
        arma::vec & t_cateVarRatioMaxMACVecInclude,
        double t_SPA_Cutoff,
        arma::mat & t_tauvec_mt,
	std::vector<std::string> & t_traitType,
        arma::mat & t_y_mt,
        std::string t_impute_method,
        bool t_flagSparseGRM,
        bool t_isFastTest,
        double t_pval_cutoff_for_fastTest,
        arma::umat & t_locationMat_mt,
        arma::mat & t_valueVec_mt,

        int t_dimNum, //to update

        bool t_isCondition,
        std::vector<uint32_t> & t_condition_genoIndex,
        bool t_is_Firth_beta,
        double t_pCutoffforFirth,
        arma::mat & t_offset_mt,
        arma::mat & t_resout_mt,
        arma::uvec & t_colXvec,
        arma::uvec & t_sampleIndexLenVec,
	arma::umat & t_sampleIndexMat); 

   void set_seed(unsigned int seed);

   void scoreTest(arma::vec & t_GVec,
                     double& t_Beta,
                     double& t_seBeta,
                     std::string& t_pval_str,
		      double& t_pval,
                     bool& t_islogp,
                     double t_altFreq,
                     double &t_Tstat,
                     double &t_var1,
                     double &t_var2,
                     arma::vec & t_gtilde,
                     arma::vec & t_P2Vec,
		     double& t_gy,
                     bool t_is_region,
		     arma::uvec & t_indexForNonZero);


void scoreTestFast(arma::vec & t_GVec,
                     arma::uvec & t_indexForNonZero,
                     double& t_Beta,
                     double& t_seBeta,
                     std::string& t_pval_str,
                     double& t_pval,
                     bool& t_islogp,
                     double t_altFreq,
                     double &t_Tstat,
                     double &t_var1,
                     double &t_var2);


     void set_flagSparseGRM_cur(bool t_flagSparseGRM_cur);

     void get_mu(arma::vec & t_mu);

     void getadjG(arma::vec & t_GVec, arma::vec & g);
     void getadjGFast(arma::vec & t_GVec, arma::vec & g,  arma::uvec & iIndex);

     void getMarkerPval(arma::vec & t_GVec,
			        arma::uvec & iIndex,
                               arma::uvec & iIndexComVec,
                               double& t_Beta,
                               double& t_seBeta,
                               std::string& t_pval,
                               std::string& t_pval_noSPA,
                               double t_altFreq,
                               double& t_Tstat,
				double& t_gy,
                               double& t_var1,
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
                                arma::rowvec & t_G1tilde_P_G2tilde,
				bool & t_isFirth,
                                bool & t_isFirthConverge, 
				bool t_isER);


    void getindices(arma::uvec & t_case_indices,
      arma::uvec & t_ctrl_indices);


    void setupSparseMat(int r, arma::umat & locationMatinR, arma::vec & valueVecinR);

    arma::sp_mat gen_sp_SigmaMat();

    bool assignVarianceRatio(double MAC, bool issparseforVR);

    void assignSingleVarianceRatio(bool issparseforVR);


    void assignSingleVarianceRatio_withinput(double t_varRatioVal);

/*
    void assignConditionFactors(
      arma::mat & t_P2Mat_cond,
      arma::mat & t_VarInvMat_cond,
      arma::mat & t_VarMat_cond,
      arma::vec & t_Tstat_cond,
      arma::vec & t_G2_Weight_cond,
      arma::vec & t_MAF_cond,
      double t_qsum_cond,
      arma::vec & t_gsum_cond,
      std::vector<std::string> & t_p_cond);

    void assignConditionFactors(
      arma::mat & t_P2Mat_cond,
      arma::mat & t_VarInvMat_cond,
      arma::mat & t_VarMat_cond,
      arma::vec & t_Tstat_cond,
      arma::mat & t_G2_Weight_cond,
      arma::vec & t_MAF_cond,
      arma::mat & t_qsum_cond,
      arma::mat & t_gsum_cond,
      std::vector<std::string> & t_p_cond
      );
*/

     void assignConditionFactors_scalefactor(
        arma::vec & t_scalefactor_G2_cond);	


    void extract_XV_XXVX_inv(arma::mat & t_XV, arma::mat & t_XXVX_inv);

    void fast_logistf_fit_simple(arma::mat & x,
                arma::vec & y,
                arma::vec & offset,
                bool firth,
        arma::vec init,
        int maxit,
        int maxstep,
        int maxhs,
        double lconv,
        double gconv,
        double xconv,
        double & beta_G,
        double & sebeta_G,
	bool & isfirthconverge);	
/*
    void assignConditionFactors(
      arma::mat & t_P2Mat_cond,
      arma::mat & t_VarInvMat_cond,
            arma::mat & t_VarMat_cond,
      arma::vec & t_Tstat_cond,
       arma::mat & t_G2_Weight_cond,
      arma::vec & t_MAF_cond,
        arma::mat & t_qsum_cond,
        arma::mat & t_gsum_cond,
      arma::vec & t_p_cond);
*/


  void assignConditionFactors(
      arma::mat & t_P2Mat_cond,
      arma::mat & t_VarInvMat_cond,
      arma::mat & t_VarMat_cond,
      arma::vec & t_Tstat_cond,
      arma::mat & t_G2_Weight_cond,
      arma::vec & t_MAF_cond,
      arma::mat & t_qsum_cond,
      arma::mat & t_gsum_cond,
      std::vector<std::string> & t_p_cond
      );




     void assignConditionFactors_scalefactor_multiTrait( arma::mat & t_scalefactor_G2_cond,
                  unsigned int oml);
void assign_for_itrait_binaryindices(unsigned int t_itrait);

void assign_for_itrait(unsigned int t_itrait);

void assignSingleVarianceRatio(bool issparseforVR, bool isnoXadj, bool issample);

void assignSingleVarianceRatio(bool issparseforVR, bool isnoXadj);

    bool assignVarianceRatio(double MAC, bool issparseforVR, bool isnoXadj);
    bool assignVarianceRatio(double MAC, bool issparseforVR, bool isnoXadj, bool issample);




};


}
#endif
