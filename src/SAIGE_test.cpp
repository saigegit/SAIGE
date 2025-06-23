#define ARMA_USE_SUPERLU 1

// [[Rcpp::depends(BH)]]

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>



#include "SAIGE_test.hpp"
#include "SPA.hpp"
#include "ER_binary_func.hpp"
#include "UTIL.hpp"
#include "getMem.hpp"
#include "getMem.hpp"
#include <thread>         // std::this_thread::sleep_for
#include <chrono>         // std::chrono::seconds
#include <boost/iostreams/filter/zstd.hpp>
#include <boost/date_time.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/chi_squared.hpp>

namespace SAIGE {

SAIGEClass::SAIGEClass(
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
        arma::mat & t_varRatio_null_noXadj_mt,

        arma::vec & t_cateVarRatioMinMACVecExclude,
        arma::vec & t_cateVarRatioMaxMACVecInclude,
        double t_SPA_Cutoff,
        arma::mat & t_tauvec_mt,
	std::vector<std::string> & t_traitType,
        //std::string t_traitType,
        arma::mat & t_y_mt,
        std::string t_impute_method,
        bool t_flagSparseGRM,
        bool t_isFastTest,
	bool t_isnoadjCov,
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
	arma::umat & t_sampleIndexMat) {

    std::cout << "SAIGEClass" << std::endl;
    //m_res_mt_vec = std::vector<arma::vec>();
    //m_XV_mt_vec = std::vector<arma::mat>();
    
    
    m_XVX_mt = t_XVX_mt;
    m_XV_mt = t_XV_mt;
    m_XXVX_inv_mt = t_XXVX_inv_mt;
    m_XVX_inv_XV_mt = t_XVX_inv_XV_mt;
    m_Sigma_iXXSigma_iX_mt = t_Sigma_iXXSigma_iX_mt;
    m_isVarPsadj = false;
    if (m_Sigma_iXXSigma_iX_mt.n_cols == 1 && m_Sigma_iXXSigma_iX_mt.n_rows == 1) {
        m_isVarPsadj = false;
    } else {
        m_isVarPsadj = true;
    }
    m_X_mt = t_X_mt;
    m_S_a_mt = t_S_a_mt;
    m_res_mt = t_res_mt;
    m_resout_mt = t_resout_mt;
    m_mu2_mt = t_mu2_mt;
    m_mu_mt = t_mu_mt;
    m_varRatio_sparse_mt = t_varRatio_sparse_mt;
    m_varRatio_null_mt = t_varRatio_null_mt;
    m_varRatio_null_noXadj_mt = t_varRatio_null_noXadj_mt;
    
    m_varRatio_null_mt.print("m_varRatio_null_mt");
    m_varRatio_sparse_mt.print("m_varRatio_sparse_mt");
    m_varRatio_null_noXadj_mt.print("m_varRatio_null_noXadj_mt");
    
    m_cateVarRatioMinMACVecExclude = t_cateVarRatioMinMACVecExclude;
    m_cateVarRatioMaxMACVecInclude = t_cateVarRatioMaxMACVecInclude;
    m_tauvec_mt = t_tauvec_mt;
    m_traitType_vec = t_traitType;
    m_y_mt = t_y_mt;
    
    m_n = t_y_mt.n_rows;
    m_p = t_XV_mt.n_rows;
    m_SPA_Cutoff = t_SPA_Cutoff;
    m_impute_method = t_impute_method;
    m_isCondition = t_isCondition;
    m_condition_genoIndex = t_condition_genoIndex;
    if (m_isCondition) {
        m_numMarker_cond = t_condition_genoIndex.size();
    } else {
        m_numMarker_cond = 0;
    }
    m_colXvec = t_colXvec;
    m_sampleIndexLenVec = t_sampleIndexLenVec;
    m_sampleindices_mt = t_sampleIndexMat-1;

    bool isbinary=false;
    for(unsigned int j = 0; j < m_traitType_vec.size(); j++){
	if(m_traitType_vec.at(j) == "binary"){
		isbinary=true;	
	}
	if(isbinary){
                m_is_Firth_beta = t_is_Firth_beta;
                m_pCutoffforFirth = t_pCutoffforFirth;
                m_offset_mt = t_offset_mt;
	}
    }
    m_dimNum = t_dimNum;
    m_flagSparseGRM = t_flagSparseGRM;
    m_isFastTest = t_isFastTest;
    m_isnoadjCov = t_isnoadjCov;
    m_pval_cutoff_for_fastTest = t_pval_cutoff_for_fastTest;
    if (m_dimNum != 0) {
        m_locationMat_mt = t_locationMat_mt;
        m_valueVec_mt = t_valueVec_mt;
    }

        m_res_sum_vec.set_size(m_traitType_vec.size());
	    m_mu2_sum_vec.set_size(m_traitType_vec.size());

    unsigned int k_startip, k_endip, k_p;
    arma::uvec k_ip;
    for(unsigned int k_itrait = 0; k_itrait < m_traitType_vec.size(); k_itrait++){
        arma::uvec sampleindices_sub_vec = m_sampleindices_mt.col(k_itrait);
	arma::uvec k_sampleindices_vec = sampleindices_sub_vec.subvec(0, (m_sampleIndexLenVec[k_itrait]-1));
	if(k_itrait == 0){
                k_startip = 0;
        }else{
                k_startip = arma::sum(m_colXvec.subvec(0, k_itrait - 1));
        }
	k_p = m_colXvec(k_itrait);
	k_endip =  k_startip + m_colXvec[k_itrait]-1;
	k_ip.set_size(k_endip - k_startip + 1);
	for (unsigned int i = 0; i < k_ip.size(); ++i) {
		k_ip(i) = k_startip + i;
	}
	arma::vec  k_y_sub = m_y_mt.col(k_itrait);
	arma::vec k_y = k_y_sub.elem(k_sampleindices_vec);
        arma::vec  k_res_sub = t_res_mt.col(k_itrait);
        arma::vec k_res = k_res_sub.elem(k_sampleindices_vec);
        arma::vec  k_mu2_sub = t_mu2_mt.col(k_itrait);
        arma::vec k_mu2 = k_mu2_sub.elem(k_sampleindices_vec);

	//m_res_mt_vec.push_back(k_res);
        m_res_sum_vec(k_itrait) = arma::accu(k_res);
        m_mu2_sum_vec(k_itrait) = arma::accu(k_mu2);

	//arma::mat k_XV_submat = m_XV_mt.submat(k_ip, k_sampleindices_vec);
	//m_XV_mt_vec.push_back(k_XV_submat);
    }
   // m_XV_mt.clear();
}

// http://thecoatlessprofessor.com/programming/set_rs_seed_in_rcpp_sequential_case/
void SAIGEClass::set_seed(unsigned int seed){
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(seed);
}

void SAIGEClass::scoreTest(arma::vec & t_GVec,
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
		     arma::uvec & t_indexForNonZero){
    arma::vec Sm, var2m;
    double S, var2;
    //std::cout << "getadjGFast0 " << std::endl;
    //std::cout << "t_GVec.n_elem " << t_GVec.n_elem << std::endl;
    //std::cout << "t_indexForNonZero.n_elem " << t_indexForNonZero.n_elem << std::endl;
    getadjGFast(t_GVec, t_gtilde, t_indexForNonZero);
    //t_gtilde.print("t_gtilde");
    //t_indexForNonZero.print("t_indexForNonZero");
    //std::cout << "getadjGFast1 " << std::endl;
    //getadjG(t_GVec, t_gtilde);
// arma::mat XV_submat = m_XV_mt.submat(m_sampleindices_vec, m_ip);

    if(t_is_region && m_traitType == "binary"){
      t_gy = dot(t_gtilde, m_y);
     }
    S = dot(t_gtilde, m_res);
    //S = dot(t_gtilde, m_res_mt_vec.at(m_itrait));
    //std::cout << "S " << S << std::endl;
    S = S/m_tauvec_mt(0,m_itrait);;

    //std::cout << "S b " << S << std::endl;

    if(!m_flagSparseGRM_cur){
      t_P2Vec = t_gtilde % m_mu2 *(m_tauvec_mt(0,m_itrait)); 

      var2m = dot(t_P2Vec , t_gtilde);
       //std::cout << "!m_flagSparseGRM_cur" << std::endl;
    }else{
      //arma::sp_mat m_SigmaMat_sp = gen_sp_SigmaMat();
      arma::sp_mat m_SigmaMat_sp = gen_sp_SigmaMat_multiTrait();
      arma::vec m_diagSigma = arma::vec(m_SigmaMat_sp.diag());
      //t_P2Vec = arma::spsolve(m_SigmaMat_sp, t_gtilde);
      t_P2Vec = getPCG1ofSigmaAndGtilde_wo_precomp(m_SigmaMat_sp, m_diagSigma, t_gtilde, 100, 0.02); 
      var2m = dot(t_P2Vec , t_gtilde);
      if(m_isVarPsadj){
	var2m = var2m - t_gtilde.t() * m_Sigma_iXXSigma_iX * m_X.t() * t_P2Vec;	
      }
      //std::cout << "m_flagSparseGRM_cur" << std::endl;
    }	      
    var2 = var2m(0,0);
	//m_mu2.print("m_mu2");
      //std::cout << "m_tauvec_mt(0,m_itrait) " <<  m_tauvec_mt(0,m_itrait) << std::endl;
    //std::cout << "var2 " << var2 << std::endl;	
    double var1 = var2 * m_varRatioVal;
    //std::cout << "m_varRatioVal " << m_varRatioVal << std::endl;	
    //std::cout << "var1 " << var1 << std::endl;	

    double stat = S*S/var1;
    if (var1 <= std::numeric_limits<double>::min()){
          t_pval = 1;
    }else{
      if(!std::isnan(stat) && std::isfinite(stat)){
          boost::math::chi_squared chisq_dist(1);
          t_pval = boost::math::cdf(complement(chisq_dist, stat));
      }else{
          t_pval = 1;
	  stat = 0.0;
      }	      
    }
    char pValueBuf[100];
 
  //if(!std::isnan(stat)){
    if (t_pval != 0){
        sprintf(pValueBuf, "%.6E", t_pval);
	t_islogp = false;
    }else{	    
        double logp = R::pchisq(stat,1,false,true);
        double log10p = logp/(log(10));
        int exponent = floor(log10p);
        double fraction = pow(10.0, log10p - exponent);
        if (fraction >= 9.95) {
          fraction = 1;
           exponent++;
        }
        sprintf(pValueBuf, "%.1fE%d", fraction, exponent);
	t_pval = logp;
	t_islogp = true;
    }
    std::string buffAsStdStr = pValueBuf;
    t_pval_str = buffAsStdStr;
    t_Beta = S/var1;
    t_seBeta = fabs(t_Beta) / sqrt(fabs(stat));
    t_Tstat = S;
    t_var1 = var1;
    t_var2 = var2;
}


void SAIGEClass::scoreTestFast(arma::vec & t_GVec,
                     arma::uvec & t_indexForNonZero,
                     double& t_Beta,
                     double& t_seBeta,
                     std::string& t_pval_str,
		     double& t_pval,
                     bool& t_islogp,
                     double t_altFreq,
                     double &t_Tstat,
                     double &t_var1,
                     double &t_var2){

    arma::vec g1 = t_GVec.elem(t_indexForNonZero);
    arma::mat m_X_submat = m_X_mt.submat(m_sampleindices_vec, m_ip);
        arma::mat X1 = m_X_submat.rows(t_indexForNonZero);
    //arma::mat X1 = m_X.rows(t_indexForNonZero);
    arma::mat m_XVX_inv_XV_submat = m_XVX_inv_XV_mt.submat(m_sampleindices_vec, m_ip);
    arma::mat A1 = m_XVX_inv_XV_submat.rows(t_indexForNonZero);
    //arma::mat A1 = m_XVX_inv_XV.rows(t_indexForNonZero);
    arma::vec mu21;
    arma::vec res1 = m_res.elem(t_indexForNonZero);
    arma::vec Z = A1.t() * g1;
    //std::cout << "e " << std::endl;
    arma::vec B = X1 * Z;
    //std::cout << "e " << std::endl;
    arma::vec g1_tilde = g1 - B;
    double var1, var2, S, S1, S2, g1tildemu2;
    arma::vec S_a2;
    double Bmu2;
	//std::cout << "m_XVX_mt " << m_XVX_mt.n_cols << " " << m_XVX_mt.n_rows << std::endl;
	//m_ip.print("m_ip");
	//std::cout << "m_p " << m_p << std::endl;
    arma::mat XVX_submat = m_XVX_mt.cols(m_ip);

    arma::mat XVX_submat2 = XVX_submat.rows(0, (m_p-1));
    //arma::mat  ZtXVXZ = Z.t() * m_XVX * Z;
    //std::cout << "e " << std::endl;

	//std::cout << "Z " << Z.n_cols << " " << Z.n_rows << std::endl;
	//std::cout << "XVX_submat2 " << XVX_submat2.n_cols << " " << XVX_submat2.n_rows << std::endl;


    arma::mat  ZtXVXZ = Z.t() * XVX_submat2 * Z;
    //std::cout << "e " << std::endl;
    if(m_traitType == "binary" || m_traitType == "survival"){
      mu21  = m_mu2.elem(t_indexForNonZero);
      g1tildemu2 = dot(square(g1_tilde), mu21);
      Bmu2 = arma::dot(square(B),  mu21);
      var2 = ZtXVXZ(0,0) - Bmu2 + g1tildemu2;
    }else if(m_traitType == "quantitative"){
      Bmu2 = dot(g1, B);
      var2 = ZtXVXZ(0,0)*(m_tauvec_mt(0,m_itrait)) +  dot(g1,g1) - 2*Bmu2;
      //var2 = ZtXVXZ(0,0)*m_tauvec[0] +  dot(g1,g1) - 2*Bmu2;
    }

    var1 = var2 * m_varRatioVal;
    S1 = dot(res1, g1_tilde);
    arma::mat res1X1_temp = (res1.t()) * X1;
    arma::vec res1X1 = res1X1_temp.t();

    arma::vec S_a_vec = m_S_a_mt.col(m_itrait);
    //m_S_a = S_a_vec.subvec(0,m_p);

    //S_a2 = m_S_a - res1X1;
    //
    //
    S_a2 = S_a_vec.subvec(0,m_p-1) - res1X1;
    S2 = - arma::dot(S_a2,  Z);
    S = S1 + S2;
    //S = S/m_tauvec[0];
    S = S/(m_tauvec_mt(0,m_itrait));

    double stat = S*S/var1;
    if (var1 <= std::numeric_limits<double>::min()){
          t_pval = 1;
    }else{
      if(!std::isnan(stat) && std::isfinite(stat)){	    
          boost::math::chi_squared chisq_dist(1);
          t_pval = boost::math::cdf(complement(chisq_dist, stat));

      }else{
          t_pval = 1;
	  stat = 0.0;
      }
    }
    char pValueBuf[100];
    if (t_pval != 0){
        sprintf(pValueBuf, "%.6E", t_pval);
	t_islogp = false;
    }else {
	double logp = R::pchisq(stat,1,false,true);
        double log10p = logp/(log(10));
        int exponent = floor(log10p);
        double fraction = pow(10.0, log10p - exponent);
        if (fraction >= 9.95) {
          fraction = 1;
           exponent++;
         }
        sprintf(pValueBuf, "%.1fE%d", fraction, exponent);
	t_pval = logp;
	t_islogp = true;
    }
    std::string buffAsStdStr = pValueBuf;
    t_pval_str = buffAsStdStr;
    t_Beta = S/var1;
    t_seBeta = fabs(t_Beta) / sqrt(fabs(stat));	
    t_Tstat = S;
    t_var1 = var1;
    t_var2 = var2;
}


void SAIGEClass::scoreTestFast_noadjCov(arma::vec & t_GVec,
		     arma::uvec & t_indexForNonZero,
                     double& t_Beta,
                     double& t_seBeta,
                     std::string& t_pval_str,
                     double& t_pval,
                     bool& t_islogp,
                     double t_altFreq,
                     double &t_Tstat,
                     double &t_var1,
                     double &t_var2){

    //arma::vec Sm, var2m;
      arma::vec g1 = t_GVec.elem(t_indexForNonZero);
      arma::vec m_mu21 = m_mu2.elem(t_indexForNonZero);
      arma::vec m_res1 = m_res.elem(t_indexForNonZero);
    double S, var2;
    //getadjGFast(t_GVec, t_gtilde, t_indexForNonZero);
    //getadjG(t_GVec, t_gtilde);


    //if(t_is_region && m_traitType == "binary"){
    //  t_gy = dot(t_gtilde, m_y);
    // }
    double var2_a = dot(m_mu21,pow(g1,2));
    double var2_b = dot(m_mu21, 2*2*t_altFreq*g1);
    double var2_c = arma::accu(m_mu2)*pow(2*t_altFreq, 2);
    var2 = var2_a - var2_b + var2_c;

    //arma::vec t_GVec_center = t_GVec-arma::mean(t_GVec);
    //var2 = dot(m_mu2, pow(t_GVec_center,2));  
    //S = dot(t_GVec_center, m_res);
    S = dot(g1, m_res1)  - arma::accu(m_res)*(2*t_altFreq);
    //std::cout << "S " << S << std::endl;
    S = S/m_tauvec[0];

    //std::cout << "S b " << S << std::endl;

    //var2 = var2m(0,0);
    double var1 = var2 * m_varRatioVal;

    double stat = S*S/var1;
    if (var1 <= std::numeric_limits<double>::min()){
          t_pval = 1;
    }else{
      if(!std::isnan(stat) && std::isfinite(stat)){
          boost::math::chi_squared chisq_dist(1);
          t_pval = boost::math::cdf(complement(chisq_dist, stat));
      }else{
          t_pval = 1;
          stat = 0.0;
      }
    }
    char pValueBuf[100];

  //if(!std::isnan(stat)){
    if (t_pval != 0){
        sprintf(pValueBuf, "%.6E", t_pval);
        t_islogp = false;
    }else{
        double logp = R::pchisq(stat,1,false,true);
        double log10p = logp/(log(10));
        int exponent = floor(log10p);
        double fraction = pow(10.0, log10p - exponent);
        if (fraction >= 9.95) {
          fraction = 1;
           exponent++;
        }
        sprintf(pValueBuf, "%.1fE%d", fraction, exponent);
        t_pval = logp;
        t_islogp = true;
    }
    std::string buffAsStdStr = pValueBuf;
    t_pval_str = buffAsStdStr;
    t_Beta = S/var1;
    t_seBeta = fabs(t_Beta) / sqrt(fabs(stat));
    t_Tstat = S;
    t_var1 = var1;
    t_var2 = var2;
}





void SAIGEClass::getadjG(arma::vec & t_GVec, arma::vec & g){
    //g = m_XV * t_GVec;
   //    arma::mat subMat = inputMat.submat(arma::span(rowStart, rowEnd), colIndices);

   
   g = m_XV_mt.submat(m_ip, m_sampleindices_vec) * t_GVec;
   
   //t_GVec.t().print("t_GVec");
   //std::cout << "sum(t_GVec) " << arma::accu(t_GVec) << std::endl;
   //      m_XV.print("m_XV");
  // g = t_GVec - m_XXVX_inv * g;
  g = t_GVec - m_XXVX_inv_mt.submat(m_sampleindices_vec, m_ip) * g;
    //m_XXVX_inv.print("m_XXVX_inv");
    //g.t().print("g");
}


void SAIGEClass::getadjGFast(arma::vec & t_GVec, arma::vec & g, arma::uvec & iIndex)
{
/*
  	std::cout << "m_XV.n_cols " << m_XV_mt.n_cols << std::endl;
  	std::cout << "m_XV.n_rows " << m_XV_mt.n_rows << std::endl;
  	std::cout << "m_XXVX_inv_mt.n_cols " << m_XXVX_inv_mt.n_cols << std::endl;
  	std::cout << "m_XXVX_inv_mt.n_rows " << m_XXVX_inv_mt.n_rows << std::endl;
 std::cout << "m_XXVX_inv(0,0) " << m_XXVX_inv_mt[0,0] << std::endl; 
*/  
  
  // To increase computational efficiency when lots of GVec elements are 0
// std::cout << "m_p " << m_p << std::endl;
 arma::mat XV_submat = m_XV_mt.submat(m_ip, m_sampleindices_vec);
	/*m_ip.print("m_ip");

	std::cout << "m_XV_submat.n_cols " << XV_submat.n_cols << std::endl;
  	std::cout << "m_XV_submat.n_rows " << XV_submat.n_rows << std::endl;

std::cout << "XV_submat(0,0) " << XV_submat[0,0] << std::endl;
std::cout << "XV_submat(0,1) " << XV_submat[0,1] << std::endl;
std::cout << "XV_submat(1,0) " << XV_submat[1,0] << std::endl;
std::cout << "XV_submat(1,1) " << XV_submat[1,1] << std::endl;
*/
  arma::vec m_XVG(XV_submat.n_rows, arma::fill::zeros);
  arma::mat m_XXVX_inv_submat = m_XXVX_inv_mt.submat(m_sampleindices_vec, m_ip);
  //std::cout << "m_XXVX_inv_submat(0,0) " << m_XXVX_inv_submat[0,0] << std::endl;
  //std::cout << "m_XXVX_inv_submat(0,1) " << m_XXVX_inv_submat[0,1] << std::endl;
  //std::cout << "m_XXVX_inv_submat(1,0) " << m_XXVX_inv_submat[1,0] << std::endl;
  //std::cout << "m_XXVX_inv_submat(1,1) " << m_XXVX_inv_submat[1,1] << std::endl;	
	//std::cout << "XV_submat.n_cols " << XV_submat.n_cols << std::endl;
  	//std::cout << "XV_submat.n_rows " << XV_submat.n_rows << std::endl;
	
  for(int i = 0; i < iIndex.n_elem; i++){
	m_XVG += XV_submat.col(iIndex(i)) *  t_GVec(iIndex(i));    
  }
  //g = t_GVec - m_XXVX_inv_mt.submat(m_sampleindices_vec, m_ip) * m_XVG; 
// m_XVG.print("m_XVG");
 
 g = t_GVec - m_XXVX_inv_submat * m_XVG; 
}


void SAIGEClass::get_mu(arma::vec & t_mu){
    t_mu = m_mu;
}

void SAIGEClass::getindices(arma::uvec & t_case_indices,
      arma::uvec & t_ctrl_indices){
     t_case_indices = m_case_indices;
     t_ctrl_indices = m_ctrl_indices;
  }


void SAIGEClass::setupSparseMat(int r, arma::umat & locationMatinR, arma::vec & valueVecinR) {
    m_locationMat = locationMatinR;
    m_valueVec = valueVecinR;
    m_dimNum = r;
}



arma::sp_mat SAIGEClass::gen_sp_SigmaMat() {
    arma::sp_mat resultMat(m_locationMat, m_valueVec, m_dimNum, m_dimNum);
    return resultMat;
}

// revised
// need to add sparse Sigma version 
// This function only uses variance ratio and does not use sparse GRM
void SAIGEClass::getMarkerPval(arma::vec & t_GVec,
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
				bool t_isER,
                                bool t_isnoadjCov,
                                 bool t_isSparseGRM) 
{



  t_isFirth = false;
  //arma::vec adjGVec = getadjGFast(t_GVec);
  std::string t_pval_str;
  double t_var2, t_SPApval;
  //iIndex = arma::find(t_GVec != 0);
  //arma::vec t_gtilde;
  bool isScoreFast = true;
  if(m_flagSparseGRM_cur){
    isScoreFast = false;
  }

 
  double pval_noadj, pval, t_qval_Firth; //can be log or not raw
  bool ispvallog;

if(!m_flagSparseGRM_cur && t_isnoadjCov){
        is_gtilde = false;
        isScoreFast = true;
	//std::cout << "scoreTestFast_noadjCov_multiTrait "  << std::endl;
        scoreTestFast_noadjCov_multiTrait(t_GVec, iIndex, t_Beta, t_seBeta, t_pval_noSPA, pval_noadj, ispvallog, t_altFreq,t_Tstat, t_var1, t_var2);

}else if(m_flagSparseGRM_cur){
        is_gtilde = true;
        isScoreFast = false;
        //std::cout << "scoreTest "  << std::endl;
        scoreTest(t_GVec, t_Beta, t_seBeta, t_pval_noSPA, pval_noadj, ispvallog, t_altFreq, t_Tstat, t_var1, t_var2, t_gtilde, t_P2Vec, t_gy, is_region, iIndex);

}else{
        is_gtilde = false;
        isScoreFast = true;
        //std::cout << "scoreTestFast "  << std::endl;
        scoreTestFast(t_GVec, iIndex, t_Beta, t_seBeta, t_pval_noSPA, pval_noadj, ispvallog, t_altFreq, t_Tstat, t_var1, t_var2);
}
  //for test
  //arma::vec timeoutput3 = getTime();
 // std::cout << "pval_noadj " << pval_noadj << std::endl;
  double StdStat = std::abs(t_Tstat) / sqrt(t_var1);

  t_isSPAConverge = false;
  
 //arma::vec timeoutput3_a = getTime();
  double q, qinv, m1, NAmu, NAsigma, tol1, p_iIndexComVecSize;

  //arma::uvec iIndexComVec = arma::find(t_GVec == 0);
  //arma::uvec iIndexVec = arma::find(t_GVec != 0);
 
  
  unsigned int iIndexComVecSize = iIndexComVec.n_elem;
  unsigned int iIndexSize = iIndex.n_elem; 
  arma::vec gNB(iIndexSize, arma::fill::none);
  arma::vec gNA(iIndexComVecSize, arma::fill::none);
  arma::vec muNB(iIndexSize, arma::fill::none);
  arma::vec muNA(iIndexComVecSize, arma::fill::none);


  double gmuNB;



if((StdStat > m_SPA_Cutoff || std::isnan(StdStat)) && m_traitType != "quantitative" && t_isER){
	t_isER = true;
}else{	
	t_isER = false;
}


if(!t_isER){


  if(!std::isnan(StdStat) && (StdStat > m_SPA_Cutoff) && m_traitType != "quantitative"){

       if(!is_gtilde){
          //t_gtilde.resize(m_n);
          getadjGFast(t_GVec, t_gtilde, iIndex);
	  is_gtilde = true;
       }
	//int t_gtilden = t_gtilde.n_elem;
        p_iIndexComVecSize = double(iIndexComVecSize)/m_n;
	//std::cout << m_mu.n_elem << std::endl;
	//std::cout << t_gtilde.n_elem << std::endl;
   	m1 = dot(m_mu, t_gtilde);
	//std::cout << "SPA 0" << std::endl;
	if(p_iIndexComVecSize >= 0.5){
		unsigned int j1 = 0;
		unsigned int j2 = 0;
/*
		for(unsigned int j = 0; j < m_n ; j++){	
			//std::cout << "j " << j << std::endl;
			if(t_GVec(j) != 0){
			//std::cout << "j1 " << j1 << std::endl;
				gNB(j1) = t_gtilde(j);
				muNB(j1) = m_mu(j);
				j1 = j1 + 1;	
			}else{
			//std::cout << "j2 " << j2 << std::endl;
				gNA(j2) = t_gtilde(j);
				muNA(j2) = m_mu(j);
				j2 = j2 + 1;
          //process_mem_usage(mem1, mem2);
//   std::cout << "VM 4 a 1.3c: " << mem1/1000 << "; RSS 4 a 1.3: " << mem2/1000 << std::endl;
			}	
		}
		*/

/*
	std::cout << "gNB.n_elem " <<  gNB.n_elem << std::endl;	
	std::cout << "gNA.n_elem " <<  gNA.n_elem << std::endl;	
*/
	gNB = t_gtilde(iIndex);
	gNA = t_gtilde(iIndexComVec);
   	muNB = m_mu(iIndex);
   	muNA = m_mu(iIndexComVec);

/*	
	    std::cout << "gNA.n_elem 2 " << gNA.n_elem << std::endl;
        std::cout << "gNB.n_elem 2 " << gNB.n_elem << std::endl;
        std::cout << "muNA.n_elem 2 " << muNA.n_elem << std::endl;
        std::cout << "muNB.n_elem 2 " << muNB.n_elem << std::endl;
*/


  	gmuNB = dot(gNB,muNB);	 
	//std::cout << "SPA 2" << std::endl;
   	NAmu= m1-gmuNB;
	//std::cout << "SPA 3" << std::endl;

   }
	/*else{
		gNA.clear();
		gNB.clear();
		muNA.clear();
		gNB.clear();

	}*/

   	if(m_traitType == "binary"){
                q = t_Tstat/sqrt(t_var1/t_var2) + m1;

                if((q-m1) > 0){
                        qinv = -1 * std::abs(q-m1) + m1;
                }else if ((q-m1) == 0){
                        qinv =  m1;
                }else{
                        qinv = std::abs(q-m1) + m1;
                }
		if(p_iIndexComVecSize >= 0.5){
           		NAsigma = t_var2 - arma::sum(muNB % (1-muNB) % arma::pow(gNB,2));
		}
        }else if(m_traitType == "survival"){
                q = t_Tstat/sqrt(t_var1/t_var2);
                qinv = -q;
  		if(p_iIndexComVecSize >= 0.5){
           		NAsigma = t_var2 - arma::sum(muNB % arma::pow(gNB,2));
		}
/*		std::cout << "NAsigma is " << NAsigma << std::endl;
		std::cout << "t_var2 is " << t_var2 << std::endl;
		muNB.print("muNB");
		gNB.print("gNB");
  */
  }
    	//bool logp=false;
	double tol0 = std::numeric_limits<double>::epsilon();
	tol1 = std::pow(tol0, 0.25);
	if(p_iIndexComVecSize >= 0.5 && !m_flagSparseGRM_cur){
	//std::cout << "SPA 3" << std::endl;
		//std::cout << "SPA_fast" << std::endl;
        	SPA_fast(m_mu, t_gtilde, q, qinv, pval_noadj, ispvallog, gNA, gNB, muNA, muNB, NAmu, NAsigma, tol1, m_traitType, t_SPApval, t_isSPAConverge);
	
	}else{
		//std::cout << "SPA" << std::endl;
		SPA(m_mu, t_gtilde, q, qinv, pval_noadj, tol1, ispvallog, m_traitType, t_SPApval, t_isSPAConverge);	
	}


    boost::math::normal ns;
    double t_qval;



    if(t_isSPAConverge){
        try {
          t_qval = R::qnorm(t_SPApval/2, 0, 1, false, ispvallog);		
          t_qval = fabs(t_qval);
          t_seBeta = fabs(t_Beta)/t_qval;
        }catch (const std::overflow_error&) {
          t_qval = std::numeric_limits<double>::infinity();
	  t_isSPAConverge = false;
        }
    }


    if(!ispvallog && t_SPApval == 0){
          t_isSPAConverge = false;
    }


  }
   //t_pval_noSPA = pval_noadj; 
   char pValueBuf_SPA[100];
   if(m_traitType!="quantitative"){
        if(t_isSPAConverge){
		if(!ispvallog){
		  sprintf(pValueBuf_SPA, "%.6E", t_SPApval);
		}else{
        	  double t_SPApval_log10 = t_SPApval/(log(10));
        	  int exponent = floor(t_SPApval_log10);
        	  double fraction = pow(10.0, t_SPApval_log10 - exponent);
        	  if (fraction >= 9.95) {
          	    fraction = 1;
           	    exponent++;
         	  }
        	  sprintf(pValueBuf_SPA, "%.1fE%d", fraction, exponent);
		}
		std::string buffAsStdStr_SPA = pValueBuf_SPA;
	        t_pval = buffAsStdStr_SPA;
		pval = t_SPApval;	
        }else{
                t_pval = t_pval_noSPA;
		pval = pval_noadj;
		//pval_log10 = pval_noadj_log10;
        }

	if(!ispvallog){
		if(m_is_Firth_beta && pval <= m_pCutoffforFirth){
			t_isFirth = true;
			t_qval_Firth = R::qnorm(pval/2, 0, 1, false, false);
		}

	}else{
		if(m_is_Firth_beta && pval <= std::log(m_pCutoffforFirth)){
			t_isFirth = true;
			t_qval_Firth = R::qnorm(pval/2, 0, 1, false, true);
		}
	}


   }else{
        t_pval = t_pval_noSPA;
	pval = pval_noadj;
	//pval_log10 = pval_noadj_log10;
   }

}else{ //if(!t_isER){

    arma::mat Z_er(t_GVec.n_elem, 1);
    Z_er.col(0) = t_GVec;
    arma::vec res_er = m_res;
    arma::vec pi1_er = m_mu;
    arma::vec resout_er = m_resout;
    uint32_t n_case = arma::sum(m_y == 1); 

    /*
    std::cout << "t_GVec(0) " << t_GVec(0) << std::endl;
    std::cout << "res_er(0) " << res_er(0) << std::endl;
    std::cout << "pi1_er(0) " << pi1_er(0) << std::endl;
    resout_er.print("resout_er");
    std::cout << "iIndex.n_elem" << iIndex.n_elem << std::endl;
    std::cout << "iIndexComVec.n_elem" << iIndexComVec.n_elem << std::endl;
    std::cout << "m_n_case " << n_case << std::endl;
    iIndex.print("iIndex");
    */
    double pval_ER =  SKATExactBin_Work(Z_er, res_er, pi1_er, n_case, iIndex, iIndexComVec, resout_er, 2e+6, 1e+4, 1e-6, 1);
    char pValueBuf_ER[100];
    sprintf(pValueBuf_ER, "%.6E", pval_ER);
    std::string buffAsStdStr_ER = pValueBuf_ER;
    t_pval = pValueBuf_ER;
    pval = pval_ER;
    //std::cout << "t_pval " << t_pval << std::endl;
    boost::math::normal ns;
    double t_qval_ER;
    try{
      t_qval_ER = boost::math::quantile(ns, pval_ER/2);
      t_qval_ER = fabs(t_qval_ER);
      t_seBeta = fabs(t_Beta)/t_qval_ER;
    }catch (const std::overflow_error&) {
      t_qval_ER = std::numeric_limits<double>::infinity();
      t_seBeta = 0;
    }

    if(m_is_Firth_beta & pval <= m_pCutoffforFirth){
	t_isFirth = true;
	t_qval_Firth = t_qval_ER;
    }	    
}

   if(t_isFirth){
	if(!is_gtilde){
                getadjGFast(t_GVec, t_gtilde, iIndex);
                is_gtilde = true;
        }
	arma::mat x(t_GVec.n_elem, 2, arma::fill::ones);	
	x.col(1) = t_gtilde;
	arma::vec init(2, arma::fill::zeros);
	fast_logistf_fit_simple(x, m_y, m_offset, true, init, 50, 15, 15, 1e-5, 1e-5, 1e-5, t_Beta ,t_seBeta, t_isFirthConverge);
	//back calculates se based on beta from firth adjustion and the p-value that accounts for case-control imbalance
	t_seBeta = fabs(t_Beta)/fabs(t_qval_Firth);
   }
   
 //arma::vec timeoutput4 = getTime();
 //printTime(timeoutput3, timeoutput3_a, "Test Marker  ScoreTest");
//printTime(timeoutput3, timeoutput4, "Test Marker 3 to 4");
//printTime(timeoutput3_a, timeoutput4, "Test Marker SPA");

   //condition
   if(t_isCondition){
	if(!is_gtilde){
        	getadjGFast(t_GVec, t_gtilde, iIndex);
        	is_gtilde = true;
        }
        t_G1tilde_P_G2tilde = sqrt(m_varRatioVal) * t_gtilde.t() * m_P2Mat_cond;
        arma::vec t_Tstat_ctemp =  t_G1tilde_P_G2tilde * m_VarInvMat_cond * m_Tstat_cond;
	arma::mat tempgP2 = t_gtilde.t() * m_P2Mat_cond;

    	t_Tstat_c = t_Tstat - t_Tstat_ctemp(0);
    	arma::vec t_varT_ctemp = t_G1tilde_P_G2tilde * m_VarInvMat_cond * (t_G1tilde_P_G2tilde.t());
    	t_varT_c = t_var1 - t_varT_ctemp(0);

    double S_c = t_Tstat_c;

    double stat_c = S_c*S_c/t_varT_c;

    double pval_noSPA_c;

     if (t_varT_c <= std::numeric_limits<double>::min()){
        pval_noSPA_c = 1;
	stat_c = 0;
     }else{
       if(!std::isnan(stat_c) && std::isfinite(stat_c)){	     
        boost::math::chi_squared chisq_dist(1);
        pval_noSPA_c = boost::math::cdf(complement(chisq_dist, stat_c));
       }else{
        pval_noSPA_c = 1;
	stat_c = 0;
       }	       
     }

    char pValueBuf_c[100];

    if (pval_noSPA_c != 0){
        sprintf(pValueBuf_c, "%.6E", pval_noSPA_c);
        ispvallog = false;
    }else {
	double logp_c = R::pchisq(stat_c,1,false,true);	
	double log10p_c = logp_c/(log(10));
	int exponent_c = floor(log10p_c);
        double fraction_c = pow(10.0, log10p_c - exponent_c);
        if (fraction_c >= 9.95) {
          fraction_c = 1;
           exponent_c++;
         }
        sprintf(pValueBuf_c, "%.1fE%d", fraction_c, exponent_c);
	ispvallog = true;
	pval_noSPA_c = logp_c;
    }
    std::string buffAsStdStr_c = pValueBuf_c;
    t_pval_noSPA_c = buffAsStdStr_c;
    t_Beta_c = S_c/t_varT_c;
    t_seBeta_c = fabs(t_Beta_c) / sqrt(stat_c);
    t_Tstat_c = S_c;


    bool t_isSPAConverge_c;
    if(m_traitType != "quantitative" && stat_c > std::pow(m_SPA_Cutoff,2)){
	double q_c, qinv_c, pval_noadj_c, SPApval_c;    
	if(m_traitType == "binary"){
                q_c = t_Tstat_c/sqrt(t_varT_c/t_var2) + m1;
                if((q_c-m1) > 0){
                        qinv_c = -1 * std::abs(q_c-m1) + m1;
                }else if ((q_c-m1) == 0){
                        qinv_c =  m1;
                }else{
                        qinv_c = std::abs(q_c-m1) + m1;
                }
        }else if(m_traitType == "survival"){
                q_c = t_Tstat_c/sqrt(t_varT_c/t_var2);
                qinv = -q_c;
        }


        if(p_iIndexComVecSize >= 0.5 && !m_flagSparseGRM_cur){
		//std::cout << "SPA_fast " << std::endl;
                SPA_fast(m_mu, t_gtilde, q_c, qinv_c, pval_noSPA_c, ispvallog, gNA, gNB, muNA, muNB, NAmu, NAsigma, tol1, m_traitType, SPApval_c, t_isSPAConverge_c);
        }else{
		//std::cout << "SPA " << std::endl;
                SPA(m_mu, t_gtilde, q_c, qinv_c, pval_noSPA_c, tol1, ispvallog, m_traitType, SPApval_c, t_isSPAConverge_c);
        }

	
        boost::math::normal ns;
        //double pval_SPA = t_SPApval;
        double t_qval_c;
     if(t_isSPAConverge_c){	
        try {
          if(!ispvallog){
           t_qval_c = boost::math::quantile(ns, SPApval_c/2);
          }else{
	   t_qval_c = R::qnorm(SPApval_c/2, 0, 1, false, true);
          }
           t_qval_c = fabs(t_qval_c);
           t_seBeta_c = fabs(t_Beta_c)/t_qval_c;
        }catch (const std::overflow_error&) {
          t_qval_c = std::numeric_limits<double>::infinity();
          //t_seBeta_c = 0;
	  t_isSPAConverge_c = false;
        }
      }

	if(!ispvallog && SPApval_c == 0){
		t_isSPAConverge_c = false;	
	}


      char pValueBuf_SPA_c[100];
        if(t_isSPAConverge_c){
                if(!ispvallog){
                  sprintf(pValueBuf_SPA_c, "%.6E", SPApval_c/2);
                }else{
                  double SPApval_c_log10 = SPApval_c/(log(10));
                  int exponent = floor(SPApval_c_log10);
                  double fraction = pow(10.0, SPApval_c_log10 - exponent);
                  if (fraction >= 9.95) {
                    fraction = 1;
                    exponent++;
                  }
                  sprintf(pValueBuf_SPA_c, "%.1fE%d", fraction, exponent);
                }
                std::string buffAsStdStr_SPA_c = pValueBuf_SPA_c;
                t_pval_c = buffAsStdStr_SPA_c;
        }else{
                t_pval_c = t_pval_noSPA_c;
        }

    }else{
	t_isSPAConverge_c = false;    
    	t_pval_c = t_pval_noSPA_c;	    
    }	    
 }


    gNA.clear();
    gNB.clear();
    muNA.clear();
    gNB.clear();



    if(is_region && !is_gtilde){
	getadjGFast(t_GVec, t_gtilde, iIndex);
	is_gtilde = true; 
    }

    if(is_region && isScoreFast){

      t_gy = dot(t_gtilde, m_y);
      if(!m_flagSparseGRM_cur){
        t_P2Vec = t_gtilde % m_mu2 *m_tauvec[0];
      }else{
        //arma::sp_mat m_SigmaMat_sp = gen_sp_SigmaMat();
	arma::sp_mat m_SigmaMat_sp = gen_sp_SigmaMat_multiTrait();
        //t_P2Vec = arma::spsolve(m_SigmaMat_sp, t_gtilde);
        arma::vec m_diagSigma = arma::vec(m_SigmaMat_sp.diag());
        t_P2Vec = getPCG1ofSigmaAndGtilde_wo_precomp(m_SigmaMat_sp, m_diagSigma, t_gtilde, 100, 0.02);

      }
    }
}


bool SAIGEClass::assignVarianceRatio(double MAC, bool issparseforVR){
    bool hasVarRatio = false;
    arma::vec m_varRatio;
    if(issparseforVR){
	m_varRatio = m_varRatio_sparse;
    }else{
	m_varRatio = m_varRatio_null;
    }
    for(unsigned int i = 0; i < m_cateVarRatioMaxMACVecInclude.n_elem; i++)
    {
        if(MAC <= m_cateVarRatioMaxMACVecInclude(i) && MAC > m_cateVarRatioMinMACVecExclude(i)){    	    
		m_varRatioVal = m_varRatio(i);
		hasVarRatio = true;
	}	
    }

    if(!hasVarRatio){	
	if(MAC <= m_cateVarRatioMinMACVecExclude(0)){
		m_varRatioVal = m_varRatio(0);
		hasVarRatio = true;
	}	
    }

    if(!hasVarRatio){
        if(MAC > m_cateVarRatioMaxMACVecInclude.back()){
		//m_varRatioVal = m_varRatio(a-1);
		m_varRatioVal = m_varRatio.back();
                hasVarRatio = true;
        }
    }

    return(hasVarRatio);    
}


bool SAIGEClass::assignVarianceRatio(double MAC, bool issparseforVR, bool isnoXadj){
    bool hasVarRatio = false;
    arma::vec m_varRatio;
    if(issparseforVR){
        //m_varRatio = m_varRatio_sparse;
        m_varRatio = m_varRatio_sparse_mt.col(m_itrait);
    }else{
        if(!isnoXadj){
            //m_varRatio = m_varRatio_null;
            m_varRatio = m_varRatio_null_mt.col(m_itrait);
        }else{
        //m_varRatio_null_noXadj.print("m_varRatio_null_noXadj");
            m_varRatio = m_varRatio_null_noXadj_mt.col(m_itrait);
            //m_varRatio.print("m_varRatio");
        }
    }


    //m_cateVarRatioMinMACVecExclude.print("m_cateVarRatioMinMACVecExclude");
    //m_cateVarRatioMaxMACVecInclude.print("m_cateVarRatioMaxMACVecInclude");
    //m_varRatio.print("m_varRatio");
    //m_varRatio_null_noXadj_mt.print("m_varRatio_null_noXadj_mt");
    for(unsigned int i = 0; i < m_cateVarRatioMaxMACVecInclude.n_elem; i++)
    {
        if(MAC <= m_cateVarRatioMaxMACVecInclude(i) && MAC > m_cateVarRatioMinMACVecExclude(i)){
                m_varRatioVal = m_varRatio(i);
                hasVarRatio = true;
        }
    }

    if(!hasVarRatio){
        if(MAC < m_cateVarRatioMinMACVecExclude(0)){
                m_varRatioVal = m_varRatio(0);
                hasVarRatio = true;
        }
    }

    if(!hasVarRatio){
        if(MAC > m_cateVarRatioMaxMACVecInclude.back()){
                //m_varRatioVal = m_varRatio(a-1);
                m_varRatioVal = m_varRatio.back();
                hasVarRatio = true;
        }
    }

   //m_varRatioVal = m_varRatio(0);
   hasVarRatio = true;
   //std::cout << "hasVarRatio " << hasVarRatio << std::endl;
   return(hasVarRatio);
}

bool SAIGEClass::assignVarianceRatio(double MAC, bool issparseforVR, bool isnoXadj, bool issample){
    bool hasVarRatio = false;
    arma::vec m_varRatio;
    if(issparseforVR){
        //m_varRatio = m_varRatio_sparse;
        m_varRatio = m_varRatio_sparse_mt.col(m_itrait);
    }else{
        if(!isnoXadj){
            //m_varRatio = m_varRatio_null;
           if(!issample){
	   std::cout << "m_varRatio_null_mt.col " << m_varRatio_null_mt.n_cols << std::endl;
	   std::cout << "m_varRatio_null_mt.n_rows " << m_varRatio_null_mt.n_rows << std::endl;
            m_varRatio = m_varRatio_null_mt.col(m_itrait);
	    m_varRatio.print("m_varRatio");	


           }else{
            m_varRatio = m_varRatio_null_sample_mt.col(m_itrait);
           }
        }else{
        //m_varRatio_null_noXadj.print("m_varRatio_null_noXadj");
            m_varRatio = m_varRatio_null_noXadj_mt.col(m_itrait);
            //m_varRatio.print("m_varRatio");
        }
    }


    //m_cateVarRatioMinMACVecExclude.print("m_cateVarRatioMinMACVecExclude");
    //m_cateVarRatioMaxMACVecInclude.print("m_cateVarRatioMaxMACVecInclude");
    //m_varRatio.print("m_varRatio");
    //m_varRatio_null_noXadj_mt.print("m_varRatio_null_noXadj_mt");
    std::cout << "m_cateVarRatioMaxMACVecInclude.n_elem" << m_cateVarRatioMaxMACVecInclude.n_elem << std::endl;
    for(unsigned int i = 0; i < m_cateVarRatioMaxMACVecInclude.n_elem; i++)
    {
        if(MAC <= m_cateVarRatioMaxMACVecInclude(i) && MAC > m_cateVarRatioMinMACVecExclude(i)){
                m_varRatioVal = m_varRatio(i);
                hasVarRatio = true;
        }
    }

    if(!hasVarRatio){
        if(MAC < m_cateVarRatioMinMACVecExclude(0)){
                m_varRatioVal = m_varRatio(0);
                hasVarRatio = true;
        }
    }

    if(!hasVarRatio){
        if(MAC > m_cateVarRatioMaxMACVecInclude.back()){
                //m_varRatioVal = m_varRatio(a-1);
                m_varRatioVal = m_varRatio.back();
                hasVarRatio = true;
        }
    }

   //m_varRatioVal = m_varRatio(0);
   hasVarRatio = true;
   std::cout << "hasVarRatio " << hasVarRatio << std::endl;
   return(hasVarRatio);
}


void SAIGEClass::assignSingleVarianceRatio(bool issparseforVR){ 
    arma::vec m_varRatio;
    if(issparseforVR){
        m_varRatio = m_varRatio_sparse;
    }else{
        m_varRatio = m_varRatio_null;
    }	
    m_varRatioVal = m_varRatio(0);
}

void SAIGEClass::assignSingleVarianceRatio(bool issparseforVR, bool isnoXadj){
    arma::rowvec m_varRatio;
        //std::cout << "issparseforVR i0 " << issparseforVR << std::endl;
    if(issparseforVR){
        m_varRatio = m_varRatio_sparse_mt.row(0);
    }else{
        if(isnoXadj){
            m_varRatio = m_varRatio_null_noXadj_mt.row(0);
        }else{
            m_varRatio = m_varRatio_null_mt.row(0);
        }
    }
    //std::cout << "assignSingleVarianceRatio" << std::endl;
    //m_varRatio.print("m_varRatio");
    m_varRatioVal = m_varRatio(m_itrait);
}

void SAIGEClass::set_isnoadjCov_cur(bool t_isnoadjCov_cur){
        m_isnoadjCov_cur = t_isnoadjCov_cur;
	}

void SAIGEClass::assignSingleVarianceRatio(bool issparseforVR, bool isnoXadj, bool issample){
    arma::rowvec m_varRatio;
    std::cout << "issparseforVR " << issparseforVR << std::endl;
    if(issparseforVR){
        m_varRatio = m_varRatio_sparse_mt.row(0);
    }else{
        if(isnoXadj){
            m_varRatio = m_varRatio_null_noXadj_mt.row(0);
        }else{
          if(!issample){
            m_varRatio = m_varRatio_null_mt.row(0);
          }else{
            m_varRatio = m_varRatio_null_sample_mt.row(0);
          }
        }
    }
    std::cout << "assignSingleVarianceRatio" << std::endl;
    m_varRatio.print("m_varRatio");
    m_varRatioVal = m_varRatio(m_itrait);
}

/*
void SAIGEClass::assignSingleVarianceRatio(bool issparseforVR, bool isnoXadj){
    arma::rowvec m_varRatio;
    //std::cout << "issparseforVR i0 " << issparseforVR << std::endl;
    if(issparseforVR){
        m_varRatio = m_varRatio_sparse;
    }else{
        if(isnoXadj){
            m_varRatio = m_varRatio_null_noXadj;
        }else{
            m_varRatio = m_varRatio_null;
        }
    }
    //std::cout << "assignSingleVarianceRatio" << std::endl;
    //m_varRatio.print("m_varRatio");
    m_varRatioVal = m_varRatio(0);
}
*/


void SAIGEClass::assignSingleVarianceRatio_withinput(double t_varRatioVal){
        m_varRatioVal = t_varRatioVal;
}





void SAIGEClass::assignConditionFactors(
      arma::mat & t_P2Mat_cond,
      arma::mat & t_VarInvMat_cond,
      arma::mat & t_VarMat_cond,
      arma::vec & t_Tstat_cond,
      arma::mat & t_G2_Weight_cond,
      arma::vec & t_MAF_cond,
      arma::mat & t_qsum_cond,
      arma::mat & t_gsum_cond,
      std::vector<std::string> & t_p_cond
      ){
	m_P2Mat_cond = t_P2Mat_cond;
	m_VarInvMat_cond = t_VarInvMat_cond;
	m_VarMat_cond = t_VarMat_cond;
	m_Tstat_cond = t_Tstat_cond;
	m_MAF_cond = t_MAF_cond;
	m_qsum_cond = t_qsum_cond;
	m_gsum_cond = t_gsum_cond;
	m_G2_Weight_cond = t_G2_Weight_cond;
	m_p_cond = t_p_cond;
	m_numMarker_cond = t_Tstat_cond.n_elem;
}

void SAIGEClass::assignConditionFactors_scalefactor(
	arma::vec & t_scalefactor_G2_cond	
		){
	m_scalefactor_G2_cond = t_scalefactor_G2_cond;
	arma::mat scalefactor_G2_cond_Mat = arma::diagmat(arma::sqrt(m_scalefactor_G2_cond));
	arma::mat weightMat_G2_G2 = m_G2_Weight_cond * m_G2_Weight_cond.t(); 
	arma::mat VarMat_cond_scaled = scalefactor_G2_cond_Mat * m_VarMat_cond * scalefactor_G2_cond_Mat;
	arma::mat VarMat_cond_scaled_weighted = VarMat_cond_scaled % weightMat_G2_G2;
	m_VarInvMat_cond_scaled_weighted = arma::pinv(VarMat_cond_scaled_weighted);
	//m_VarInvMat_cond_region_binary = (1/scalefactor_G2_cond_Mat) * m_VarInvMat_cond	* (1/scalefactor_G2_cond_Mat);
	
}

void SAIGEClass::extract_XV_XXVX_inv(arma::mat & t_XV, arma::mat & t_XXVX_inv){
	t_XV = m_XV;
	t_XXVX_inv = m_XXVX_inv;	
}



void SAIGEClass::fast_logistf_fit_simple(arma::mat & x,
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
	bool & isfirthconverge){
  isfirthconverge = false;	
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
  while(iter <= maxit){
        beta_old = beta;
        arma::vec wpi = pi % (1 - pi);
        arma::vec W2 = arma::sqrt(wpi);
        //arma::vec wpi_sqrt = arma::sqrt(wpi);
        //arma::vec W2 = weight % wpi_sqrt;
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
                ypih = (y - pi) + (h % (0.5 - pi));
        }else{
                ypih = (y - pi);
        }
        //ypih.print();
        arma::vec xcol(n, arma::fill::zeros);
        U_star = x.t() * ypih;

        arma::mat XX_XW2(n, k, arma::fill::zeros);
        for(int j = 0; j < k; j++){
                xcol = x.col(j);
                XX_XW2.col(j) = xcol % W2;
        }
        arma::mat XX_Fisher = XX_XW2.t() * (XX_XW2);
        bool isinv = arma::inv_sympd (XX_covs, XX_Fisher);
        
	if(!isinv){
                break;
        }
        //}
        arma::vec delta = XX_covs * U_star;
        //delta.replace(arma::datum::nan, 0);

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
                beta_G = arma::datum::nan;
                sebeta_G = arma::datum::nan;
        }else{
                beta_G = beta(1);
                sebeta_G = sqrt(XX_covs(1,1));
        }

	//std::cout << "beta_G " << beta_G << std::endl;
	//std::cout << "sebeta_G " << sebeta_G << std::endl;
        //return beta;
}

void SAIGEClass::set_flagSparseGRM_cur(bool t_flagSparseGRM_cur){
	m_flagSparseGRM_cur = t_flagSparseGRM_cur;
}

void SAIGEClass::assign_for_itrait(unsigned int t_itrait){
        m_itrait = t_itrait;
        m_traitType = m_traitType_vec.at(m_itrait);
	arma::uvec sampleindices_sub_vec = m_sampleindices_mt.col(t_itrait);
	m_sampleindices_vec = sampleindices_sub_vec.subvec(0, (m_sampleIndexLenVec[t_itrait]-1));

	if(t_itrait == 0){
		m_startip = 0;
	}else{
		m_startip = arma::sum(m_colXvec.subvec(0, t_itrait - 1));
	}
        m_p = m_colXvec(t_itrait);
	//std::cout << "assign_for_itrait 1 " << std::endl;	
	m_endip =  m_startip + m_colXvec[t_itrait]-1;
        
	m_ip.set_size(m_endip - m_startip + 1);
        unsigned int diff = m_endip - m_startip + 1;
	    
	//std::cout << "assign_for_itrait 2 " << std::endl;	
        // Populate m_ip with values from m_startip to m_endip
        for (unsigned int i = 0; i < m_ip.size(); ++i) {
            m_ip(i) = m_startip + i;
        }	
	
            //std::cout << "assign_for_itrait5 " << std::endl;	
            //std::cout << "m_y_mt.n_cols " << m_y_mt.n_cols << " " << m_itrait <<std::endl;	
	
       arma::vec  m_y_sub = m_y_mt.col(m_itrait);
          //  std::cout << "m_y_mt.n_cols " << m_y_mt.n_cols << " " << m_itrait <<std::endl;
	//m_sampleindices_vec.print("m_sampleindices_vec");
	//std::cout << "m_sampleindices_vec.n_elem " << m_sampleindices_vec.n_elem << std::endl;
       m_y = m_y_sub.elem(m_sampleindices_vec);
       
            //       std::cout << "assign_for_itrait5a " << std::endl;
       
       
       arma::vec  m_res_sub = m_res_mt.col(m_itrait);
          //         std::cout << "assign_for_itrait5a1 " << std::endl;
       m_res = m_res_sub.elem(m_sampleindices_vec);
       //m_res = m_res_mt_vec.at(m_itrait);
	//std::cout << "assign_for_itrait5b " << std::endl;

       arma::vec  m_mu2_sub = m_mu2_mt.col(m_itrait);
       m_mu2 = m_mu2_sub.elem(m_sampleindices_vec);
       arma::vec  m_mu_sub = m_mu_mt.col(m_itrait);
       m_mu = m_mu_sub.elem(m_sampleindices_vec);
       //std::cout << "assign_for_itrait5c " << std::endl;	

	m_startic = m_itrait*m_numMarker_cond;
        m_endic = m_startic + m_numMarker_cond - 1;
	//m_resout_mt.print("m_resout_mt");
	//if(m_traitType == "binary"){
	//	m_resout = m_resout_mt.col(m_itrait);
	//}

            //std::cout << "assign_for_itrait6 " << std::endl;	
        //m_startip = m_itrait*m_p;
        //m_endip = m_startip + m_p - 1;
        //std::cout << "assign_for_itrait m_is_gxe " << m_is_gxe << std::endl;
        /*if(m_is_gxe){
        //std::cout << "m_XV_gxe_mt.n_rows " << m_XV_gxe_mt.n_rows << std::endl;
        //std::cout << "m_itrait " << m_itrait << std::endl;
        m_p_gxe = (m_XV_gxe_mt.n_rows) / (m_traitType_vec.size());
        //std::cout << "m_p_gxe " << m_p_gxe << std::endl;
        m_startip_gxe = m_itrait*m_p_gxe;
        m_endip_gxe = m_startip_gxe + m_p_gxe - 1;
        m_n_gxe = (m_XV_gxe_mt.n_cols) / (m_traitType_vec.size());
        m_startin_gxe = m_itrait*m_n_gxe;
        m_endin_gxe = m_startin_gxe + m_n_gxe - 1;
        }*/
}
/*
void SAIGEClass::assignConditionFactors(
      arma::mat & t_P2Mat_cond,
      arma::mat & t_VarInvMat_cond,
      arma::mat & t_VarMat_cond,
      arma::vec & t_Tstat_cond,
      arma::mat & t_G2_Weight_cond,
      arma::vec & t_MAF_cond,
      arma::mat & t_qsum_cond,
      arma::mat & t_gsum_cond,
      std::vector<std::string> & t_p_cond
      ){
        m_P2Mat_cond = t_P2Mat_cond;
        m_VarInvMat_cond = t_VarInvMat_cond;
        m_VarMat_cond = t_VarMat_cond;
        m_Tstat_cond = t_Tstat_cond;
        m_MAF_cond = t_MAF_cond;
        m_qsum_cond_Mat = t_qsum_cond;
        m_gsum_cond_Mat = t_gsum_cond;
        m_G2_Weight_cond = t_G2_Weight_cond;

        m_G2_Weight_cond.print("m_G2_Weight_cond");

        m_p_cond = t_p_cond;
        m_numMarker_cond = t_Tstat_cond.n_elem;
}
*/
void SAIGEClass::assignConditionFactors_scalefactor_multiTrait(
        arma::mat & t_scalefactor_G2_cond,
        unsigned int oml
                ){
        unsigned int startoml = oml * (m_VarMat_cond.n_rows);
        unsigned int endoml = (oml+1) * (m_VarMat_cond.n_rows) - 1;
        m_scalefactor_G2_cond_Mat = t_scalefactor_G2_cond;

        m_scalefactor_G2_cond = t_scalefactor_G2_cond.col(oml);
        arma::mat scalefactor_G2_cond_Mat = arma::diagmat(arma::sqrt(m_scalefactor_G2_cond));
        arma::mat weightMat_G2_G2 = m_G2_Weight_cond * m_G2_Weight_cond.t();
        arma::mat VarMat_cond_scaled = scalefactor_G2_cond_Mat * (m_VarMat_cond.cols(startoml, endoml)) * scalefactor_G2_cond_Mat;
        arma::mat VarMat_cond_scaled_weighted = VarMat_cond_scaled % weightMat_G2_G2;
        m_VarInvMat_cond_scaled_weighted.cols(startoml, endoml) = VarMat_cond_scaled_weighted.i();
        //m_VarInvMat_cond_region_binary = (1/scalefactor_G2_cond_Mat) * m_VarInvMat_cond       * (1/scalefactor_G2_cond_Mat);
}


//to update

void SAIGEClass::assign_for_itrait_binaryindices(unsigned int t_itrait){
	m_y = m_y_mt.col(t_itrait);
	arma::vec m_y_sub = m_y.elem(m_sampleindices_vec);
        m_case_indices = arma::find(m_y_sub == 1);
        m_ctrl_indices = arma::find(m_y_sub == 0);
        m_n_case = m_case_indices.n_elem;
        m_n_ctrl = m_ctrl_indices.n_elem;
}


void SAIGEClass::scoreTestFast_noadjCov_multiTrait(arma::vec & t_GVec,
                     arma::uvec & t_indexForNonZero,
                     double& t_Beta,
                     double& t_seBeta,
                     std::string& t_pval_str,
                     double& t_pval,
                     bool& t_islogp,
                     double t_altFreq,
                     double &t_Tstat,
                     double &t_var1,
                     double &t_var2){

    //arma::vec t_GVec_center = t_GVec-arma::mean(t_GVec);
    //double var2 = dot(m_mu2, pow(t_GVec_center,2));
    //arma::vec Sm, var2m;
      arma::vec g1 = t_GVec.elem(t_indexForNonZero);
      arma::vec m_mu21 = m_mu2.elem(t_indexForNonZero);
      //arma::vec g_res = m_res_mt_vec.at(m_itrait);
      arma::vec m_res1 = m_res.elem(t_indexForNonZero);
    double S, var2;
    //getadjGFast(t_GVec, t_gtilde, t_indexForNonZero);
    //getadjG(t_GVec, t_gtilde);


    //if(t_is_region && m_traitType == "binary"){
    //  t_gy = dot(t_gtilde, m_y);
    // }
    double var2_a = dot(m_mu21,pow(g1,2));
    double var2_b = dot(m_mu21, 2*2*t_altFreq*g1);
    double var2_c = m_mu2_sum_vec[m_itrait]*pow(2*t_altFreq, 2);
    var2 = var2_a - var2_b + var2_c;

    //arma::vec t_GVec_center = t_GVec-arma::mean(t_GVec);
    //var2 = dot(m_mu2, pow(t_GVec_center,2));
    //S = dot(t_GVec_center, m_res);
    S = dot(g1, m_res1)  - m_res_sum_vec[m_itrait]*(2*t_altFreq);
    //std::cout << "S " << S << std::endl;


    double var1 = var2 * m_varRatioVal;
    //std::cout << "var2 " << var2 << std::endl;
    //std::cout << "var1 " << var1 << std::endl;


    //double S = dot(m_res_mt_vec.at(m_itrait), t_GVec_center);
    S = S/m_tauvec_mt(0,m_itrait);
    double stat = S*S/var1;
    //double t_pval;
    //std::cout << "S " << S << std::endl;
    //if (var1 <= std::pow(std::numeric_limits<double>::min(), 2)){
    //
    //
    if (var1 <= std::numeric_limits<double>::min()){
          t_pval = 1;
    }else{
      if(!std::isnan(stat) && std::isfinite(stat)){
          boost::math::chi_squared chisq_dist(1);
          t_pval = boost::math::cdf(complement(chisq_dist, stat));

      }else{
          t_pval = 1;
          stat = 0.0;
      }
    }
    char pValueBuf[100];
    if (t_pval != 0){
        sprintf(pValueBuf, "%.6E", t_pval);
        t_islogp = false;
    }else {
        double logp = R::pchisq(stat,1,false,true);
        double log10p = logp/(log(10));
        int exponent = floor(log10p);
        double fraction = pow(10.0, log10p - exponent);
        if (fraction >= 9.95) {
          fraction = 1;
           exponent++;
         }
        sprintf(pValueBuf, "%.1fE%d", fraction, exponent);
        t_pval = logp;
        t_islogp = true;
    }

    std::string buffAsStdStr = pValueBuf;
    t_pval_str = buffAsStdStr;
    t_Beta = S/var1;
    t_seBeta = fabs(t_Beta) / sqrt(fabs(stat));
    t_Tstat = S;
    t_var1 = var1;
    t_var2 = var2;
    //std::cout << "t_pval_str scoreTestFast " << t_pval_str << std::endl;
    //std::cout << "end of scoreTestFast" << std::endl;
}

void SAIGEClass::assign_for_itrait_sampleIndices(unsigned int t_itrait){
        //m_itrait = t_itrait;
        //m_traitType = m_traitType_vec.at(m_itrait);
        arma::uvec sampleindices_sub_vec = m_sampleindices_mt.col(t_itrait);
        m_sampleindices_vec = sampleindices_sub_vec.subvec(0, (m_sampleIndexLenVec[t_itrait]-1));
}

arma::vec SAIGEClass::getPCG1ofSigmaAndGtilde_wo_precomp(arma::sp_mat & m_spSigmaMat, arma::vec & m_diagSigma, arma::vec& bVec, int maxiterPCG, double tolPCG) {
    int Nnomissing = m_spSigmaMat.n_rows;
    arma::vec xVec(Nnomissing, arma::fill::zeros); // Initialize xVec to zeros
    arma::vec rVec = bVec; // Residual vector
    arma::vec zVec(Nnomissing);
    //arma::vec minvVec = 1.0 / spsigma.diag(); // Preconditioner (reciprocal of diagonal elements)
    arma::vec minvVec = 1.0 / m_diagSigma; // Convert diagonal view to dense vector
    //m_diagSigma.print("m_diagSigma");
    zVec = minvVec % rVec; // Apply preconditioner
    double sumr2 = arma::dot(rVec, rVec); // Initial residual norm
    arma::vec pVec = zVec; // Search direction

    int iter = 0;
    while (sumr2 > tolPCG && iter < maxiterPCG) {
        iter++;
        arma::vec ApVec = m_spSigmaMat * pVec; // Sparse matrix-vector multiplication
        //ApVec.print("ApVec");
        double alpha = arma::dot(rVec, zVec) / arma::dot(pVec, ApVec); // Step size
        //std::cout << "alpha" << alpha << std::endl;
        //pVec.print("pVec");
        //xVec += alpha * pVec; // Update solution
        xVec = xVec + alpha * pVec;

        arma::vec r1Vec = rVec - alpha * ApVec; // Update residual

        arma::vec z1Vec = minvVec % r1Vec; // Apply preconditioner to new residual
        double beta = arma::dot(z1Vec, r1Vec) / arma::dot(zVec, rVec); // Update beta

        pVec = z1Vec + beta * pVec; // Update search direction
        zVec = z1Vec; // Update preconditioned residual
        rVec = r1Vec; // Update residual
        //xVec.print("xVec");
        sumr2 = arma::dot(rVec, rVec); // Update residual norm
    }

    if (iter >= maxiterPCG) {
        std::cout << "PCG in getPCG1ofSigmaAndGtilde did not converge. You may increase maxiter number." << std::endl;
    }

    //std::cout << "Iterations from getPCG1ofSigmaAndGtilde: " << iter << std::endl;
    return xVec;
}


void  SAIGEClass::setsparseSigmaMtxMultiTraits(arma::umat & sparseSigmaLocationMtx, arma::vec & sparseSigmaValueVec, arma::umat & sparseSigmaIndiceMtx, arma::ivec & dimNumVec){
	m_sparseSigmaLocationMtx = sparseSigmaLocationMtx;
	m_sparseSigmaValueVec = sparseSigmaValueVec;
	m_sparseSigmaIndiceMtx = sparseSigmaIndiceMtx;
	m_dimNumVec = dimNumVec;
}

arma::sp_mat SAIGEClass::gen_sp_SigmaMat_multiTrait() {
    unsigned int starts = m_sparseSigmaIndiceMtx(m_itrait, 0);
    unsigned int ends = m_sparseSigmaIndiceMtx(m_itrait, 1);
    m_locationMat = m_sparseSigmaLocationMtx.rows(starts, ends);
    m_valueVec = m_sparseSigmaValueVec.subvec(starts, ends-starts+1);
    m_dimNum = m_dimNumVec[m_itrait];
    arma::sp_mat resultMat(m_locationMat, m_valueVec, m_dimNum, m_dimNum);
    return resultMat;
}

}
