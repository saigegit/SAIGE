#ifndef SKAT_HPP
#define SKAT_HPP

#include <RcppArmadillo.h>
#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/non_central_chi_squared.hpp>


#include "qfc_rcpp.hpp"


// [[Rcpp::depends(RcppArmadillo)]]
// // [[Rcpp::plugins(cpp11)]]

arma::vec Get_Lambda(arma::mat& K, bool isFast = false, int maxK = 100);
Rcpp::List SKAT_META_Optimal_Param(arma::mat& Phi, arma::vec& r_all);
Rcpp::List Get_Liu_Params_Mod(arma::vec& c1);
Rcpp::List Get_Liu_Params_Mod_Lambda(arma::vec& lambda, arma::ivec& df1);
arma::vec Get_Liu_PVal_MOD_Lambda(arma::vec& Q_all, arma::vec& lambda, arma::ivec& df1, bool log_p = false);
std::string Get_Liu_PVal_MOD_Lambda_Zero(double Q, double muQ, double muX, double sigmaQ, double sigmaX, double l, double d);
arma::vec SKAT_META_Optimal_Get_Q(arma::vec& Score, arma::vec& r_all);
arma::mat SKAT_META_Optimal_Get_Q_Res(arma::mat& Score_res, arma::vec& r_all);
double daviesPValue(arma::vec& eigenvalues, double q, double tol = 1e-25);
Rcpp::List Get_Davies_PVal(arma::mat & Q, arma::mat & W, arma::mat & Q_resampling, bool isFast = false);
Rcpp::List SKAT_davies(double q, arma::vec& lambda, arma::ivec& h, arma::vec& delta, double sigma = 0.0, int lim = 10000, double acc = 0.0001);
Rcpp::List Get_PValue_Lambda(arma::vec &lambda, arma::vec &Q, arma::ivec &df1);
Rcpp::List SKAT_Optimal_Each_Q(Rcpp::List &param_m, arma::mat &Q_all, arma::vec &r_all, Rcpp::List & lambda_all, std::string & method);
arma::vec SKAT_Optimal_Integrate_Func_Davies(arma::vec & x,  arma::mat &pmin_q, Rcpp::List &param_m,  arma::vec &r_all);
arma::vec integrate_SKAT_Optimal_Davies(arma::vec &pmin_q, Rcpp::List &param_m, arma::vec &r_all, double lower, double upper, int subdivisions, double abs_tol);
double SKAT_Optimal_PValue_Davies(arma::vec &pmin_q, Rcpp::List &param_m,
                                     arma::vec &r_all, double pmin);
arma::vec SKAT_Optimal_Integrate_Func_Liu(arma::vec & x,  arma::mat &pmin_q,  Rcpp::List &param_m,  arma::vec &r_all);
arma::vec integrate_SKAT_Optimal_Liu(arma::mat &pmin_q, Rcpp::List &param_m, arma::vec &r_all,double lower, double upper, int subdivisions, double abs_tol);
double SKAT_Optimal_PValue_Liu( arma::mat &pmin_q,  Rcpp::List &param_m,
                                   arma::vec &r_all, double pmin);
Rcpp::List SKAT_META_Optimal_Get_Pvalue( arma::mat &Q_all,  arma::mat &Phi,  arma::vec &r_all, std::string &method, bool isFast);
Rcpp::List SKAT_META_Optimal( arma::vec &Score,  arma::mat &Phi,  arma::vec &r_all,
                              std::string method,  arma::mat &Score_Resampling, bool isFast);
Rcpp::List Met_SKAT_Get_Pvalue( arma::vec &Score,  arma::mat &Phi,  arma::vec &r_corr, std::string &method,  arma::mat &Score_Resampling, bool isFast = false);
Rcpp::List Get_Liu_Params( arma::vec& c1);
Rcpp::List Get_Liu_PVal( arma::vec& Q,  arma::mat& W,  arma::mat& Q_resampling);



#endif
