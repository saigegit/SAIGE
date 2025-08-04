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


//#include "qfc_rcpp.hpp"


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
double integrate_SKAT_Optimal_Davies(arma::vec &pmin_q, Rcpp::List &param_m, arma::vec &r_all, double lower, double upper, int subdivisions, double abs_tol);
double SKAT_Optimal_PValue_Davies(arma::vec &pmin_q, Rcpp::List &param_m,
                                     arma::vec &r_all, double pmin);
arma::vec SKAT_Optimal_Integrate_Func_Liu(arma::vec & x,  arma::mat &pmin_q,  Rcpp::List &param_m,  arma::vec &r_all);
double SKAT_Optimal_PValue_Liu(arma::vec &pmin_q, Rcpp::List &param_m,
                                     arma::vec &r_all, double pmin);

Rcpp::List SKAT_META_Optimal_Get_Pvalue( arma::mat &Q_all,  arma::mat &Phi,  arma::vec &r_all, std::string &method, bool isFast);
Rcpp::List SKAT_META_Optimal( arma::vec &Score,  arma::mat &Phi,  arma::vec &r_all,
                              std::string method,  arma::mat &Score_Resampling, bool isFast);
Rcpp::List Met_SKAT_Get_Pvalue( arma::vec &Score,  arma::mat &Phi,  arma::vec &r_corr, std::string &method, bool isFast = false);
Rcpp::List Get_Liu_Params( arma::vec& c1);
Rcpp::List Get_Liu_PVal( arma::vec& Q,  arma::mat& W,  arma::mat& Q_resampling);
Rcpp::List call_qfc(const arma::vec& lambdas, const arma::vec& noncentral,
                    const arma::ivec& df, int r, double sigma, double q,
                    int lim, double acc);

arma::mat forceSymmetric(const arma::mat& K);

void get_SKAT_pvalue_cpp(arma::vec& Score,
                               arma::mat& Phi,
                               arma::vec& r_corr,
                               double& Pvalue_SKATO,
                               double& Pvalue_Burden,
                               double& Pvalue_SKAT,
                               double& BETA_Burden,
                               double& SE_Burden,
                               int& error_code
                               );
double qchisq_log(double pval, double df);

void SPA_ER_kernel_related_Phiadj_fast_new_cpp(arma::vec& p_new,
                                                 arma::vec& Score,
                                                 arma::mat& Phi,
                                                 double p_value_burden,
                                                 std::string regionTestType,
                                                 arma::vec& scaleFactor
                                                 );
void get_newPhi_scaleFactor_cpp(double q_sum,
                                  arma::vec& mu_a,
                                  arma::vec& g_sum,
                                  arma::vec& p_new,
                                  arma::vec& Score,
                                  arma::mat& Phi,
                                  std::string regionTestType,
                                  arma::vec& scaleFactor);


double get_jointScore_pvalue(arma::vec& Score, arma::mat& Phi);

double integrate_SKAT_Optimal_Liu(arma::vec &pmin_q, Rcpp::List &param_m, arma::vec &r_all, double lower, double upper, int subdivisions, double abs_tol);

double integrate_SKAT_Optimal_Liu_v2(arma::vec &pmin_q, Rcpp::List &param_m, arma::vec &r_all, double lower, double upper, int subdivisions, double abs_tol);

double integrate_SKAT_Optimal_Liu_v3(arma::vec &pmin_q, Rcpp::List &param_m, arma::vec &r_all, double lower, double upper, int subdivisions, double abs_tol);

double integrate_SKAT_Optimal_Davies_v2(arma::vec &pmin_q, Rcpp::List &param_m, arma::vec &r_all, double lower, double upper, int subdivisions, double abs_tol);

/*
Rcpp::List integrate_SKATliu(arma::vec &pmin_q, Rcpp::List &param_m, arma::vec &r_all, double lower, double upper, int subdivisions, double abs_tol);

class Func
{
public:
    virtual double operator()(const double& x) const = 0;
    virtual void eval(double* x, const int n) const
    {
        for(int i = 0; i < n; i++)
            x[i] = this->operator()(x[i]);
    }
    
    virtual ~Func() {}
};



class SKATLiu : public Numer::Func
{
private:
    arma::mat pmin_q;
    Rcpp::List param_m ;
    arma::vec r_all;
public:
    // Constructor to initialize shape parameters
    SKATLiu(arma::mat &pmin_q_,  Rcpp::List &param_m_,  arma::vec &r_all_) : pmin_q(pmin_q_), param_m(param_m_), r_all(r_all_){}

    // Overload the function call operator to compute the PDF value
    double operator()(const double& x) const override
    {
        // Wrap the scalar x into an Armadillo vector
        arma::vec x_vec = {x};

        // Use calculate_beta_pdf to compute the value
	const arma::mat& pmin_q_c = pmin_q;
	const arma::vec& r_all_c = r_all;
	const Rcpp::List& param_m_c = param_m;
	arma::vec integratevec  = SKAT_Optimal_Integrate_Func_Liu(x_vec, pmin_q_c, param_m_c, r_all_c);


        // Return the first (and only) result
        return integratevec(0);
    }
};
*/



#endif
