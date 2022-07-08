#ifndef ER_HPP
#define ER_HPP



// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>


void GetProb_new(int k, int ngroup, int ncase, int* group, double* weight, double* prob);
void SKATExactBin_ComputeProb_Group(arma::uvec & idx, arma::uvec & idxCompVec, arma::vec & pi1, uint32_t n, uint32_t ncase, int type_group, std::vector<double> & prob);

int fact(int n);
int n_choose_r(int n, int r);

void Get_Total_K(int k, std::vector<int> & n_total_k);

void SKATExactBin_ComputProb_New(arma::uvec & idx, arma::uvec & idxCompVec, arma::vec & pi1, uint32_t n, uint32_t ncase, int NResampling, int ExactMax, int test_type, int type_group, std::vector<double> & prob, std::vector<int> & IsExactVec, std::vector<int> & n_total_k, int & n_total, bool & Is_ExactP);


void  Get_Res_Arrays(arma::mat & res_out, arma::uvec & idx, std::vector<int> & resarray, int & nres, std::vector<int> & nres_k);

std::vector< std::vector<double> > mat_to_std_vec(arma::mat &A);


double SKATExactBin_Work(arma::mat & Z, arma::vec & res, arma::vec & pi1, uint32_t ncase, arma::uvec & idx, arma::uvec & idxCompVec, arma::mat & res_out, int NResampling, int ExactMax, double epsilon, int test_type);

#endif
