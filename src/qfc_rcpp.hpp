#ifndef QFC_HPP
#define QFC_HPP

#include <RcppArmadillo.h>
#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>

// [[Rcpp::depends(RcppArmadillo)]]


namespace QUADFORM {

class QuadFormClass {
private:
    double sigsq, lmax, lmin, mean, c;
    double intl, ersm;
    int count, r, lim;
    static bool ndtsrt, fail;
    std::vector<int> n, th;
    std::vector<double> lb, nc;

public:
    QuadFormClass(int r, int lim);
    double exp1(double x);
    void counter();
    double log1(double x, bool first);
    void order();
    double errbd(double u, double& cx);
    double ctff(double accx, double& upn);
    double truncation(double u, double tausq);
    void findu(double& utx, double accx);
    void integrate(int nterm, double interv, double tausq, bool mainx);
    double cfe(double x);
    
   void qfc_1(
        arma::vec& lb1, arma::vec& nc1, arma::ivec& n1, int r,
           double  sigma, double c1, int lim,
           double acc, arma::vec& trace, int ifault,
           double & res
); 
    
};

}
#endif
