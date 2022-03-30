// [[Rcpp::depends(BH)]]

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "getMem.hpp"
#include <chrono>         // std::chrono::seconds
#include <boost/iostreams/filter/zstd.hpp>
#include <boost/date_time.hpp>
#include <boost/math/distributions/cauchy.hpp>

#define _USE_MATH_DEFINES
#include <math.h> 


// [[Rcpp::export]]
double CCT_cpp(arma::vec & pval){
    double cauchyp, cctstat;
    uint np = pval.n_elem;    
    //check if there is NA	
    bool isna = pval.has_nan();
    if(isna){
        std::cout << "Cannot have NAs in the p-values!" << std::endl;
    }

    // check if all p-values are between 0 and 1

    arma::uvec between0and1Indice = arma::find(pval <= 1 && pval >= 0);
    if (between0and1Indice.n_elem != np){
        exit(EXIT_FAILURE); 
    }

    bool isnonzero = arma::all(pval);
    if(!isnonzero){
       cauchyp = 0;  
    }

    //if any of the p-value is one
    arma::uvec notoneIndice = arma::find(pval != 1.0);
    if(notoneIndice.n_elem < np){
	double minpval = pval.min();
	minpval = minpval * np;
	cauchyp = std::min(1.0, minpval);
    }else{	    
    double weight = 1.0/double(np);
    arma::vec weights(np, arma::fill::zeros);
    weights.replace(0, weight);

    arma::uvec issmallIndice = arma::find(pval < 1e-16);
    arma::uvec isNotSmallIndice = arma::find(pval >= 1e-16);
    if(issmallIndice.n_elem == 0){


	arma::vec tantemp = arma::tan((0.5 - pval) * M_PI);
	arma::vec weighttantemp = tantemp % weights;
        cctstat = arma::sum(weighttantemp);
    }else{
       arma::vec weightsub = weights.elem(issmallIndice);
       arma::vec pvalsub = pval.elem(issmallIndice);
       cctstat = arma::sum((weightsub/pvalsub)/M_PI);
       arma::vec pvalsub_b = pval.elem(isNotSmallIndice);
       arma::vec weightsub_b = weights.elem(isNotSmallIndice); 
       arma::vec tantemp_b = (0.5 - pvalsub_b) * M_PI;	
	tantemp_b = arma::tan(tantemp_b);
	arma::vec weighttantemp_b = tantemp_b % weightsub_b;
	cctstat = cctstat + arma::sum(weighttantemp_b);
    }

    if(cctstat > 1e+15){
       cauchyp = (1/cctstat)/M_PI;
    }else{
       boost::math::cauchy cauchy_dist(0,1);
       cauchyp = boost::math::cdf(complement(cauchy_dist, cctstat));
    }
    }//if(notoneIndice.n_elem < np){ 
    return(cauchyp);    
}




	
