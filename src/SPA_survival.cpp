// [[Rcpp::depends(BH)]]

//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <iostream>
#include <fstream>
#include <cassert>
#include <stdexcept>
#include <memory>
#include <sstream>
#include <time.h>
#include <Rcpp.h>
#include <stdint.h>
#include <boost/math/distributions/normal.hpp>
#include "SPA_survival.hpp"
#include "UTIL.hpp"


// [[Rcpp::export]]
double Korg_Poi(double t1, arma::vec & mu, arma::vec & g)
{

	arma::vec temp0 = arma::exp(g * t1) - g * t1 - 1;
        double out = arma::dot(mu, temp0);
	return(out);
}


// [[Rcpp::export]]
double K1_adj_Poi(double t1, arma::vec & mu, arma::vec & g, double q)
{
	arma::vec temp1;
	arma::vec temp2;

	temp2 = mu % g;
	temp1 = temp2 % arma::exp(g * t1) - temp2;
	double out  = arma::sum(temp1) - q;
	return(out);
}


// [[Rcpp::export]]
double K2_Poi(double t1, arma::vec & mu, arma::vec & g)
{
        arma::vec temp0;
        arma::vec temp1;
        arma::vec temp2;
        arma::vec temp3;
		
	temp0 = arma::exp(g * t1);
	temp1 = arma::pow(g,2);
        temp2 = mu % temp1 % temp0;
        double out = sum_arma1(temp2);
        return(out);
}



// [[Rcpp::export]]
Rcpp::List getroot_K1_Poi(double init, arma::vec & mu, arma::vec & g, double q, double tol, int maxiter){
	Rcpp::List result;
	double root;
	int niter;
	bool Isconverge;
        int rep = 1;
	double K1_eval, K2_eval, t, tnew, newK1;
	double prevJump;
	
	t = init;
	K1_eval = K1_adj_Poi(t,mu,g,q);
	prevJump = std::numeric_limits<double>::infinity();
	bool conv = true;
	while(rep <= maxiter){
		K2_eval = K2_Poi(t,mu,g);
		tnew = t-K1_eval/K2_eval;
		if(tnew == NA_REAL){
			conv = false;
			break;
		}
			
		if(std::abs(tnew-t)<tol){	
			conv = true;
			break;	
		}
			
		if(rep == maxiter)
                {
                	conv = false;
                        break;
                }

		newK1 = K1_adj_Poi(tnew,mu,g,q);
                if(arma::sign(K1_eval) != arma::sign(newK1))
                {
                	if(std::abs(tnew-t) > (prevJump-tol))
                        {
                        	tnew = t + (arma::sign(newK1-K1_eval))*prevJump/2;
                                newK1 = K1_adj_Poi(tnew,mu,g,q);
                                prevJump = prevJump/2;
                         } else {
                                prevJump = std::abs(tnew-t);
                         }
                 }

		rep = rep + 1;
		t = tnew;
                K1_eval = newK1;
	}

	root=t;
	niter=rep;
	Isconverge=conv;	
	//std::cout << "here13" << std::endl;
	result["root"]	= root;
	result["niter"] = niter;
	result["Isconverge"] = Isconverge;
	return(result);
}



// [[Rcpp::export]]
Rcpp::List Get_Saddle_Prob_Poi(double zeta,  arma::vec & mu, arma::vec & g, double q, bool logp)
{
	double k1 = Korg_Poi(zeta, mu, g);
	double k2 = K2_Poi(zeta, mu, g);
	double temp1, w, v, Ztest, pval;
	double negative_infinity = - std::numeric_limits<double>::infinity();

	temp1 = zeta * q - k1;
	Rcpp::List result;
	bool isSaddle = false;
/*	std::cout << "k1 " << k1 << std::endl;
	std::cout << "k2 " << k2 << std::endl;
	std::cout << "temp1 " << temp1 << std::endl;
	std::cout << "zeta " << zeta << std::endl;
	std::cout << "q " << q << std::endl;
*/
	bool flagrun=false;
        if(std::isfinite(k1) && std::isfinite(k2) && temp1 >= 0 && k2 >= 0){
                 w = arma::sign(zeta) * std::sqrt(2 *temp1);
                 v = zeta *  std::sqrt(k2);
                 if(w != 0){
                        flagrun = true;
                 }
        }



	//if(std::isfinite(k1) && std::isfinite(k2) && temp1 > 0)
	if(flagrun)
	{
		//temp1 = zeta * q - k1;
		//w = arma::sign(zeta) * std::sqrt(2 *temp1);
		//v = zeta *  std::sqrt(k2);
		Ztest = w + (1/w) * std::log(v/w);

		/*std::cout << "w: " << w << std::endl;
		std::cout << "v: " << v << std::endl;
		std::cout << "Ztest: " << Ztest << std::endl;
*/
   		//boost::math::normal norm_dist(0,1);
		double pval0;
		

		if(Ztest > 0){
			//pval0 = boost::math::cdf(complement(norm_dist,Ztest));
			pval0 = R::pnorm(Ztest,0,1,false,logp);
			pval= pval0;
			//if(logp){
			//	pval = std::log(pval);
			//}	
		} else {
			//pval0 = boost::math::cdf(norm_dist,Ztest);
			pval0 = R::pnorm(Ztest,0,1,true,logp);
			pval= -1*pval0;
			
			//	-Rcpp::pnorm( Ztest, mean = 0.0, sd = 1.0, lower = true, log = logp );
		}
		isSaddle = true;
	} else {
		if(logp)
		{
			pval =  negative_infinity;
		}else {
			pval= 0;
		}
	}
	result["pval"] = pval;
	result["isSaddle"] = isSaddle;	
	return(result);
}



// [[Rcpp::export]]
Rcpp::List SPA_survival(arma::vec & mu, arma::vec & g, double q, double qinv, double pval_noadj, double tol, bool logp){
	using namespace Rcpp;
	List result;
	double p1, p2, pval;
	bool Isconverge = true;
	Rcpp::List outuni1 = getroot_K1_Poi(0, mu, g, q, tol);
	Rcpp::List outuni2 = getroot_K1_Poi(0, mu, g, qinv, tol);
	double outuni1root = outuni1["root"];
	double outuni2root = outuni2["root"];
	bool Isconverge1 = outuni1["Isconverge"];
	bool Isconverge2 = outuni2["Isconverge"];

	//std::cout << "outuni1root" << outuni1root << std::endl;
        //std::cout << "outuni2root" << outuni2root << std::endl;
        //std::cout << "Isconverge1" << Isconverge1 << std::endl;
        //std::cout << "Isconverge2" << Isconverge2 << std::endl;


	Rcpp::List getSaddle;
        Rcpp::List getSaddle2;
	if(outuni1["Isconverge"]  && outuni2["Isconverge"])
	{
		//std::cout << "q is " << q << " 3 qinv is " << qinv << std::endl;
		getSaddle = Get_Saddle_Prob_Poi(outuni1["root"], mu, g, q, logp);

		if(getSaddle["isSaddle"]){
			p1 = getSaddle["pval"];
		}else{	

			if(logp){
				p1 = pval_noadj-std::log(2);
			}else{
				p1 = pval_noadj/2;	
			}	
		}
		getSaddle2 = Get_Saddle_Prob_Poi(outuni2["root"], mu, g, qinv, logp);
		if(getSaddle2["isSaddle"]){
                        p2 = getSaddle2["pval"];
                }else{
			if(logp){
                                p2 = pval_noadj-std::log(2);
                        }else{
                                p2 = pval_noadj/2;
                        }
		}
		//std::cout << "p1_nofast " << p1 << "p2 " << p2 << std::endl;
		
		if(logp)
		{
			pval = add_logp(p1,p2);
		} else {
			pval = std::abs(p1)+std::abs(p2);
		}
		Isconverge=true;
	}else {
 			//std::cout << "Error_Converge" << std::endl;
			pval = pval_noadj;
			Isconverge=false;	
		}		
	result["pvalue"] = pval;
	result["Isconverge"] = Isconverge;
	return(result);
}


// [[Rcpp::export]]
double Korg_fast_Poi(double t1, arma::vec & mu, arma::vec & g,  arma::vec & gNA, arma::vec & gNB, arma::vec & muNA, arma::vec & muNB, double NAmu, double NAsigma)
{

    	arma::vec temp0 = arma::exp(gNB*t1) - gNB*t1 - 1; 
	arma::vec temp = muNB % temp0;

        double out = arma::sum(temp) + 0.5*NAsigma*pow(t1,2);
	return(out);
}


// [[Rcpp::export]]
double K1_adj_fast_Poi(double t1, arma::vec & mu, arma::vec & g, double q, arma::vec & gNA, arma::vec & gNB, arma::vec & muNA, arma::vec & muNB,  double NAmu, double NAsigma)
{

	arma::vec temp0;
	arma::vec temp1;
	double temp2;

	temp0 = muNB % gNB;
	temp1 = temp0 % arma::exp(gNB * t1) - temp0;
	temp2 = NAsigma*t1;
	double out  = arma::sum(temp1) + temp2 - q;
	return(out);
}



// [[Rcpp::export]]
double K2_fast_Poi(double t1, arma::vec & mu, arma::vec & g, arma::vec & gNA, arma::vec & gNB, arma::vec & muNA, arma::vec & muNB,  double NAmu, double NAsigma)
{

        arma::vec temp0;
        arma::vec temp1;
        arma::vec temp2;
	//std::cout << "t1 " << t1 << std::endl;		
	temp0 = arma::exp(gNB * t1);
       	temp1 = arma::pow(gNB,2) % temp0;			
        temp2 = muNB % temp1;
	//temp2.print("temp2");
	//std::cout << "NAsigma " << NAsigma << std::endl;		
	
        double out = sum_arma1(temp2)+NAsigma;
        return(out);
}





// [[Rcpp::export]]
Rcpp::List getroot_K1_fast_Poi(double init, arma::vec & mu, arma::vec & g, double q, arma::vec & gNA, arma::vec & gNB, arma::vec & muNA, arma::vec & muNB,  double NAmu, double NAsigma, double tol, int maxiter){
	Rcpp::List result;
	double root;
	int niter;
	bool Isconverge;
        int rep = 1;
	double K1_eval, K2_eval, t, tnew, newK1;
	double prevJump;

	t = init;
	K1_eval = K1_adj_fast_Poi(t,mu,g,q,gNA,gNB,muNA,muNB,NAmu, NAsigma);
	//std::cout << "K1_eval " << K1_eval << std::endl; 
	prevJump = std::numeric_limits<double>::infinity();
	bool conv = true;
	while(rep <= maxiter){
		K2_eval = K2_fast_Poi(t,mu,g, gNA,gNB,muNA,muNB,NAmu, NAsigma);
		//std::cout << "K2_eval " << K2_eval << std::endl; 
		tnew = t-K1_eval/K2_eval;
		if(tnew == NA_REAL){
			conv = false;
			break;

		}
			
		if(std::abs(tnew-t)<tol){	
			conv = true;
			break;	
		}
			
		if(rep == maxiter)
                {
                	conv = false;
                        break;
                }

		newK1 = K1_adj_fast_Poi(tnew,mu,g,q, gNA,gNB,muNA,muNB,NAmu, NAsigma);
                if((K1_eval * newK1) < 0)
                {
                	if(std::abs(tnew-t) > (prevJump-tol))
                        {
                        	tnew = t + (arma::sign(newK1-K1_eval))*prevJump/2;
                                newK1 = K1_adj_fast_Poi(tnew,mu,g,q, gNA,gNB,muNA,muNB,NAmu, NAsigma);
                                prevJump = prevJump/2;
                        } else {
                                prevJump = std::abs(tnew-t);
                        }
                }

		rep = rep + 1;
		t = tnew;
                K1_eval = newK1;
	}
	root=t;
	niter=rep;
	Isconverge=conv;	
	result["root"]	= root;
	result["niter"] = niter;
	result["Isconverge"] = Isconverge;
	return(result);
}



// [[Rcpp::export]]
Rcpp::List Get_Saddle_Prob_fast_Poi(double zeta,  arma::vec & mu, arma::vec & g, double q, arma::vec & gNA, arma::vec & gNB, arma::vec & muNA, arma::vec & muNB,  double NAmu, double NAsigma, bool logp)
{
	double k1 = Korg_fast_Poi(zeta, mu, g,  gNA,gNB,muNA,muNB,NAmu, NAsigma);
	double k2 = K2_fast_Poi(zeta, mu, g,  gNA,gNB,muNA,muNB,NAmu, NAsigma);
	double temp1, w, v, Ztest, pval;
	double negative_infinity = - std::numeric_limits<double>::infinity();
	temp1 = zeta * q - k1;
	Rcpp::List result;
	bool isSaddle = false; 

	bool flagrun=false;
        if(std::isfinite(k1) && std::isfinite(k2) && temp1 >= 0 && k2 >= 0){
                 w = arma::sign(zeta) * std::sqrt(2 *temp1);
                 v = zeta *  std::sqrt(k2);
                 if(w != 0){
                        flagrun = true;
                 }
        }



	//if(std::isfinite(k1) && std::isfinite(k2) && temp1 > 0 )
	if(flagrun)
	{
		//temp1 = zeta * q - k1;
		//w = arma::sign(zeta) * std::sqrt(2 *temp1);
		//v = zeta *  std::sqrt(k2);
		Ztest = w + (1/w) * std::log(v/w);
	
		boost::math::normal norm_dist(0,1);
                double pval0;
		if(Ztest > 0){
			//pval0 = boost::math::cdf(complement(norm_dist,Ztest));
			pval0 = R::pnorm(Ztest,0,1,false,logp);
			pval=pval0;
			//if(logp){
			//	pval = std::log(pval);
			//}	
		} else {
			//pval0 = boost::math::cdf(norm_dist,Ztest);
			pval0 = R::pnorm(Ztest,0,1,true,logp);
			pval= -pval0;
			
			//	-Rcpp::pnorm( Ztest, mean = 0.0, sd = 1.0, lower = true, log = logp );
		}
		isSaddle = true;
		} else {
			if(logp)
			{
				pval =  negative_infinity;
			}else {
				pval=0;
			}
		}
	result["pval"] = pval;
        result["isSaddle"] = isSaddle;
        return(result);


}



// [[Rcpp::export]]
Rcpp::List SPA_survival_fast(arma::vec & mu, arma::vec & g, double q, double qinv, double pval_noadj, bool logp, arma::vec & gNA, arma::vec & gNB, arma::vec & muNA, arma::vec & muNB,  double NAmu, double NAsigma, double tol){
	using namespace Rcpp;
	List result;
	double p1, p2, pval;
	bool Isconverge = true;
	Rcpp::List outuni1 = getroot_K1_fast_Poi(0, mu, g, q, gNA,gNB,muNA,muNB,NAmu, NAsigma, tol);
	//double qinv = -1*q;
	Rcpp::List outuni2 = getroot_K1_fast_Poi(0, mu, g, qinv, gNA,gNB,muNA,muNB,NAmu, NAsigma, tol);
	//std::cout << "outuni1root" << outuni1["root"] << std::endl;
	//std::cout << "outuni2root" << outuni2["root"] << std::endl;

	Rcpp::List getSaddle;
	Rcpp::List getSaddle2;
	if(outuni1["Isconverge"]  && outuni2["Isconverge"])
	{
		
		getSaddle  = Get_Saddle_Prob_fast_Poi(outuni1["root"], mu, g, q, gNA,gNB,muNA,muNB,NAmu, NAsigma, logp);
		if(getSaddle["isSaddle"]){
			p1 = getSaddle["pval"];
		}else{	

		        if(logp){
				p1 = pval_noadj-std::log(2);
			}else{
				p1 = pval_noadj/2;	
			}	
		}

		getSaddle2 = Get_Saddle_Prob_fast_Poi(outuni2["root"], mu, g, qinv, gNA,gNB,muNA,muNB,NAmu, NAsigma, logp);
		if(getSaddle2["isSaddle"]){
			p2 = getSaddle2["pval"];
		
		}else{
			if(logp){
				p2 = pval_noadj-std::log(2);
			}else{
				p2 = pval_noadj/2;	
			}	
		}

		if(logp){
			pval = add_logp(p1,p2);
		} else {
			pval = std::abs(p1)+std::abs(p2);
			//std::cout << "p1 " << p1 << "p2 " << p2 << std::endl; 
		}
		Isconverge=true;
	}else {
 			//std::cout << "Error_Converge" << std::endl;
			pval = pval_noadj;
			Isconverge=false;	
		}		
	result["pvalue"] = pval;
	result["Isconverge"] = Isconverge;
	return(result);
}

