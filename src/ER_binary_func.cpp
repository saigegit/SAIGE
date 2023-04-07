//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cmath>
#include "Binary_HyperGeo.hpp"
#include "Binary_ComputeExact.hpp"
#include "Binary_global.hpp"

//https://github.com/cran/SKAT/blob/ebad53e501f67070f4858c76eb420d9e0ba28219/src/Binary_global.cpp#L118

void GetProb_new(int k, int ngroup, int ncase, int* group, double* weight, double* prob){    
    HyperGeo geo;
    geo.Run(k,ngroup,ncase, group, weight);
    geo.Get_lprob(prob);
    
}


void SKATExactBin_ComputeProb_Group(arma::uvec & idx, arma::uvec & idxCompVec, arma::vec & pi1, uint32_t n, uint32_t ncase, int type_group, std::vector<double> & prob){

	uint32_t k = idx.n_elem;
	int ngroup1 = 10;  //use the default value as ER will only be used for variants with MAC <= 10;
	arma::vec p1 = pi1(idx);
	arma::vec p2 = pi1(idxCompVec);
		
	arma::uvec id_temp = arma::find(p1 >= 1);
	
	if(id_temp.n_elem > 0){
		for(unsigned int j = 0; j < id_temp.n_elem; j++){
			unsigned int id_temp_j = id_temp(j);
			p1(id_temp_j) = 0.999;
		}
	}

	std::vector<double> weight;
	std::vector<int> group;
	arma::uvec a1Vec, a2Vec, IDX;
	double a1, a2, p1temp, p2temp, oddtemp, p2oddtemp;
	arma::vec p1tempVec;

	for(int i = 0; i < ngroup1; i++){
		a1 = double(i) / ngroup1;
		a2 = double(i + 1) / ngroup1;
		if((i+1) < ngroup1){
			a1Vec = arma::find(p1 >= a1);				
			a2Vec = arma::find(p1 < a2);				
		}else{	
			a1Vec = arma::find(p1 >= a1);				
			a2Vec = arma::find(p1 <= a2);					
		}
		IDX = arma::intersect(a1Vec, a2Vec);
		
		if(IDX.n_elem > 0){
			p1tempVec = p1(IDX);
			p1temp = arma::mean(p1tempVec);
			oddtemp = p1temp/(1-p1temp);
			weight.push_back(oddtemp);
			group.push_back(IDX.n_elem);
		}	
        }
	p2temp = arma::mean(p2);
	p2oddtemp = p2temp / (1-p2temp);
	weight.push_back(p2oddtemp);

	for(int i = 0; i < weight.size(); i++){	
		weight[i] = weight[i] / p2oddtemp;
	}
	group.push_back(n-k);
	//std::vector<double> prob_k(k+1, 0.0);
	int ngroup = group.size();
	int ncasei = int(ncase);
	GetProb_new(k, ngroup, ncasei, &group[0], &weight[0], &prob[0]);
}


int fact(int n) {
	if (n == 0 || n == 1){
   		return 1;
	}else{
   		return n * fact(n - 1);
   	}
}

int n_choose_r(int n, int r){
	int comb;
	comb = fact(n) / (fact(r) * fact(n-r));	
	return(comb);
}


void Get_Total_K(int k, std::vector<int> & n_total_k){
	//std::vector<int> n_total_k(k+1, 0);
	for(int i = 0; i <= k; i++){
		n_total_k[i] = n_choose_r(k, i);
	}	
	//Rcpp::List re;
	//int ntotal = std::accumulate(n_total_k.begin(),n_total_k.end(),0);
	//re["ntotal"] = ntotal;
	//re["n_total_k"] = n_total_k;
        //return(re);
}

void SKATExactBin_ComputProb_New(arma::uvec & idx, arma::uvec & idxCompVec, arma::vec & pi1, uint32_t n, uint32_t ncase, int NResampling, int ExactMax, int test_type, int type_group, std::vector<double> & prob, std::vector<int> & IsExactVec, std::vector<int> & n_total_k, int & n_total, bool & Is_ExactP){
	
        uint32_t k = idx.n_elem;
        //arma::vec p1 = pi1(idx);

	//std::vector<double> prob(k+1, 0.0);
	SKATExactBin_ComputeProb_Group(idx, idxCompVec, pi1, n, ncase, type_group, prob);


	//std::vector<int> IsExactVec(k+1, 1);
		
	//Rcpp::List obj_total = Get_Total_K(k);
	Get_Total_K(k, n_total_k);
	n_total = std::accumulate(n_total_k.begin(),n_total_k.end(),0);
	//std::vector<int> n_total_k = obj_total["n_total_k"];
	Is_ExactP = true;
	if(n_total > NResampling){
		for(int i = 0; i <= k; i++){
			if(n_total_k[i] > ExactMax){
				n_total_k[i] = int(ceil(NResampling * prob[i]));
				IsExactVec[i] = 0;
			}	
		}
		Is_ExactP = false;
	}	
	n_total = std::accumulate(n_total_k.begin(),n_total_k.end(),0);

	//Rcpp::List re;
        //re["prob"] = prob;
        //re["k"] = k;
        //re["n"] = n;
        //re["ntotal"] = n_total;
        //re["n_total_k"] = n_total_k;
        //re["IsExact"] = IsExactVec;
        //re["p1"] = p1;
        //re["Is_ExactP"] = Is_ExactP;

	//return(re);
}	



void  Get_Res_Arrays(arma::mat & res_out, arma::uvec & idx, std::vector<int> & resarray, int & nres, std::vector<int> & nres_k){
	//nres = res_out.n_cols;
	//std::vector<int> nres_k(nres, 0);
	arma::vec res_out_colvec;
	//std::vector<int> resarray;
	for(int i = 0; i < nres; i++){
		 res_out_colvec = res_out.col(i);
		 arma::uvec res_out_i = arma::find( res_out_colvec > 0);
		 arma::uvec res_out_i_s = arma::sort(res_out_i);
		 int res_out_i_s_k = res_out_i_s.n_elem;
		 nres_k[i] = res_out_i_s.n_elem;
		 for(int k = 0; k < res_out_i_s_k; k++){
			resarray.push_back(res_out_i_s(k));
		 }	 
	}

	//Rcpp::List re;
	//re["resarray"] = resarray;
	//re["nres"] = nres;
	//re["nres_k"] = nres_k;
	//return(re);
}


std::vector< std::vector<double> > mat_to_std_vec(arma::mat &A) { 

    std::vector< std::vector<double> > V(A.n_rows);
    typedef std::vector<double> stdvec;
    for (size_t i = 0; i < A.n_rows; ++i) {
        V[i] = arma::conv_to< stdvec >::from(A.row(i));
    };
    return V;


}


double SKATExactBin_Work(arma::mat & Z, arma::vec & res, arma::vec & pi1, uint32_t ncase, arma::uvec & idx, arma::uvec & idxCompVec, arma::mat & res_out, int NResampling, int ExactMax, double epsilon, int test_type){

        uint32_t n = res.n_elem;
        arma::vec p1 = pi1(idx);
        arma::vec p2 = pi1(idxCompVec);

	arma::mat Z_1 = Z.rows(idx);
	arma::mat Z1temp = (Z_1 % (-p1)).t();
	arma::vec Z0 = arma::vectorise(Z1temp);

	arma::mat Z1temp2 = (Z_1 % (1-p1)).t();
	arma::vec Z1 =  arma::vectorise(Z1temp2);

	uint32_t m = Z_1.n_cols;
	int k = idx.n_elem;
	std::vector<int> n_total_k(k+1, 0);

	std::vector<double> prob(k+1, 0.0);
	std::vector<int> IsExactVec(k+1, 1);

	int n_total = 0;
	bool Is_ExactP = false;
	SKATExactBin_ComputProb_New(idx, idxCompVec, pi1, n, ncase, NResampling, ExactMax, test_type, 2, prob, IsExactVec, n_total_k, n_total, Is_ExactP); 
			
	
	double p1mean = arma::mean(p1);
	arma::vec p1_adj = p1 / p1mean;
	arma::vec odds = p1 / (1 - p1);
	int test_type_new=1;

	if(res_out.is_empty()){
		res_out = res(idx);
        }else{
		arma::vec res_out2 = arma::join_cols(res(idx), res_out(idx));
		res_out.resize(res_out2.n_elem);
		res_out = res_out2;
	}	


	std::vector<int> resarray;
	int nres = res_out.n_cols;
	std::vector<int> nres_k(nres, 0);
	
	Get_Res_Arrays(res_out, idx, resarray, nres, nres_k);
	std::vector<double> pval(nres, 0.0); 
	std::vector<double> pval1(nres, 0.0); 
	double minP = 100;
	
	typedef std::vector<double> stdvec;
	stdvec Z1std = arma::conv_to< stdvec >::from(Z1);
	stdvec Z0std = arma::conv_to< stdvec >::from(Z0);
	stdvec oddsstd = arma::conv_to< stdvec >::from(odds);
	stdvec p1_adjstd = arma::conv_to< stdvec >::from(p1_adj);

	typedef std::vector<int> istdvec;
	//istdvec nres_kstd = arma::conv_to< istdvec >::from(re_arr["nres_k"]);	
	//istdvec n_total_kstd = arma::conv_to< istdvec >::from(pr["n_total_k"]);
	//istdvec IsExactstd = arma::conv_to< istdvec >::from(pr["IsExact"]);


	SKAT_Exact(resarray, nres, &nres_k[0], &Z0std[0], &Z1std[0], k, m, n_total, &n_total_k[0], &prob[0], &oddsstd[0], &p1_adjstd[0], &IsExactVec[0], &pval[0], &pval1[0], &minP, test_type_new, epsilon);

	//arma::mat pval_re(pval.size(),2);

	//arma::vec pval_b =  arma::conv_to< arma::vec >::from(pval);
	//arma::vec pval1_b =  arma::conv_to< arma::vec >::from(pval1);

	//pval_re.col(0) = pval_b;
	//pval_re.col(1) = pval1_b;

	//minP = minP/2;
	//arma::vec prob_k_b =  arma::conv_to< arma::vec >::from(probstd);


	//arma::mat pval_re1(pval.size(),3);
	//pval_re1.col(0) = pval_b;
	//pval_re1.col(1) = pval_b - pval1_b/2;
	//pval_re1.col(3) = pval1_b;


  	//if(!pr["Is_ExactP"]){
        //        minP=-1
        //}

        //p.value.resampling=NULL
        //p.value.standard.resampling=NULL
        //Q.resampling=NULL

	//std::cout << "pval[0] " << pval[0] << std::endl;
	//std::cout << "pval1[0] " << pval1[0]/2 << std::endl;
	double pvalue = pval[0] - pval1[0]/2;
	//std::cout << "pvalue " << pvalue << std::endl;
	return(pvalue);	
}
