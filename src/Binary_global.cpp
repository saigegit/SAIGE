#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cmath>
#include "Binary_ComputeExact.h"
#include "Binary_ComputeExactMC.h"
#include "Binary_HyperGeo.h"
#include "Binary_global.h"
#include "Binary_Permu_SKAT.h"

#ifdef _STAND_ALONE_
    #define Rprintf printf
#else
    #include <R.h>
    #include <Rmath.h>
#endif

/***************************************************
 
 Global functions
 
 ***************************************************/

void * SL_calloc(size_t num, size_t size){
	
	void * re =calloc(num, size);
	if(re == NULL){
		/*error("memory allocation error!");*/
        Rprintf("memory allocation error!");
        return(re) ;
	}
	
	return(re);
}


void SL_free(void * ptr){
	
	if(ptr != NULL){
		free(ptr);
        ptr=NULL;
	}
    
}

void SL_Sample(int k, int n, int *y, int *x){
    
    int i, j;
    
    for (i = 0; i < n; i++){
        x[i] = i;
    }
    for (i = 0; i < k; i++) {
        
        //temp = rand();
        //j = n * temp / (RAND_MAX+1);
        //j = temp % n;
        j = SL_runif_INT(n);
        y[i] = x[j] ;
        x[j] = x[--n];
    }
}


/*********************************************************
    Interface
 *********************************************************/

void SKAT_Exact(int * resarray, int nres, int * nres_k, double * Z0, double *Z1, int k, int m, int total, int * total_k, double *prob_k, double * odds, double * p1, int * IsExact, double * pval, double *pval_same, double *minP, int test_type, double epsilon)
{
    
    class ComputeExact exact;
    
    exact.Init(resarray, nres, nres_k, Z0, Z1, k, m, total, total_k, prob_k, odds, p1, IsExact, epsilon);
    exact.Run(test_type);
    //exact.PrintPval();
    exact.GetPvalues(pval, pval_same, prob_k, minP);
    
    
}


void SKATO_Exact(int * resarray, int nres, int * nres_k, double * Z0, double *Z1, double * r_corr, int n_r,
                 double * param, int k, int m, int total, int * total_k, double *prob_k, double * odds, double * p1, int * IsExact, double * pval, double *pval_same, double *minP, int test_type, double epsilon)
{
    
    
    class ComputeExactSKATO exact;
    
    
    exact.Init(resarray, nres, nres_k, Z0, Z1, r_corr, n_r, param, k,  m,  total,  total_k, prob_k,  odds, p1, IsExact, epsilon);

    
    exact.Run(test_type);
    exact.GetPvalues(pval, pval_same, prob_k, minP);
    
    //exact.PrintPval();    
}

/*
void SKAT_ExactMC(double * Q, int nQ, int ncohort, double * aZ0, double * aZ1, int * ak, int * am
                  , int total_m, int * amarker_idx, int * atotal, int * atotal_k, double * aprob_k
                  , double * ap1, int * aIsExact, double * pval, double *pval_same)
{
    
    class ComputeExactMC exact;
    
    
    exact.Init(Q, nQ, ncohort, aZ0, aZ1, ak, am, total_m, amarker_idx, atotal, atotal_k, aprob_k, ap1, aIsExact);
    exact.Run();
    exact.GetPvalues(pval, pval_same);
    
    
}

*/

void   GetProb(int k, int ngroup, int ncase, int * group, double * weight, double * prob){
    
    HyperGeo geo;
    geo.Run(k,ngroup,ncase, group, weight);
    geo.Get_lprob(prob);
    
}


void SKAT_Permu(double *Z, int *Y, int nSNP, int nSample, int nPermu, double * pval, double * pval_same, double epsilon)
{
    
    class Binary_Permu_SKAT permu1;
    
    permu1.Init(Z, Y, nSNP, nSample, nPermu, epsilon);

    permu1.Run();
    //exact.PrintPval();
    permu1.GetPvalues(pval, pval_same);
    
    
}

/*
extern "C" {
    
 
    void RSKATExactMC(double * Q, int * nQ, int * ncohort, double * aZ0, double * aZ1, int * ak, int * am
                      , int * total_m, int * amarker_idx, int * atotal, int * atotal_k, double * aprob_k
                      , double * ap1, int * aIsExact, double * pval, double *pval_same)
    {
        
        
        SKAT_ExactMC( Q, nQ[0], ncohort[0], aZ0, aZ1, ak, am
                     , total_m[0], amarker_idx, atotal, atotal_k, aprob_k
                     , ap1, aIsExact, pval, pval_same);
        
        
        
    }

    
    void RSKATExact(int * resarray, int * nres, int * nres_k, double * Z0, double *Z1, int * k, int * m, int * total, int * total_k, double *prob_k, double * odds, double * p1, int * IsExact, double * pval, double *pval_same, double *minP, int * test_type, double * epsilon)
    {
        
        SKAT_Exact(resarray, nres[0], nres_k, Z0, Z1, k[0], m[0], total[0], total_k, prob_k, odds, p1, IsExact, pval, pval_same, minP, test_type[0], epsilon[0]);
        
        
    }

    void RSKATOExact(int * resarray, int * nres, int * nres_k, double * Z0, double *Z1, double * r_corr, int * n_r, double * param, int * k, int * m, int * total, int * total_k, double *prob_k, double * odds, double * p1, int * IsExact, double * pval, double *pval_same, double *minP, int * test_type, double * epsilon)
    {
        
        SKATO_Exact(resarray, nres[0], nres_k, Z0, Z1, r_corr, n_r[0], param, k[0], m[0], total[0], total_k, prob_k, odds, p1, IsExact, pval, pval_same, minP, test_type[0], epsilon[0]);
        
        
    }
    
    void   RGetProb(int* k, int* ngroup, int* ncase, int * group, double * weight, double * prob){
        
        GetProb(k[0], ngroup[0], ncase[0], group, weight, prob);
    }

    
    void   RSKATPermu(double *Z, int *Y, int *nSNP, int *nSample, int *nPermu,double * pval, double *pval_same,  double * epsilon){
        
        SKAT_Permu(Z, Y, nSNP[0], nSample[0], nPermu[0], pval, pval_same, epsilon[0]);
        
    }
    
} // extern "C"
*/






