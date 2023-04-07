#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>

#include "Binary_ComputeExact.hpp"
#include "Binary_global.hpp"

#ifdef _STAND_ALONE_
#define Rprintf printf
#else
#include <R.h>
#include <Rmath.h>
#endif

/***************************************************

   Protected functions
 
***************************************************/


double     ComputeExactSKATO::Cal_Pvalue_Rcorr(double test_skat, double test_C, int idx){
    
    double pval=0;
    double Qnorm= 0;
    double rho = m_rcorr[idx];
    double muQ = m_muQ[idx];
    double varQ = m_varQ[idx];
    double df = m_df[idx];
    
    double Q = (1-rho) * test_skat + rho * test_C ;
    //Rprintf("Param: %d-[%f][%f][%f]\n",idx, Q, test_skat,test_C );
    
#ifdef _STAND_ALONE_
    // if there is no chisq function, just normalize it
    pval = (muQ  - Q)/sqrt(varQ);
#else

    pval=1;
    
    if(varQ > 0){

        Qnorm = (Q  - muQ)/sqrt(varQ) * sqrt(2 * df) + df;
        pval = pchisq(Qnorm, df, 0,1);
    }

#endif
    //Rprintf("Param1: %d-[%f][%f][%f][%f][%f][%f][%f][%f]\n",idx, Qnorm, pval,Q, muQ, varQ, df, test_skat,test_C );
    
    return pval;
}


// int     ComputeExactSKATO::CalTestStat(int k, int * array){
double ComputeExactSKATO::CalTestStat
(int k, int * array, bool is_save, bool is_minIdx, int * minIdx ){
    

	int i, j, l, temp;
    double test_C  = m_teststat_Z0_C;
    double test_skat = 0;
    double teststat;
    memcpy(m_teststat_one.data(), m_teststat_Z0.data(), sizeof(double) *m_m);
    
	
    for(i=0;i< k;i++){
        l = array[i];
        temp = l*m_m;
        for(j=0;j< m_m;j++){
            m_teststat_one[j]+=m_Z1[temp+j] - m_Z0[temp+j] ;
        }
        test_C += m_Z1_C[l] - m_Z0_C[l] ;
    }
    
	for(j=0;j< m_m;j++){
		test_skat+=m_teststat_one[j] * m_teststat_one[j];
	}
    test_C = test_C * test_C;
    
    for(i=0; i<m_rcorr.size(); i++){
        
        double stat = Cal_Pvalue_Rcorr(test_skat, test_C, i);
        /*if(!is_save){
            Rprintf("%d-[%d][%e][%e][%e]\n",k,i, stat, test_skat, test_C);
        }*/

        if(i==0){
            teststat = stat;
        } else {
            
            teststat = fmin(teststat, stat);
        }
        
        if(is_minIdx){
            if(teststat == stat){
                *minIdx = i;
            }
        }
        
    }
    
    /* since larger Q has a smaller p-value, it should be -teststat, not teststat */
    if(is_save){
        m_teststat[m_idx] = -teststat;  
        
    }
/*    
    Rprintf("[%d]-",k);
    for(i=0;i< k;i++){
        l = array[i];
        Rprintf("[%d]",l);
    }
    Rprintf("=[%e]",-teststat);
    if(is_save){
        Rprintf("S\n");
    } else {
        Rprintf("N\n");
    }
*/    
    return -teststat;
    
}

//int      ComputeExactSKATO::CalTestStat_INV(int k, int * array)
double      ComputeExactSKATO::CalTestStat_INV(int k, int * array, bool is_save, bool is_minIdx, int * minIdx ){
    
    
	int i, j, l, temp;
    double test_C  = m_teststat_Z1_C;
    double test_skat = 0;
    double teststat;
    // There are easier ways to copy vectors
    memcpy(m_teststat_one.data(), m_teststat_Z1.data(), sizeof(double) *m_m);
    
    for(i=0;i< k;i++){
        l = array[i];
        temp = l*m_m;
        for(j=0;j< m_m;j++){
            m_teststat_one[j]+=m_Z0[temp+j] - m_Z1[temp+j] ;
        }
        test_C += m_Z0_C[l] - m_Z1_C[l] ;
    }
    test_C = test_C * test_C;
    
	for(j=0;j< m_m;j++){
		test_skat+=m_teststat_one[j] * m_teststat_one[j];
	}
    
    for(i=0; i<m_rcorr.size(); i++){
        
        double stat = Cal_Pvalue_Rcorr(test_skat, test_C, i);
        if(i==0){
            teststat = stat;
        } else {
            
            teststat = fmin(teststat, stat);
        }
        if(is_minIdx){
            if(teststat == stat){
                *minIdx = i;
            }
        }

    }
    
    if(is_save){
        m_teststat[m_idx] = -teststat;        
    }
    
    
/*    
    Rprintf("[%d]-",k);
    for(i=0;i< k;i++){
        l = array[i];
        Rprintf("[%d]",l);
    }
    Rprintf("=[%e]",-teststat);
    if(is_save){
        Rprintf("S\n");
    } else {
        Rprintf("N\n");
    }
*/

    return -teststat;
    
}


/***************************************************
 
    Public function
 
 **************************************************/


ComputeExactSKATO::ComputeExactSKATO(){
    
    m_pprod=1;

}
/*
ComputeExactSKATO::~ComputeExactSKATO(){
    
    

}
*/
/************************************************
*
*	Init function
*
*	Q : Test statistic
*	Z0 : Z * (0 - muhat) matrix
*	Z1 : Z * (1 - muhat) matrix
*	k : number of samples with at least one minor allele
*	m : 
*	total : total possible combination
*	total_k : possible combination for each k (k=0,..., K)
*	prob_k : probability to be K=k
*	p1 : prob to be case
*	pval : outcome p-value
*
******************************************************/

int     ComputeExactSKATO::Init(int * resarray, int nres, int * nres_k, double * Z0, double *Z1, double * r_corr, int n_r, double * param, int k, int m, int total, int * total_k, double *prob_k, double * odds, double * p1, int * IsExact, double epsilon, bool IsSmallmemory){

    int i,j;
    SaveParam(Z0, Z1, k, m, total, total_k, prob_k, odds, p1, IsExact, epsilon, IsSmallmemory);
    
    
    /* For SKAT-O */
    if(n_r > 0){
        for(i=0;i<n_r;i++){
            m_rcorr.push_back(r_corr[i]);
        }
        m_teststat_Z0_C = 0;
        m_teststat_Z1_C = 0;
        
        m_Z0_C = (double *) SL_calloc(m_k * 1, sizeof(double));
        m_Z1_C = (double *) SL_calloc(m_k * 1, sizeof(double));	
        memset(m_Z0_C,0, sizeof(double)*m_k);
        memset(m_Z1_C,0, sizeof(double)*m_k);
        
        for(i=0;i<m_k;i++){
            int idx = i*m_m;
            for(j=0;j<m_m;j++){
                m_Z0_C[i]+=m_Z0[idx+j];
                m_Z1_C[i]+=m_Z1[idx+j];
            }
            
        }
        for(i=0;i<m_k;i++){
            m_teststat_Z0_C += m_Z0_C[i];
            m_teststat_Z1_C += m_Z1_C[i];
        }
        
    }
    
    int idx;
    for(i=0;i<n_r;i++){
        for(j=0;j<3; j++){
            
            idx = i*3 + j;
            if(j==0){
                m_muQ.push_back(param[idx]);
            } else if (j==1){
                m_varQ.push_back(param[idx]);
            } else if (j==2){
                m_df.push_back(param[idx]);
            }
            
            //Rprintf("%d-[%f]\n",j, param[idx]);
            
        }
    }
	
    idx = 0;
    int minIdx;
    for(i=0;i<nres;i++){
        
        int * array = (resarray + idx) ;
        int k1 = nres_k[i];
        idx += k1;
        
        double Q = CalTestStat(k1, array, false, true, &minIdx);
        m_Q.push_back(Q);
        m_minIdx.push_back(minIdx);
        
        //Rprintf("[%d][%e]\n", i, Q);
    }
    
    return 1;
}
