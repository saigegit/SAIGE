#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>

#include "Binary_ComputeExactMC.h"
#include "Binary_global.h"

#ifdef _STAND_ALONE_
#define Rprintf printf
#else
#include <R.h>
#include <Rmath.h>
#endif


/*****************************************************
 Cohort Info
 ***************************************************/



CohortInfo::CohortInfo(){
    
    m_fprob=NULL;  
    m_Z0=NULL;
    m_Z1=NULL;
    
    m_teststat_one=NULL;
    m_teststat_Z0=NULL;
    m_teststat_Z1=NULL;
    m_teststat_all=NULL;
    
    m_pprod=1;
    m_fprob=NULL;
    m_array=NULL;
    
}

CohortInfo::~CohortInfo(){

	
	SL_free(m_fprob);
	SL_free(m_teststat_all);

    
	SL_free(m_Z0);
	SL_free(m_Z1);
	SL_free(m_teststat_Z0);
    SL_free(m_teststat_Z1);
	SL_free(m_teststat_one);
    
    SL_free(m_array);
   
    
}



int     CohortInfo::CalTestStat(int k, int * array, int is_case){
    
	int i, j, l, temp;
    
    if(is_case == 1){
        memcpy(m_teststat_one, m_teststat_Z0, sizeof(double) *m_m);
	
        for(i=0;i< k;i++){
            l = array[i];
            temp = l*m_m;
            for(j=0;j< m_m;j++){
                m_teststat_one[j]+=m_Z1[temp+j] - m_Z0[temp+j] ;
            }
        }
    } else {
        memcpy(m_teststat_one, m_teststat_Z1, sizeof(double) *m_m);
        
        for(i=0;i< k;i++){
            l = array[i];
            temp = l*m_m;
            for(j=0;j< m_m;j++){
                m_teststat_one[j]+=m_Z0[temp+j] - m_Z1[temp+j] ;
            }
        }
    }
    
    i = (m_idx ) * m_m;
    memcpy((m_teststat_all + i), m_teststat_one, sizeof(double) *m_m); 
    return 0;
    
}

int     CohortInfo::CalFisherProb(int k, int * array, int is_case){
    
	int i,l,k1;
	double temp = 1;
    
    if(is_case == 1){
        temp = 1;
        k1=k;
        for(i=0;i< k;i++){
            l = array[i];
            temp = temp * m_p1[l];
        }
    } else {
        temp = m_pprod;
        k1 = m_k - k;
        for(i=0;i< k;i++){
            l = array[i];
            temp = temp / m_p1[l];
        }        
    }
	m_fprob[m_idx] = temp;
    m_denomi[k1] = m_denomi[k1]+temp;
    
    return 0;
    
}



int     CohortInfo::Exact_Recurse(int k, int * array, int cell, int start, int end, int is_case){
    
	int i;
	
	/* node */
	if(k == cell){
		CalTestStat(k, array,is_case);
		CalFisherProb(k, array,is_case);
		m_idx++;
	} else {
        for(i=start; i<end; i++) {
            array[cell] = i;
            Exact_Recurse(k, array, cell+1, i+1, end, is_case);
        }
        
	}
    return 0;
}




int      CohortInfo::Init(double * Z0, double *Z1, int k, int m, int total_m, int * marker_idx, int total
              , int * total_k, double *prob_k, double * p1, int * IsExact){
    
	int i,j;
	m_idx=0;
	//Rprintf("0\n");
	
	m_k = k;
	m_m = m;
    m_total_m = total_m;
    
	m_total = total;
    
    // Init _k vectors
    m_pprod=1;
    for(i=0;i<=k;i++){
        
        m_total_k.push_back(total_k[i]);
        m_prob_k.push_back(prob_k[i]);
        m_denomi.push_back(0);
        m_IsExact.push_back(IsExact[i]);
        
        // # of samples = k so i should be smaller than k.
        if(i<k){
            m_p1.push_back(p1[i]);
            m_pprod = m_pprod*p1[i];
        }
        
    }
    for(i=0;i<m;i++){
        
        m_marker_idx.push_back(marker_idx[i]);
    }
	
	m_Z0 = (double *) SL_calloc(m_k * m_m, sizeof(double));
	m_Z1 = (double *) SL_calloc(m_k * m_m, sizeof(double));	
	m_teststat_Z0 = (double *) SL_calloc(m_m, sizeof(double));
    m_teststat_Z1 = (double *) SL_calloc(m_m, sizeof(double));

    
    
    m_array = (int *) SL_calloc(m_k+1, sizeof(int));
	
	memcpy(m_Z0, Z0, sizeof(double) * m_k * m_m);
	memcpy(m_Z1, Z1, sizeof(double) * m_k * m_m);
    
    memset(m_teststat_Z0,0, sizeof(double)*m_m);
    memset(m_teststat_Z1,0, sizeof(double)*m_m);
	/* generate prob matrix */
	for(i=0;i< m_k;i++){
		int idx = i*m_m;
		for(j=0;j< m_m;j++){
			m_teststat_Z0[j]+=m_Z0[idx+j];
            m_teststat_Z1[j]+=m_Z1[idx+j];
		}
	}
    
	m_fprob = (double *) SL_calloc(m_total, sizeof(double));
	m_teststat_one = (double *) SL_calloc(m_m, sizeof(double));
    m_teststat_all = (double *)SL_calloc(m_m * m_total, sizeof(double));
    
    
    m_idx = 0;
    int idx = 0;
    for(i=0;i< m_k+1;i++){
        if(i <= m_k/2 +1){
            m_combi_case.push_back(1);
            Exact_Recurse(i, m_array, 0, 0, m_k, 1);
        } else {
            m_combi_case.push_back(0);
            Exact_Recurse(m_k - i, m_array, 0, 0, m_k, 0);
        }
        
        for(j=idx;j< idx + m_total_k[i];j++){
            m_fprob[j] = m_fprob[j] / m_denomi[i] * m_prob_k[i];
        }
        idx = idx + m_total_k[i];
    }
    
    
    return 1;
}

//
//  idx start from 0
int    CohortInfo::Sum_TestStat(int idx, double* teststat){
    
    if(idx >= m_total){
        return -1;
    }

    int i, l;
    int idx1; 

    // m_marker_idx start from 1
    for(i=0;i<m_m;i++){
        l = m_marker_idx[i]-1;
        idx1 = idx * m_m;
        teststat[l]+=m_teststat_all[idx1 + i] ;
    }
    return 1;
}

int    CohortInfo::Delete_TestStat(int idx, double* teststat){
    
    if(idx >= m_total){
        return -1;
    }
    if(idx < 0){
        return 1;
    }
    
    int i, l;
    int idx1; 
    
    for(i=0;i<m_m;i++){
        
        l = m_marker_idx[i]-1;
        idx1 = idx * m_m;
        teststat[l]-=m_teststat_all[idx1 + i] ;
    }
    return 1;
}


double  CohortInfo::GetProb(int idx){
    
    
    if(idx >= m_total){
        return -1;
    }
    
    return m_fprob[idx];
    
}


/***************************************************
 
    Public function
 
 **************************************************/


ComputeExactMC::ComputeExactMC(){
    m_teststat_buf = NULL;   
    m_idx=0;
}

ComputeExactMC::~ComputeExactMC(){
    
    int n = m_cohortinfo.size();
    CohortInfo * info;
    for(int i=0;i<n;i++){
        
        info = m_cohortinfo[i];
        delete info;
    }
    
    SL_free(m_teststat_buf);
    
}

int      ComputeExactMC::Init(double * Q, int nQ, int ncohort, double * aZ0, double * aZ1, int * ak, int * am, int total_m, int * amarker_idx, int * atotal, int * atotal_k, double * aprob_k, double * ap1, int * aIsExact){
    
    int i;
    int k,m, total;
    int * marker_idx, * total_k, * IsExact;
    double * Z0, * Z1, * prob_k, * p1;
    CohortInfo * info;
    m_total_m = total_m;
    m_ncohort = ncohort;
    m_total = 1;

    int idx1=0, idx2=0, idx3 = 0, idx4 = 0;
    
    for(i=0;i<nQ;i++){
        
        m_Q.push_back(Q[i]);
    }
    
    for(i=0;i<ncohort;i++){
        
        info = new CohortInfo;
        m=am[i];
        k=ak[i];
        total = atotal[i];
        
        Z0 = aZ0 + idx1;
        Z1 = aZ1 + idx1;

        total_k = atotal_k + idx3;
        prob_k = aprob_k + idx3;
        IsExact = aIsExact + idx3;
        p1 = ap1 + idx2;
        
        marker_idx = amarker_idx + idx4;
        
        idx1 += m*k;
        idx2 += k;
        idx3 += k+1;
        idx4 += m;
        
        info->Init(Z0, Z1, k, m, m_total_m, marker_idx, total, total_k, prob_k, p1, IsExact);
        
        m_cohortinfo.push_back(info);    
        m_total = m_total * total;
        
    }
    
    // Buffers for test
    m_teststat_buf = (double *) SL_calloc(m_total_m * m_ncohort, sizeof(double));
    m_fprob = (double *) SL_calloc(m_total, sizeof(double));
    m_teststat = (double *) SL_calloc(m_total, sizeof(double));
    
    return 1;
    
}

int     ComputeExactMC::Run(){

    
    double prob = 1;
    double * teststat = NULL;
    teststat = (double *) SL_calloc(m_total_m, sizeof(double));
    
    Recurse_GetTestStat(0, teststat, prob);
    
    
    SL_free(teststat);
    

    // Compute p-values 
    unsigned long i;
    double total_prob_sum = 0;
	for(i=0;i < m_total;i++){
        total_prob_sum += m_fprob[i];
	}
    
	for(i=0;i < m_total;i++){
        m_fprob[i] =  m_fprob[i]/total_prob_sum;
	}
    
    // Compute p-values 
	for(int l=0;l< m_Q.size(); l++){
		double n_num   =0;
		double n_same  =0;
        //int     n1=0;
        
		for(i=0;i<m_total ; i++){
            
            if(m_Q[l] <= m_teststat[i]){
				n_num += m_fprob[i] ;
                //Rprintf("[%e]\n", m_fprob[i]);
				if( m_teststat[i]<=m_Q[l] ){
					n_same += m_fprob[i] ;
				}
			}
		}
		
		m_pval.push_back(n_num);
		m_pval_same.push_back(n_same);
	}
    return 1;
} 





int    ComputeExactMC::Recurse_GetTestStat(int idx_cohort, double * teststat, double prob){
    
	int i, end;
    double prob_prev;
    //int  idx;
    CohortInfo * info;

	if(idx_cohort >= m_ncohort){
		put_results(teststat, prob);
        return 0;
	} else {
        info =m_cohortinfo[idx_cohort];
        end = info->GetTotal();
        //idx = idx_cohort * m_total_m;
        for(i=0; i<end; i++) {
            
            //memcpy(&(teststat_buf[idx]), teststat, m_total_m* sizeof(double));
            prob_prev = prob;
            info->Sum_TestStat(i, teststat);
            prob = prob * info->GetProb(i);
            Recurse_GetTestStat(idx_cohort+1, teststat, prob);
            
            info->Delete_TestStat(i, teststat);
            prob = prob_prev;
            
            //memcpy(teststat, &(teststat_buf[idx]), m_total_m* sizeof(double));
        }
        
	}
    return 0;
}

int  ComputeExactMC::put_results(double * teststat, double prob){
    
    int j; 
    
    m_teststat[m_idx] = 0;
	for(j=0;j< m_total_m;j++){
		m_teststat[m_idx]+=teststat[j] * teststat[j];
	}
    m_fprob[m_idx] = prob;
    m_idx++;
    
    return 1;
    
}

int     ComputeExactMC::GetPvalues(double * pval, double * pval_same){
    
    int i;
    for(i=0;i< m_pval.size(); i++){
        
        
        pval[i]= m_pval[i];
        pval_same[i] = m_pval_same[i];
        
    }
    
    return 1;
}


int     ComputeExactMC::PrintPvals(){
    
    for(int i=0;i< m_pval.size(); i++){
        //for(int i=0;i<5; i++){
        
        Rprintf("[%e][%e]\n", m_pval[i],m_pval_same[i]);        
    }
    
    return 1;
}








