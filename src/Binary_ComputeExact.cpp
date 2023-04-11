#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <float.h>

#include "Binary_ComputeExact.hpp"
#include "Binary_global.hpp"

#ifdef _STAND_ALONE_
#define Rprintf printf
#else
#include <R.h>
#include <Rmath.h>
#endif

#define MIN_SIM1    0
//#define EPSILON    (1.0E-7)
//#define ERR1 0
/***************************************************

   Protected functions
 
***************************************************/

double     ComputeExact::CalTestStat(int k, int * array, bool is_save,bool is_minIdx , int * minIdx){

	int i, j, l, temp;
    double stat = 0;
    m_teststat_one = m_teststat_Z0;

    for(i=0;i< k;i++){
        l = array[i];
        temp = l*m_m;
        for(j=0;j< m_m;j++){
            m_teststat_one.at(j) += m_Z1.at(temp +j) - m_Z0.at(temp + j);
        }
    }

	for(j=0;j< m_m;j++){
        stat +=m_teststat_one.at(j)* m_teststat_one.at(j);
	}
    
    //stat=log(stat);
    if(is_save){
        m_teststat.at(m_idx) = stat ;
    }
    return stat;

}


double     ComputeExact::CalTestStat_INV(int k, int * array, bool is_save, bool is_minIdx, int * minIdx ){
    
	int i, j, l, temp;
    double stat = 0;
    m_teststat_one = m_teststat_Z1;
	
    for(i=0;i< k;i++){
        l = array[i];
        temp = l*m_m;
        for(j=0;j< m_m;j++){
            m_teststat_one.at(j) += m_Z0.at(temp+j) - m_Z1.at(temp+j);
        }
    }
    
	for(j=0;j< m_m;j++){
		stat+=m_teststat_one.at(j) * m_teststat_one.at(j);
	}
    
    //stat=log(stat);
    if(is_save){
        m_teststat.at(m_idx)=stat ;
    }
    return stat;
    
}


int     ComputeExact::CalFisherProb(int k, vector<int> & array){

	int i,l;
	double temp = 1;
	for(i=0;i< k;i++){
		l = array.at(i);
		temp = temp * m_odds.at(l);
	}
	m_fprob.at(m_idx) = temp;
    m_denomi.at(k) = m_denomi.at(k) + temp;
	
    return 0;

}


int     ComputeExact::CalFisherProb_INV(int k, vector<int> & array){
    
	int i,l, k1;
    k1 = m_k - k;
	double temp = m_pprod;
	for(i=0;i< k;i++){
        l = array.at(i);
		temp = temp / m_odds.at(l);
	}
	m_fprob.at(m_idx) = temp;
    m_denomi.at(k1) = m_denomi.at(k1)+temp;
	
    return 0;
    
}


int     ComputeExact::SKAT_Exact_Recurse(int k, vector<int> & array, int cell, int start, int end){

	int i;
	
	/* node */
	if(k == cell){
		CalTestStat(k, array.data());
		CalFisherProb(k, array);
		m_idx++;
	} else {
		 for(i=start; i<end; i++) {
            array.at(cell) = i;
            SKAT_Exact_Recurse(k, array, cell+1, i+1, end);
        }
	
	}
    return 0;
}


int     ComputeExact::SKAT_Exact_Recurse_INV(int k, vector<int> & array, int cell, int start, int end){
    
	int i;

	/* node */
	if(k == cell){
		CalTestStat_INV(k, array.data());
		CalFisherProb_INV(k, array);
		m_idx++;
	} else {
        for(i=start; i<end; i++) {
            array.at(cell) = i;
            SKAT_Exact_Recurse_INV(k, array, cell+1, i+1, end);
        }
        
	}
    return 0;
}


int     ComputeExact::SKAT_Resampling(int k, vector<int> & array){
    
    int k1 = m_k-k;
    if(k <= m_k/2 +1){
        for(int i=0;i< m_total_k.at(k);i++){
            SL_Sample(k, m_k, m_temp_x, array);
            CalTestStat(k, m_temp_x.data());
            CalFisherProb(k, m_temp_x);
            m_idx++;
        }
    } else {
        
        for(int i=0;i< m_total_k.at(k);i++){
            SL_Sample(k1, m_k, m_temp_x, array);
            CalTestStat_INV(k1, m_temp_x.data());
            CalFisherProb_INV(k1, m_temp_x);
            m_idx++;
        }
    }
    
    
    return 1;
}


int     ComputeExact::SKAT_Resampling_Random(int k, std::vector<int> & array){
    
    int err;
    int k1 = m_k-k;
    if(k <= m_k/2 +1){
        for(int i=0;i< m_total_k.at(k);i++){
            SL_Binary_Boot1(m_k, k, m_p1, array, m_temp_x1, m_temp_x, &err);
            CalFisherProb(k, m_temp_x);

            m_fprob.at(m_idx) = 1;
            m_denomi.at(k) = m_denomi.at(k)+1;

            m_idx++;
        }
    } else {

        for(int i=0;i< m_total_k.at(k);i++){

            SL_Binary_Boot1(m_k, k1, m_p1_inv, array, m_temp_x1, m_temp_x, &err);
            CalFisherProb_INV(k1, m_temp_x);

            m_fprob.at(m_idx) = 1;
            m_denomi.at(k) = m_denomi.at(k)+1;
            
            m_idx++;
        }
    }
    
    
    return 1;
}


/***************************************************
 
    Public function
 
 **************************************************/


ComputeExact::ComputeExact(){

    m_pprod=1;

}

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
*   odds : odds to be case
*	p1 : prob to be case
*	pval : outcome p-value
*
******************************************************/

int     ComputeExact::Init(vector<int> resarray, int nres, int * nres_k, double * Z0, double *Z1, int k, int m, int total, int * total_k, double *prob_k, double * odds, double * p1, int * IsExact, double epsilon, bool IsSmallmemory){


    int i, idx, k1;
    int * array;
    double Q;

    SaveParam(Z0, Z1, k, m, total, total_k, prob_k, odds, p1, IsExact, epsilon, IsSmallmemory);
    
    idx = 0;
    for(i=0;i<nres;i++){

        // Can this pointer go stale? Does resarray get appended to / reallocated past here?
        array = (resarray.data() + idx) ;
        k1 = nres_k[i];
        idx += k1;
        
        Q = CalTestStat(k1, array, false);
        //Rprintf("T-[%e]\n", Q);
        m_Q.push_back(Q);
            
    }


    return 1;
}



int     ComputeExact::SaveParam(double * Z0, double *Z1, int k, int m, int total, int * total_k, double *prob_k, double * odds, double * p1, int * IsExact, double epsilon, bool IsSmallmemory){
    
	int i,j;
	m_idx=0;
	
	m_k = k;
	m_m = m;
	m_total = total;
    m_IsSmallmemory = IsSmallmemory;
    m_epsilon=epsilon;
    
    // Init _k vectors
    m_pprod=1;
    for(i=0;i<=k;i++){
        
        m_total_k.push_back(total_k[i]);
        m_prob_k.push_back(prob_k[i]);
        m_denomi.push_back(0);
        m_IsExact.push_back(IsExact[i]);
        
        if(i<k){
            m_p1.push_back(p1[i]);
            m_odds.push_back(odds[i]);
            m_pprod = m_pprod*odds[i];
            
            m_p1_inv.push_back(1 - p1[i]);
        }
    }
 

	
	m_Z0.resize(m_k * m_m);
	m_Z1.resize(m_k * m_m);
	m_teststat_Z0.resize(m_m);
    m_teststat_Z1.resize(m_m);

    memcpy(m_Z0.data(), Z0, sizeof(double) * m_k * m_m);
    memcpy(m_Z1.data(), Z1, sizeof(double) * m_k * m_m);

	/* generate prob matrix */
	for(i=0;i< m_k;i++){
		int idx = i*m_m;
		for(j=0;j< m_m;j++){
			m_teststat_Z0.at(j)+=m_Z0.at(idx+j);
            m_teststat_Z1.at(j)+=m_Z1.at(idx+j);
		}
        
        // debug
        m_pr1_debug.push_back(0);
	}
	
    if(!m_IsSmallmemory){
        m_fprob.resize(m_total);
        m_teststat.resize(m_total);
    } else {
        
        m_fprob.clear();
        m_teststat.clear();
        
    }
    
	m_teststat_one.resize(m_m);
    m_temp_x.resize(m_k);
    m_temp_x1.resize(m_k);
    
    return 1;
}



/*
 *  test_type == 1 : original 
    test_type == 2 : new
    test_type == 3 : use SKAT_Resampling_Random when count < 5000
 */

int     ComputeExact::Run(int test_type){
   
    int i, j, idx, l;
    vector<int> array(m_k);
    //SL_setseed(time(NULL));
	SL_setseed(1);
  
	for(i=0;i < m_k+1;i++){
        
        if(m_IsExact.at(i) == 1){
            if(i <= m_k/2 +1){
                SKAT_Exact_Recurse(i, array, 0, 0, m_k);
            } else {
                SKAT_Exact_Recurse_INV(m_k - i, array, 0, 0, m_k);
            }
        } else if( m_total_k.at(i) < MIN_SIM1 && test_type == 3) {
            
            SKAT_Resampling_Random(i, array);

        } else {
            

            SKAT_Resampling(i, array);
            
        }
	}

	//Rprintf("2, m_idx[%d], total[%d]\n", m_idx, m_total);
	
   
    //Get probability 
    idx=0;
    double total_prob_sum = 0;
	for(i=0;i < m_k+1;i++){

		for(j=idx;j< idx + m_total_k.at(i);j++){
            
            m_fprob.at(j) = m_fprob.at(j) / m_denomi.at(i) * m_prob_k.at(i);
            total_prob_sum += m_fprob.at(j);
		}	
		idx = idx + m_total_k.at(i);
	}
    idx=0;
	for(i=0;i < m_k+1;i++){
        m_prob_k.at(i) = 0;
		for(j=idx;j< idx + m_total_k.at(i);j++){
            m_fprob.at(j) = m_fprob.at(j) / total_prob_sum;
            m_prob_k.at(i) += m_fprob.at(j); // for debugging
		}	
        //Rprintf(".at(%d:%e)", i, m_prob_k.at(i));
		idx = idx + m_total_k.at(i);
	}

    
	//Rprintf("3\n");
    double temp1=0;
    for(l=0;l< m_Q.size(); l++){
		double n_num	=0;
		double n_same  =0;
       
        
        //Rprintf("Start[%d]\n", l);
        //Rprintf("[%e]\n",  m_Q[l]);
		for(i=0;i<m_total ; i++){
            
                        
            temp1 = m_Q.at(l) - m_teststat.at(i);
            //Rprintf("[%e][%e][%e][%e]\n", m_Q[l],m_teststat[i], temp1);
            
            if(fabs(temp1) <= m_epsilon){
                temp1=  0;
            }
            //Rprintf("[%e][%e][%e]\n",  temp1, fabs(temp1), m_epsilon);
            if(temp1 <= 0){
				n_num += m_fprob.at(i) ;
				//Rprintf("C\n");
				if( temp1 == 0 ){
					n_same += m_fprob.at(i) ;
					//Rprintf("D\n");
				}
			}  
            
		}

        
		m_pval.push_back(n_num);
		m_pval_same.push_back(n_same);
		//Rprintf("[%e][%e][%d]\n", m_pval[l], m_pval_same[l], l);
	}

    m_LargestQ=m_teststat.at(0);
    m_minP= 0;
    for(i=0;i<m_total ; i++){
        
        temp1 = m_LargestQ - m_teststat.at(i);
        if(fabs(temp1) <= m_epsilon){
            temp1=  0;
        }

        if(temp1 < 0){
            m_LargestQ = m_teststat.at(i);
            m_minP = m_fprob.at(i) ;
        } else if( temp1 == 0 ){
            m_minP += m_fprob.at(i) ;
           
        }
    }
    //Rprintf("LargeQ[%e][%e]\n", m_LargestQ, m_minP);
	return 1;
	
} 

int     ComputeExact::GetPvalues(double * pval, double * pval_same, double * prob_k, double * minP){
    
    int i;
    for(i=0;i< m_pval.size(); i++){
    
        
        pval[i]= m_pval[i];
        pval_same[i] = m_pval_same[i];
        
    }


	for(i=0;i < m_k+1;i++){
        prob_k[i] = m_prob_k[i];
	}
    
    if(minP != NULL){
        *minP = m_minP;
    }
    
    return 1;
}


int     ComputeExact::PrintPval(){
    
    for(int i=0;i< m_pval.size(); i++){
    //for(int i=0;i<5; i++){
        
       Rprintf("[%e][%e]\n", m_pval[i],m_pval_same[i]);        
    }
    
    Rprintf("MinP: [%e]\n", m_minP); 
    return 1;
}

