#ifndef _STAND_ALONE_

    #include <R.h>
    #include <Rmath.h>
    
#endif
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "Binary_global.hpp"
#include "Binary_Permu_SKAT.h"

/*********************************************************************


************************************************************************/

Binary_Permu_SKAT::Binary_Permu_SKAT(){
       
}


Binary_Permu_SKAT::~Binary_Permu_SKAT(){

}


int     Binary_Permu_SKAT::Run_With_Dummy(int nSNP, int nSample, int nPermu){
    
    double *dZ;
    int *dY;
    
    dZ = new double[nSNP * nSample];
    dY = new int[nSample];
    
    int i,j,k;
	for(i=0;i<nSample;i++){
		if(i<nSample/2){
			dY[i]=1;
		} else {
			dY[i]=0;
		}
	}
	k=0;
	for(i=0;i<nSNP;i++){
		for(j=0;j<nSample;j++){
            
			dZ[k]=0;
			if(i==j){
				dZ[k]=1;
			}
			k++;
			
		}
	}
    

    Init(dZ, dY, nSNP, nSample, nPermu, (1.0E-6));
    Run();
    return 1;
}


int Binary_Permu_SKAT::Init(double * Z, int * Y, int nSNP, int nSample, int nPermu, double epsilon){
    
    int total = nSNP * nSample;
    
    m_nSample=nSample;
    m_nSNP = nSNP;
    m_nPermu = nPermu;
    
    m_Z.resize(total);
    m_Y.resize(nSample);
    m_buf.resize(nSample);
    m_buf1.resize(nSample);
    m_TestStat.resize(nPermu);
    
    m_epsilon = epsilon;
    m_nCase =0;
    
    
    
    int k=0;
    for(int i=0;i<nSNP;i++){
        for(int j=0;j<nSample;j++){
            
            m_Z[k] = Z[k];
            k++;
            
        }
        
        if(i==0){
            m_Y[i] = Y[i];
            m_MeanY += Y[i];
            
            if(m_Y[i] == 1){
                
                m_nCase++;
            }
        }
    }
    
    m_MeanY = m_MeanY / (double) nSample;
    
    return 1;
}

int Binary_Permu_SKAT::Run(){
   
    int i;
    Get_TestStat(0, true);
    
    for(i=0;i<m_nPermu;i++){
        
         Get_TestStat(i, false);
        
    }
    
    
    int n_num	=0;
    int n_same  =0;
    double temp1;
    
    
    //Rprintf("Start[%d]\n", l);
    //Rprintf("[%e]\n",  m_Q[l]);
    for(i=0;i<m_nPermu;i++){
        
        
        temp1 = m_OrgTestStat - m_TestStat[i];
        if(fabs(temp1) <= m_epsilon){
            temp1=  0;
        }
        //Rprintf("[%e][%e][%e]\n",  temp1, fabs(temp1), m_epsilon);
        if(temp1 <= 0){
            n_num +=1;
            if( temp1 == 0 ){
                n_same += 1;
            }
        }  
        
    }
    
    
    m_pval = (double) (n_num +1) / (double) (m_nPermu +1); 
    m_pval_same = (double) n_same / (double) m_nPermu; 
    
    //printf("%f, %f, %d, %d, %d\n", m_pval, m_pval_same, n_num, n_same, m_nPermu);
    
    return(1);
    
    
}


int Binary_Permu_SKAT::Get_TestStat(int idx, bool is_org){
    
    int i, j,k, *pY;
    double test, test1, test2, temp;
    
    if(is_org){
        pY = (int *) m_Y.data();
        
        
    } else {
        SL_GetPermu(m_nSample, m_Y, m_buf);
        pY = (int *) m_Y.data();
        
    }
    
    k=0;
    
    test=0;
    for( i=0;i<m_nSNP;i++){
        test1=0;
        test2=0;
        for( j=0;j<m_nSample;j++){
                
            if(pY[j] == 1){
                test1 += m_Z[k];
            } else {
                test2 += m_Z[k];
            }
            k++;
                
        }
        
        temp = test1 * m_MeanY + test2 * m_MeanY * (-1);
        test += temp * temp;
    }
    
    if(is_org){
        m_OrgTestStat = test;
    } else {
        m_TestStat[idx] = test;
    }
    return 1;

}

int     Binary_Permu_SKAT::GetPvalues(double * pval, double * pval_same){
    
    *pval = m_pval;
    *pval_same = m_pval_same;
    return 1;
    
}


