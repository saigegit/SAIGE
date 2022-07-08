#ifndef _STAND_ALONE_

    #include <R.h>
    #include <Rmath.h>
    
#endif
#include <math.h>
#include <stdlib.h>
#include <string.h>

/*********************************************************************


************************************************************************/
double SL_runif_double(){

    double val ;
    #ifdef _STAND_ALONE_
        val=((double)rand()/(double)RAND_MAX);
    #else
        val = unif_rand();
    #endif
    return val ;
 
}

int SL_runif_INT(int max){
    
    int val ;
#ifdef _STAND_ALONE_
    val=rand() % max;
#else
    val = floor(unif_rand() * (double) max);
#endif
    return val ;
    
}

int SL_runif_int(int max){

	int re;
	re = SL_runif_INT(max);
	return re;
}


void  SL_setseed(int seed){
  
    #ifdef _STAND_ALONE_
        srand ( seed );
    #else
        GetRNGstate();
    #endif
}

void  SL_out(){
    
    #ifndef _STAND_ALONE_
        PutRNGstate();
    #endif
}



/*
	n : number of samples
	k : number of samples to be sampled
	y : output
	x : buffer
*/

void SL_GetSample(int n, int k, int *y, int *x){
    int i, j;
    
    for (i = 0; i < n; i++){
        x[i] = i;
    }
    for (i = 0; i < k; i++) {
        j = SL_runif_INT(n);
        y[i] = x[j] ;
        x[j] = x[--n];
    }
}

void SL_GetPermu(int n, int *y, int *x){
    int i, j;
    
    for (i = 0; i < n; i++){
        x[i] = y[i];
    }
    for (i = 0; i < n; i++) {
        
        j = SL_runif_INT(n);
        y[i] = x[j] ;
        x[j] = x[--n];
    }
}


/*
    n : sample size
    ncase : number of cases
    pcase : prob to be case
    buf1 : buffer of length n
    buf2 : buffer of length n 
    z_one : output of length n
    err: error code
 
 */

void SL_Binary_Boot1(int n, int ncase, double * pcase, int * buf1, int * buf2, int * z_one, int *err){

	int i, i1,j,k, n1, ncase1;
	double temp;
    
	SL_GetSample(n, n, buf1, buf2);
	ncase1 = 0;
	n1 = n;
	for(k=0;k<500;k++){
		i1=0;
		for(i=0;i<n1;i++){
			j=  buf1[i];
			temp = SL_runif_double();
			if(temp <= pcase[j]){
				z_one[j]=1;
				ncase1++;
			} else {
				buf2[i1] = j;
				i1++;
			}
			if(ncase1 == ncase)
				break;
		} 
		
		if(ncase1 == ncase){
			break;
		} else if(ncase1 > ncase) {
			*err = -1;
			return;
		} else {
			n1 = n - ncase1;
			memcpy(buf1, buf2, sizeof(int) * n1);		
		}
	}
	
	if(ncase != ncase1){
		*err = -1;
        /*
        printf("Error! [%d][%d][%d]\n", n, ncase, ncase1);
        for(i=0;i<n;i++){
            
            printf("[%e]",pcase[i]);
        }
        printf("\n");
        */
        
		return;
	}
	
	*err=1;
	return;

}



void SL_Binary_Boot_2(int * pn, int * pm, int *pncase, double * pcase, int * buf1, int * buf2, int *Z, int *err ){

	int i, n, m, ncase, idx;
	int * z_one;
	n = *pn;
	m = *pm;
	ncase = *pncase;
	
    SL_setseed(100);

	
	for(i=0;i<m;i++){
	
		idx = i * n ;
		z_one = &(Z[idx]);
		SL_Binary_Boot1(n, ncase, pcase, buf1, buf2, z_one, err);
		if(*err == -1){
			SL_out();
			return;
		}
	}
	
	SL_out();
	return;
}


void Test1(int * pn, int * pm, int *pncase, double * pcase, int * buf1, int * buf2, int *buf3){

	int i, n, m, ncase, idx;
	n = *pn;
	m = *pm;
	ncase = *pncase;
	
	for(i=0;i<m;i++){
	
		idx = i * n ;
		/*SL_Binary_Boot1(n, ncase, pcase, buf1, buf2, buf3);*/
		SL_GetSample(n, n, buf1, buf2);
	}
}



int  CalTestStat(double * Z0, double *Z1, double * teststat_Z0, 
                     double * teststat_one, int m, int n, int * array, 
                     double * pQ, int is_inverse){
    
    
	int i, j, temp, arr1;
    double test_skat = 0;

    arr1=1;
    if(is_inverse > 0){
        arr1=0;
    }
    
	
    memcpy(teststat_one, teststat_Z0, sizeof(double) * m); 
     
    for(i=0;i< n;i++){
        
        if(array[i]== arr1){
            temp = i*m;
            for(j=0;j< m;j++){
                /*printf("[%d][%d][%d][%e][%e]\n",i, array[i], j, teststat_one[j], Z1[temp+j] - Z0[temp+j]);*/
                teststat_one[j]+=Z1[temp+j] - Z0[temp+j] ;
                
            }
        }
    }
   
	for(j=0;j< m;j++){
		test_skat+=teststat_one[j] * teststat_one[j];
	}
    
    pQ[0]= test_skat; 
    
    /*printf("|%e|\n",test_skat);*/
    return 1;
}


// int     ComputeExactSKATO::CalTestStat(int k, int * array){
int  CalTestStat_O(double * Z0, double *Z1, double * Z0_C, double * Z1_C, 
                     double * teststat_Z0, double teststat_Z0_C, 
                     double * teststat_one, int m, int n, int * array, 
                     double * r_corr, int n_r, double * pQ, int is_inverse){
    
    
	int i, j, temp, arr1;
    double test_C  = teststat_Z0_C;
    double test_skat = 0;
   
    
    if(n_r == 1){
        int re;
        re = CalTestStat(Z0, Z1, teststat_Z0, teststat_one, m, n, array, pQ, is_inverse);
        return re;
    }
    
    memcpy(teststat_one, teststat_Z0, sizeof(double) * m);
    
    arr1=1;
    if(is_inverse > 0){
        arr1=0;
    }
	
    for(i=0;i< n;i++){
        
        if(array[i]== arr1){
            temp = i*m;
            for(j=0;j< m;j++){
                teststat_one[j]+=Z1[temp+j] - Z0[temp+j] ;
            }
            test_C += Z1_C[i] - Z0_C[i] ;
        }
    }
    
	for(j=0;j< m;j++){
		test_skat+=teststat_one[j] * teststat_one[j];
	}
    test_C = test_C * test_C;
    
    for(i=0;i<n_r;i++){
        pQ[i] = (1-r_corr[i]) * test_skat + r_corr[i] * test_C;
        
    }
    
    return n_r;
}

void ResampleSTAT_1(double * Z0, double *Z1, double * Z0_C, double * Z1_C, 
                  double * teststat_Z0, double *teststat_Z1, double *pteststat_Z0_C, double *pteststat_Z1_C,
                  double * r_corr, int *pn_r, int *pk, int *pm, int * pn,
                  int * total_k, int * ncase_k, double * p1,
                  int *buf1, int * buf2, int *buf3, double * teststat_one, /* buffers */
                  double * Q, int *err)
{
  
    double teststat_Z0_C, teststat_Z1_C;
    int n_r, k, m, n;
    int i,j, idx, idx1;
    double *pQ = Q;
    
    teststat_Z0_C = pteststat_Z0_C[0]; 
    teststat_Z1_C = pteststat_Z1_C[0]; 
    n_r = pn_r[0]; k=pk[0]; m=pm[0]; n=pn[0];
    
    *err=1;
    SL_setseed(100);
    
    /* printf("Start[%d]!\n",k);*/
    idx=0;
    for(i=0;i<=k;i++){
        
        int is_sampling=1;
        int is_inverse =0;
        
        double * Z0_1=Z0;
        double * Z1_1=Z1;
        double * Z0_C_1=Z0_C;
        double * Z1_C_1=Z1_C;
        double * teststat_Z0_1=teststat_Z0;
        double teststat_Z0_C_1 = teststat_Z0_C;
        
        
        /*printf("[%d][%d]!\n",i, total_k[i]);*/
        if(ncase_k[i] == 0){
            for(j=0;j< n;j++){
                buf3[j]= 0;
            }
            is_sampling=0;
                
        } else if(ncase_k[i] == n){
            for(j=0;j< n;j++){
                buf3[j]= 1;
            }
            is_sampling=0;
            is_inverse=1;
                
        } else { 
            if(ncase_k[i] *2 > n) {
                is_inverse=1;
            }
            double sum1=0;
            for(j=0;j<n;j++){
                sum1 += p1[j];
            }
            for(j=0;j<n;j++){
                p1[j] = p1[j] /sum1 * ((double) ncase_k[i]);
            }           
        } 
        
        if(is_inverse > 0){
            
            Z0_1=Z1;
            Z1_1=Z0;
            Z0_C_1=Z1_C;
            Z1_C_1=Z0_C;
            teststat_Z0_1=teststat_Z1;
            teststat_Z0_C_1 = teststat_Z1_C;
            
        }
        
        for(j=0;j< total_k[i]; j++){
            
            if(is_sampling){
                memset(buf3, 0, sizeof(int) * n);
                SL_Binary_Boot1(n, ncase_k[i], p1, buf1, buf2, buf3, err);
                if(*err == -1){
                    SL_out();
                    return;
                }
            }
            
            idx1 = CalTestStat_O(Z0_1, Z1_1, Z0_C_1, Z1_C_1, teststat_Z0_1, teststat_Z0_C_1
                                 , teststat_one, m, n, buf3, r_corr, n_r, &(pQ[idx]), is_inverse);
            
            idx+=idx1;
             
        }
         

    }
    
    SL_out();
    return;
    
    
}

