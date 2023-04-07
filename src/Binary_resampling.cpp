#ifndef _STAND_ALONE_

    #include <R.h>
    #include <Rmath.h>

#endif
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "Binary_global.hpp"

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

void SL_GetSample(int n, int k, std::vector<int> & y, std::vector<int> & x){
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

void SL_GetPermu(int n, std::vector<int> & y, std::vector<int> & x){
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

void SL_Binary_Boot1(int n, int ncase, std::vector<double> & pcase, std::vector<int> & buf1, std::vector<int> & buf2, std::vector<int> & z_one, int *err){

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
			memcpy(buf1.data(), buf2.data(), sizeof(int) * n1);
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

void Test1(int * pn, int * pm, int *pncase, double * pcase, std::vector<int> & buf1, std::vector<int> & buf2, int *buf3){

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
