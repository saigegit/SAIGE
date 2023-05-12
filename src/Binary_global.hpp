#ifndef GLOBAL_001 
#define GLOBAL_001 

#include <vector>

void * SL_calloc(size_t num, size_t size);
void SL_free(void * ptr);
void SL_Sample(int k, int n, std::vector<int> & y, std::vector<int> & x);


void SL_Binary_Boot1(int n, int ncase, std::vector<double> & pcase, std::vector<int> & buf1, std::vector<int> & buf2, std::vector<int> & z_one, int *err);
double SL_runif_double();
int SL_runif_INT(int max);
void  SL_setseed(int seed);
void  SL_out();
void SL_GetSample(int n, int k, std::vector<int> & y, std::vector<int> & x);
void SL_GetPermu(int n, std::vector<int> & y, std::vector<int> & x);

void SKAT_Exact(std::vector<int> & resarray, int nres, int * nres_k, double * Z0, double *Z1, int k, int m, int total, int * total_k, double *prob_k, double * odds, double * p1, int * IsExact, double * pval, double *pval_same, double *minP, int test_type, double epsilon);

void SKATO_Exact(int * resarray, int nres, int * nres_k, double * Z0, double *Z1, double * r_corr, int n_r,
                 double * param, int k, int m, int total, int * total_k, double *prob_k, double * odds, double * p1, int * IsExact, double * pval, double *pval_same, double *minP, int test_type, double epsilon);

#endif 

