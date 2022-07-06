#ifndef COMPUTE_EXACT_001 
#define COMPUTE_EXACT_001 

#include <vector>



using namespace std;


class ComputeExact {
 
public: 
    ComputeExact();
    ~ComputeExact();
    
    int     Init(int * resarray, int nres, int * nres_k, double * Z0, double *Z1, int k, int m, int total
             , int * total_k, double *prob_k, double * odds, double * p1, int * IsExact, double epsilon, bool IsSmallmemory=0);

    
    int     Run(int test_type);
    int     GetPvalues(double * pval, double * pval_same, double * prob_k, double * minP = NULL);
    int     PrintPval();
protected:
    
    int     SaveParam(double * Z0, double *Z1, int k, int m, int total, int * total_k, double *prob_k, double * odds, double * p1, int * IsExact, double epsilon, bool IsSmallmemory);
    virtual double     CalTestStat(int k, int * array, bool is_save=true, bool is_minIdx = false, int * minIdx = NULL);
    virtual double     CalTestStat_INV(int k, int * array, bool is_save=true, bool is_minIdx = false, int * minIdx = NULL);
    
    int     CalFisherProb(int k, int * array);
    int     SKAT_Exact_Recurse(int k, int * array, int cell, int start, int end);
    int     SKAT_Resampling(int k, int * array);
    int     SKAT_Resampling_Random(int k, int * array);    

    int     CalFisherProb_INV(int k, int * array);
    int     SKAT_Exact_Recurse_INV(int k, int * array, int cell, int start, int end);
    
protected:

    double * m_fprob ;	/* fisher hg probability */
    double * m_teststat;	/* test statistic */
    double * m_Z0;		/* Z0 and Z1 matrix */
    double * m_Z1;

    double * 	m_teststat_one;
    double * 	m_teststat_Z0;
    double *    m_teststat_Z1;

    int m_k;				/* # of samples to have at least one minor allele */
    int m_m;				/* # of markers */
    int m_total;
    double m_pprod;         /* product of prob */
    vector<int> m_total_k;
    vector<double> m_prob_k;
    vector<double> m_odds;    
    vector<double> m_p1;		/* probability to be case */
    vector<double> m_p1_inv;
    vector<double> m_Q;
    vector<double> m_denomi;
    vector<int>     m_IsExact;
    
    vector<double> m_logOdds;
    
    int m_idx;
    int * m_temp_x;
    int * m_temp_x1;
    
    // OUTput
    vector<double>  m_pval;
    vector<double>  m_pval_same;
    
    bool m_IsSmallmemory ;
    
    
    // Investigate possible min p-value
    double m_LargestQ;
    double m_minP;
    
    // for debug
    vector<double> m_pr1_debug;
    
    // precision
    double m_epsilon;
    

};

class ComputeExactSKATO : public ComputeExact {

public: 
    ComputeExactSKATO() ;
    ~ComputeExactSKATO(){};
    
    int     Init(int * resarray, int nres, int * nres_k, double * Z0, double *Z1, double * r_corr, int n_r, double * param, int k, int m, int total, int * total_k, double *prob_k, double * odds, double * p1, int * IsExact, double epsilon, bool IsSmallmemory=0);

protected:
    
    virtual double     CalTestStat(int k, int * array, bool is_save=true, bool is_minIdx = false, int * minIdx = NULL);
    virtual double     CalTestStat_INV(int k, int * array, bool is_save=true, bool is_minIdx = false, int * minIdx = NULL);
    
    double  Cal_Pvalue_Rcorr(double test_skat, double test_C, int idx);

protected:
    /* for SKAT-O */
    vector<double>  m_rcorr;
    double * m_Z0_C;		
    double * m_Z1_C;    
    double 	m_teststat_Z0_C;
    double  m_teststat_Z1_C;
    
    vector<double> m_muQ;
    vector<double> m_varQ;
    vector<double> m_df;
    
    vector<int> m_minIdx;

    
};

#endif 

