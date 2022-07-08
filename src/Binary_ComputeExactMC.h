#ifndef COMPUTE_EXACT_MC_001 
#define COMPUTE_EXACT_MC_001 

#include <vector>
using namespace std;

class CohortInfo {
    
public: 
    CohortInfo();
    ~CohortInfo();
    int     Init(double * Z0, double *Z1, int k, int m, int total_m, int * marker_idx, int total, int * total_k, double *prob_k, double * p1, int * IsExact);
    
    int     Exact_Recurse(int k, int * array, int cell, int start, int end, int is_case);
    int     Sum_TestStat(int idx, double * teststat);
    int     Delete_TestStat(int idx, double* teststat);
    double  GetProb(int idx);
    
    int     CalTestStat(int k, int * array, int is_case);
    int     CalFisherProb(int k, int * array, int is_case);
    
    int     GetTotal(){ return m_total; } ;
    

public:
    
    double * m_Z0;		/* Z0 and Z1 matrix */
    double * m_Z1;
    
    double * 	m_teststat_Z0;
    double *    m_teststat_Z1;
    double *    m_teststat_one;
    
    double  * m_teststat_all;
    double * m_fprob ;	/* fisher hg probability */
    // buffers 
    int * m_array;
    
    
    int m_total_m;			/* # of total markers over all cohorts  */
    int m_m;                /* # of polymorphic markers in the cohort */
    
    int m_k;				/* # of samples to have at least one minor allele */
    int m_total;
    
    vector<int> m_marker_idx;   /* marker idx */
    vector<int> m_total_k;      /* # of all combinations that #of case = k*/
    vector<double> m_prob_k;    /* probability (# of case =k)*/
    vector<double> m_p1;		/* probability to be case */
    vector<double> m_odds;		/* probability to be case */
    vector<double> m_denomi;
    vector<int>     m_IsExact;
    

    vector<int>     m_combi_case;    /* 1 = combination for case, 0 = combination for control */
    int     m_idx;


    double m_pprod;
    
    
};


// multi cohort version 
class ComputeExactMC {
 
public: 
    ComputeExactMC();
    ~ComputeExactMC();

    int     PrintPvals()  ;
    int     Init(double * Q, int nQ, int ncohort, double * aZ0, double * aZ1, int * ak, int * am, int total_m, int * amarker_idx, int * atotal, int * atotal_k, double * aprob_k, double * ap1, int * aIsExact);
    int     Recurse_GetTestStat(int idx_cohort, double * teststat, double prob);
    int     Run();
    
    int     put_results(double * teststat, double prob);
    int     GetPvalues(double * pval, double * pval_same);


protected:
    vector<CohortInfo *>    m_cohortinfo;
    int m_ncohort;
    int m_total_m;
    vector<double> m_Q;
    
    int m_idx;
    double * m_teststat_buf;
    
    double * m_fprob ;	/* fisher hg probability */
    double * m_teststat;	/* test statistic */
    
    unsigned long m_total;
    
    // OUTput
    vector<double>  m_pval;
    vector<double>  m_pval_same;
    
};

#endif 

