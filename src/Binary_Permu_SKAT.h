#ifndef Binary_Permu_SKAT_00 
#define Binary_Permu_SKAT_00



#include <vector>



using namespace std;


class Binary_Permu_SKAT {
    
public: 
    Binary_Permu_SKAT();
    ~Binary_Permu_SKAT();
    
    int     Init(double * Z, int * Y, int nSNP, int nSample, int nPermu, double epsilon);
    int     Run();    
    int     Get_TestStat(int idx, bool is_org=true);
    int     GetPvalues(double * pval, double * pval_same);
    int     Run_With_Dummy(int nSNP, int nSample, int nPermu);
    
protected:

    int m_nSNP;
    int m_nSample;
    int m_nPermu;
    int m_nCase;
  
    double m_nMin;
    
    vector<double>  m_Z;
    vector<int>     m_Y;
    vector<int>     m_buf;
    vector<int>     m_buf1;
    vector<double>  m_TestStat;
    
    double      m_OrgTestStat;
    double      m_MeanY;
    
    double  m_pval;
    double  m_pval_same;

    
    // precision
    double m_epsilon;
    
    
    
};


#endif 

