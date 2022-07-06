#ifndef HyperGeo_001 
#define HyperGeo_001 



#include <vector>
#include <cstdio>         // std::remove


using namespace std;


class HyperGeo {
    
public: 
    HyperGeo(){};
    ~HyperGeo();
    int     HyperGeo::Run(int k, int ngroup, int ncase, std::vector<int> & group,  std::vector<double> & weight);
    int     HyperGeo::Get_lprob(std::vector<double> & prob); 
    int     Print();
    double  lCombinations(int n, int k);
    
protected:
    double  GetLogProb(int idx, int i);
    int     SaveProb(double lprob, int ncase_used);
    int     Recursive(double lprob, int idx, int ncase_used);

    int m_ngroup;
    int m_ncase;
    vector<int> m_group;   
    vector<double> m_lweight; 
    vector<double> m_kprob;
    int m_k;
    
    vector<double * > m_probtbl;
    
    double  m_ref;
   
    
};


#endif 

