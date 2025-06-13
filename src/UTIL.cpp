
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "UTIL.hpp"
#include <sys/time.h>

arma::vec nb(unsigned int n){
  return(Rcpp::rbinom(n,1,0.5));
}



double getWeights(std::string t_kernel, 
                  double t_freq, 
                  arma::vec t_wBeta)
{
  if(t_wBeta.size() != 2)
    Rcpp::stop("The size of argument t_wBeta should be 2.");
  
  double weights;
  if(t_kernel == "linear")
    weights = 1;
  
  if(t_kernel == "linear.weighted"){
    Rcpp::NumericVector freq = {t_freq};
    Rcpp::NumericVector temp = Rcpp::dbeta(freq, t_wBeta(0), t_wBeta(1));
    weights = temp(0);
  }
  
  return weights;
}

void imputeGeno(arma::vec& t_GVec, 
                double t_altFreq, 
                std::vector<uint32_t> t_indexForMissing,
                std::string t_imputeMethod) 
{
  int nMissing = t_indexForMissing.size();
  
  double imputeG = 0;
  
  if(t_imputeMethod == "mean")
    imputeG = 2 * t_altFreq;
  
  if(t_imputeMethod == "none")
    imputeG = arma::datum::nan;
  
  if(t_imputeMethod == "bestguess")
    imputeG = std::round(2 * t_altFreq);
  
  for(int i = 0; i < nMissing; i++){
    uint32_t index = t_indexForMissing.at(i);
    t_GVec.at(index) = imputeG;
  }
}

// used in Main.cpp::mainMarkerInCPP
bool imputeGenoAndFlip(arma::vec& t_GVec, 
                       double & t_altFreq,
		       double & t_altCount, 
                       std::vector<uint32_t> & t_indexForMissing,
                       std::string t_impute_method,
		       double t_dosage_zerod_cutoff,
		       double t_dosage_zerod_MAC_cutoff, 
		       double & t_MAC,
		       std::vector<uint> & t_indexZero,
		       std::vector<uint> & t_indexNonZero)   
{
  bool flip = false;
  t_indexNonZero.clear();
  t_indexZero.clear();
  int nMissing = t_indexForMissing.size();
  uint dosagesSize = t_GVec.size();
  double imputeG = 0;
  if(t_altFreq > 0.5){
    flip = true;
    t_GVec = 2 - t_GVec;
    t_altFreq = 1 - t_altFreq; 
  }

if(nMissing > 0){

switch(string_to_case.at(t_impute_method)) {
  case 1:
    imputeG = std::round(2 * t_altFreq);
    //std::cout << "t_impute_method " << t_impute_method << std::endl;
    break;
  case 2:
    imputeG = 2 * t_altFreq;
    //std::cout << "t_impute_method " << t_impute_method << std::endl;
    break;
  case 3:
    imputeG = 0;
    //std::cout << "t_impute_method " << t_impute_method << std::endl;
    break;
}


  for(int i = 0; i < nMissing; i++){
    uint32_t j = t_indexForMissing.at(i);
    t_GVec.at(j) = imputeG;
  }

 t_MAC = t_MAC + imputeG * nMissing;

}

  
  if(t_dosage_zerod_cutoff > 0){ 
    if(t_MAC <= t_dosage_zerod_MAC_cutoff){
      t_GVec.clean(t_dosage_zerod_cutoff);
    }	  
  } 

  //if(nMissing > 0 || t_dosage_zerod_cutoff > 0){
     t_altCount = arma::sum(t_GVec);
     t_altFreq = t_altCount / (2*dosagesSize);
     if(flip){
	t_altFreq = 1 - t_altFreq;
        t_altCount = 2*dosagesSize - t_altCount;	
	//t_altCount = 2 * t_altFreq * dosagesSize;
     }	     
  //}

 for(unsigned int i = 0; i < dosagesSize; i++){
 	if(t_GVec(i) == 0){   
 		t_indexZero.push_back(i);
	}else{
		t_indexNonZero.push_back(i);
	}	
 }

  return flip;
}


double getInnerProd(arma::mat& x1Mat, arma::mat& x2Mat)
{
  double innerProd = arma::accu(x1Mat % x2Mat);
  return(innerProd);
}

arma::vec getTime(){
  arma::vec Time(2, arma::fill::zeros);
  struct timeval time;
  Time(0) = 0;
  if(!gettimeofday(&time,NULL))
    Time(0) = (double)time.tv_sec + (double)time.tv_usec * .000001;
  Time(1) = (double)clock() / CLOCKS_PER_SEC;
  return Time;
}

void printTime(arma::vec t1, arma::vec t2, std::string message){
  double wallTime = t2(0) - t1(0);
  double cpuTime = t2(1) - t1(1);
  if(wallTime < 60){
    Rprintf ("It took %f seconds (%f CPU seconds) to %s.\n",
             wallTime, cpuTime, message.c_str());
  }else if(wallTime < 3600){
    Rprintf ("It took %f minutes (%f CPU minutes) to %s.\n",
             wallTime/60, cpuTime/60, message.c_str());
  }else{
    Rprintf ("It took %f hours (%f CPU hours) to %s.\n",
             wallTime/3600, cpuTime/3600, message.c_str());
  }
}

double getinvStd(double t_freq)
{
  double Std = sqrt(2 * t_freq * (1-t_freq));
  if(Std == 0)
    return 0;
  else
    return 1/Std;
}

// [[Rcpp::export]]
double sum_arma1(arma::vec& X) {
    double sum = 0;
    for (uint i = 0; i < X.n_elem; ++i) {
        if (arma::is_finite(X(i)))
            sum += X(i);
    }
    return sum;
}


// [[Rcpp::export]]
double add_logp(double p1, double p2)
{
        using namespace Rcpp;
        p1 = -std::abs(p1);
        p2 = -std::abs(p2);
        double maxp = std::max(p1,p2);
        double  minp = std::min(p1,p2);
        double result = maxp+std::log(1+std::exp(minp-maxp));
        return(result);
}


bool imputeGenoAndFlip_sub(arma::vec& t_GVec,
                       arma::vec& t_GVec_sub, // New parameter to store the subset
                       double & t_altFreq,
                       double & t_altCount,
		       arma::uvec & t_indexForMissing_arma_sub,
		       arma::uvec & t_indexForMissing_arma,
                       std::string t_impute_method,
                       double t_dosage_zerod_cutoff,
                       double t_dosage_zerod_MAC_cutoff,
                       double & t_MAC,
                       arma::uvec & t_indexZero,
                       arma::uvec & t_indexNonZero,
                       arma::uvec & t_sampleSubIndices) // Parameter for subsetting indices
{
    bool flip = false;
    //bool ischange = false;
    t_indexNonZero.clear();
    t_indexZero.clear();
    int nMissing = t_indexForMissing_arma.size();
    //uint dosagesSize = t_GVec.size();
    //double imputeG = 0;
    //double t_altFreq_0 = t_altFreq;
    // Subset t_GVec using t_sampleSubIndices and store in t_GVec_sub
    if (!t_sampleSubIndices.is_empty()) {
        t_GVec_sub = t_GVec.elem(t_sampleSubIndices); // Store subset in t_GVec_sub
	//ischange = true;
    } else {
        t_GVec_sub = t_GVec; // If no subsetting, copy the original vector
    }
    
    arma::uvec m_indexForMissing_arma_sub;
    double imputeG = 0;
    if (nMissing > 0) {
	if((!t_sampleSubIndices.is_empty())){
    		for (size_t i = 0; i < t_sampleSubIndices.n_elem; i++) {
        		if (arma::any(t_indexForMissing_arma == t_sampleSubIndices[i])) {
            			m_indexForMissing_arma_sub.insert_rows(m_indexForMissing_arma_sub.n_elem, arma::uvec{static_cast<unsigned int>(i)});
        		}
    		}
	}else{
		m_indexForMissing_arma_sub = t_indexForMissing_arma; 		         }





    if(m_indexForMissing_arma_sub.n_elem > 0){
		//ischange = true;
    		arma::uvec valid_indices = arma::regspace<arma::uvec>(0, t_GVec_sub.n_elem - 1); // All indices
    		valid_indices.shed_rows(m_indexForMissing_arma_sub); // Remove indices for missing elements
    		arma::vec valid_elements = t_GVec_sub.elem(valid_indices); // Extract valid elements
    		t_altFreq = arma::mean(valid_elements) / 2.0; // Compute mean and divide by 2
	    if(t_altFreq > 0.5){
		t_GVec_sub = 2 - t_GVec_sub;
		t_altFreq = 1 - t_altFreq;
		flip = true;
	    }

            switch (string_to_case.at(t_impute_method)) {
                case 1:
                    imputeG = std::round(2 * t_altFreq);
                    break;
                case 2:
                    imputeG = 2 * t_altFreq;
                    break;
                case 3:
                    imputeG = 0;
                    break;
            }

	    t_GVec_sub.elem(m_indexForMissing_arma_sub).fill(imputeG);
	    //if(imputeG != 0){
	    //}
	}//if(m_indexForMissing_arma_sub.n_elem > 0){
        //}	
	//else{	    	
        t_altCount = arma::accu(t_GVec_sub); 
        t_altFreq = t_altCount / (2 * (t_GVec_sub.n_elem));
        t_MAC = arma::min(arma::vec({t_altCount, 2*(t_GVec_sub.n_elem)-t_altCount}));	    
     }else{//if (nMissing > 0) {

std::cout << "t_altCount inner " << t_altCount << std::endl;
std::cout << "t_altFreq inner " << t_altFreq << std::endl;
        t_altCount = arma::accu(t_GVec_sub);
	t_altFreq = t_altCount / (2 * (t_GVec_sub.n_elem));
	 if(t_altFreq > 0.5){
                t_GVec_sub = 2 - t_GVec_sub;
                t_altFreq = 1 - t_altFreq;
                flip = true;
            }

        t_altCount = arma::accu(t_GVec_sub); 
        t_altFreq = t_altCount / (2 * (t_GVec_sub.n_elem));
        t_MAC = arma::min(arma::vec({t_altCount, 2*(t_GVec_sub.n_elem)-t_altCount}));	

std::cout << "t_altCount inner 2 " << t_altCount << std::endl;
std::cout << "t_altFreq inner 2 " << t_altFreq << std::endl;

std::cout << "flip " << flip << std::endl;


     }
    if (t_dosage_zerod_cutoff > 0) {
        if (t_MAC <= t_dosage_zerod_MAC_cutoff) {
            t_GVec_sub.clean(t_dosage_zerod_cutoff); // Clean the subsetted vector
        // Update t_altCount and t_altFreq based on the subsetted vector
            t_altCount = arma::sum(t_GVec_sub); // Sum of t_GVec_sub elements
            t_altFreq = arma::mean(t_GVec_sub) / 2; // Mean of t_GVec_sub divided by 2
	    t_MAC = arma::min(arma::vec({t_altCount, 2*(t_GVec_sub.n_elem)-t_altCount}));	
	    //ischange = true;
        }
    }

/*
    if (t_altFreq > 0.5) {
        if(flip == true){
		flip = false;
	}else{
		flip = true;
	}
        t_GVec_sub = 2 - t_GVec_sub; // Flip the subsetted vector
        t_altFreq = 1 - t_altFreq;
        t_altCount = 2 * t_GVec_sub.n_elem - t_altCount;
	ischange = true;	
    }
*/
    // Update t_indexZero and t_indexNonZero based on the subsetted vector
    
    /*
    for (unsigned int i = 0; i < t_GVec_sub.size(); i++) {
        if (t_GVec_sub(i) == 0) {
            t_indexZero.push_back(i);
        } else {
            t_indexNonZero.push_back(i);
        }
    }*/
    //if(ischange){
    	t_indexZero = arma::find(t_GVec_sub == 0);
    	t_indexNonZero = arma::find(t_GVec_sub != 0);
    //}

    return flip;
}
