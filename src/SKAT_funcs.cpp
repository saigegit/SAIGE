#include <RcppArmadillo.h>
#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/non_central_chi_squared.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>


//#include "qfc_rcpp.hpp"
#include "SKAT_funcs.hpp"
#include "SPA.hpp"


// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]


// Extern declaration of the `qfc` function from qfc.cpp
extern "C" {
 void qfc_1(double* lb1, double* nc1, int* n1, int *r1, double *sigma, double *c1, int *lim1, double *acc, double* trace, int* ifault, double *res);
 
 //void qfc_1(double* lambdas, double* noncentral, int* df, int* r, 
   //          double* sigma, double* q, int* lim, double* acc, 
     //        double* trace, int* ifault, double* res);
}

// [[Rcpp::export]]
Rcpp::List call_qfc(arma::vec& lambdas, arma::vec& noncentral,
                    arma::ivec& df, int r, double sigma, double q, 
                    int lim, double acc 
		   		    
		    ) {
    // Allocate output vectors
    arma::vec trace(7, arma::fill::zeros);
    int ifault = 0;
    double res = 0.0;

    //std::cout << "call_qfc before" << std::endl;
    //lambdas.print("lambdas");
    arma::Col<int> dfllivec(df.n_elem);

    for (arma::uword i = 0; i < df.n_elem; ++i) {
        dfllivec[i] = static_cast<int>(df[i]);
    }


    // Call the qfc function
    //qfc_1(lambdas.memptr(), noncentral.memptr(), (int*)df.memptr(), &r, &sigma, &q, &lim, &acc, trace.memptr(), &ifault, &res);
    qfc_1(lambdas.memptr(), noncentral.memptr(), dfllivec.memptr(), &r, &sigma, &q, &lim, &acc, trace.memptr(), &ifault, &res);
    
    
    
    //std::cout << "call_qfc after" << std::endl;
    //trace.print("trace");
    //std::cout << "ifault" << ifault << std::endl;
    //std::cout << "call_qfc res " << res << std::endl;
    

    // Return results as a list
    return Rcpp::List::create(
        Rcpp::Named("trace") = trace,
        Rcpp::Named("ifault") = ifault,
        Rcpp::Named("res") = res
    );
    
}


// Function to compute filtered eigenvalues from a matrix K
arma::vec Get_Lambda(arma::mat& K, bool isFast, int maxK) {
  // Eigenvalue decomposition
  arma::vec lambda1;
  arma::mat K2 = forceSymmetric(K);
  arma::eig_sym(lambda1, K2);  // Symmetric eigenvalue decomposition
  
  // Find indices of positive eigenvalues
  arma::uvec IDX1 = arma::find(lambda1 >= 0);
  
  // If no positive eigenvalues, return empty
  if (IDX1.n_elem == 0) {
    Rcpp::stop("No Eigenvalue is bigger than 0!!");
  }
  
  // Find indices of eigenvalues greater than mean(lambda1[IDX1])/100000
  double threshold = mean(lambda1.elem(IDX1)) / 100000;
  arma::uvec IDX2 = find(lambda1 > threshold);
  
  // If no eigenvalues satisfy this condition, throw an error
  if (IDX2.n_elem == 0) {
    Rcpp::stop("No Eigenvalue is bigger than the threshold!");
  }
  
  // Extract filtered eigenvalues
  arma::vec lambda = lambda1.elem(IDX2);
  
  // Return filtered eigenvalues as arma::vec
  return lambda;
}


// Function to compute the optimal parameters
Rcpp::List SKAT_META_Optimal_Param(arma::mat& Phi, arma::vec& r_all) {
  int p_m = Phi.n_cols;
  int r_n = r_all.n_elem;

  //std::cout << "SKAT_META_Optimal_Param 1" << std::endl;
  // ZMZ calculation
  arma::vec Z_item1_1 = Phi * arma::ones<arma::vec>(p_m);
  arma::mat ZZ = Phi;
  arma::mat ZMZ = Z_item1_1 * (Z_item1_1.t()) / arma::accu(ZZ);
  //std::cout << "SKAT_META_Optimal_Param 2" << std::endl;

  // W3.2 Term: mixture chi-squared
  arma::mat W3_2_t = ZZ - ZMZ;
  arma::vec lambda = Get_Lambda(W3_2_t);
  //std::cout << "SKAT_META_Optimal_Param 3" << std::endl;
  
  // W3.3 Term: variance of remaining
  double W3_3_item = arma::accu(ZMZ % (ZZ - ZMZ)) * 4;
  //std::cout << "SKAT_META_Optimal_Param 4" << std::endl;
  
  // tau term
  double z_mean_2 = arma::accu(ZZ) / std::pow(static_cast<double>(p_m), 2);

  double tau1 = arma::accu(ZZ * ZZ) / std::pow(static_cast<double>(p_m), 2) / z_mean_2;
  //std::cout << "SKAT_META_Optimal_Param 5" << std::endl;

  // Mixture Parameters
  double MuQ = arma::accu(lambda);
  double VarQ = arma::accu(arma::square(lambda)) * 2 + W3_3_item;
  double KerQ = arma::accu(arma::pow(lambda, 4)) / std::pow(arma::accu(arma::square(lambda)), 2) * 12;
  double Df = 12 / KerQ;
  //std::cout << "SKAT_META_Optimal_Param 6" << std::endl;

  // W3.1 Term: tau1 * chi-squared_1
  arma::vec tau(r_n, arma::fill::zeros);
  double r_corr, term1;
  for (int i = 0; i < r_n; i++) {
    r_corr = r_all(i);
    term1 = std::pow(p_m, 2) * r_corr * z_mean_2 + tau1 * (1 - r_corr);
    tau(i) = term1;
  }
  //std::cout << "SKAT_META_Optimal_Param 7" << std::endl;

  // Return as a list
  return Rcpp::List::create(
    Rcpp::Named("MuQ") = MuQ,
    Rcpp::Named("VarQ") = VarQ,
    Rcpp::Named("KerQ") = KerQ,
    Rcpp::Named("lambda") = arma::sort(lambda, "descend"),
    Rcpp::Named("VarRemain") = W3_3_item,
    Rcpp::Named("Df") = Df,
    Rcpp::Named("tau") = tau,
    Rcpp::Named("z_mean_2") = z_mean_2,
    Rcpp::Named("p_m") = p_m,
    Rcpp::Named("tau_1") = tau1,
    Rcpp::Named("tau_2") = p_m * z_mean_2
  );
}


// Function to compute Liu's parameters
Rcpp::List Get_Liu_Params_Mod(arma::vec& c1) {
  
  // Extract parameters from c1
  double muQ = c1(0);
  double sigmaQ = std::sqrt(2 * c1(1));
  double s1 = c1(2) / std::pow(c1(1), 1.5);
  double s2 = c1(3) / std::pow(c1(1), 2);

  double beta1 = std::sqrt(8) * s1;
  double beta2 = 12 * s2;
  
  int type1 = 0;
  double a, d, l, muX, sigmaX;

  // Conditional check for s1^2 > s2
  if (std::pow(s1, 2) > s2) {
    a = 1 / (s1 - std::sqrt(std::pow(s1, 2) - s2));
    d = s1 * std::pow(a, 3) - std::pow(a, 2);
    l = std::pow(a, 2) - 2 * d;
  } else {
    type1 = 1;
    l = 1 / s2;
    a = std::sqrt(l);
    d = 0;
  }

  muX = l + d;
  sigmaX = std::sqrt(2) * a;

  // Return results as a list
  return Rcpp::List::create(
    Rcpp::Named("l") = l,
    Rcpp::Named("d") = d,
    Rcpp::Named("muQ") = muQ,
    Rcpp::Named("muX") = muX,
    Rcpp::Named("sigmaQ") = sigmaQ,
    Rcpp::Named("sigmaX") = sigmaX
  );
}

// Function to compute Liu's parameters for the null approximation
Rcpp::List Get_Liu_Params_Mod_Lambda(arma::vec& lambda, arma::ivec& df1) {
  
  // Set default df1 to a vector of ones if not provided
  arma::ivec df1_local = df1;
  if (df1_local.n_elem == 0) {
    df1_local = arma::ones<arma::ivec>(lambda.n_elem);
  }

  // Compute the coefficients c1
  arma::vec c1(4, arma::fill::zeros);
  for (size_t i = 0; i < 4; ++i) {
    c1(i) = arma::accu(arma::pow(lambda, i+1) % df1_local);
  }

  // Extract parameters
  double muQ = c1(0);
  double sigmaQ = std::sqrt(2 * c1(1));
  double s1 = c1(2) / std::pow(c1(1), 1.5);
  double s2 = c1(3) / std::pow(c1(1), 2);

  double beta1 = std::sqrt(8) * s1;
  double beta2 = 12 * s2;
  
  int type1 = 0;
  double a, d, l, muX, sigmaX;

  // Conditional check for s1^2 > s2
  if (std::pow(s1, 2) > s2) {
    a = 1 / (s1 - std::sqrt(std::pow(s1, 2) - s2));
    d = s1 * std::pow(a, 3) - std::pow(a, 2);
    l = std::pow(a, 2) - 2 * d;
  } else {
    type1 = 1;
    l = 1 / s2;
    a = std::sqrt(l);
    d = 0;
  }

  muX = l + d;
  sigmaX = std::sqrt(2) * a;

  // Return results as a list
  return Rcpp::List::create(
    Rcpp::Named("l") = l,
    Rcpp::Named("d") = d,
    Rcpp::Named("muQ") = muQ,
    Rcpp::Named("muX") = muX,
    Rcpp::Named("sigmaQ") = sigmaQ,
    Rcpp::Named("sigmaX") = sigmaX
  );
}

//Get_Liu_PVal.MOD.Lambda
// Function to compute p-value based on Liu's null approximation
arma::vec Get_Liu_PVal_MOD_Lambda(arma::vec& Q_all, arma::vec& lambda, arma::ivec& df1, bool log_p) {
  
  // Get parameters from Get_Liu_Params_Mod_Lambda function
  Rcpp::List param = Get_Liu_Params_Mod_Lambda(lambda, df1);
  
  // Extract parameters
  double muQ = Rcpp::as<double>(param["muQ"]);
  double sigmaQ = Rcpp::as<double>(param["sigmaQ"]);
  double muX = Rcpp::as<double>(param["muX"]);
  double sigmaX = Rcpp::as<double>(param["sigmaX"]);
  double l = Rcpp::as<double>(param["l"]);
  double d = Rcpp::as<double>(param["d"]);
  
  // Normalize Q.all
  arma::vec Q_norm = (Q_all - muQ) / sigmaQ;
  arma::vec Q_norm1 = Q_norm * sigmaX + muX;
  
  // Calculate the p-values using the chi-squared distribution
  arma::vec p_value(Q_all.n_elem);
 
  //R::pchisq(x, df, lower_tail, log_p)
  //R::pchisq(x, df, ncp, lower_tail, log_p)
  //p.value<- pchisq(Q.Norm1,  df = param$l,ncp=param$d, lower.tail=FALSE, log.p=log.p)
  for (size_t i = 0; i < Q_all.n_elem; ++i) {
    boost::math::non_central_chi_squared dist(l, d);
    if(Q_norm1(i) > 0){
    	p_value(i) = 1 - boost::math::cdf(dist, Q_norm1(i));
    }else{
	p_value(i) = 1;
    }
  }
  
  return p_value; // Return as arma::vec
}

// Function to calculate p-value with the chi-squared distribution and compare it to a predefined set of quantiles
std::string Get_Liu_PVal_MOD_Lambda_Zero(double Q, double muQ, double muX, double sigmaQ, double sigmaX, double l, double d) {
  
  // Normalize the Q values
  double Q_norm = (Q - muQ) / sigmaQ;
  double Q_norm1 = Q_norm * sigmaX + muX;
  
  // Predefined quantiles (as in the R version)
  arma::vec temp = {0.05, 1e-10, 1e-20, 1e-30, 1e-40, 1e-50, 1e-60, 1e-70, 1e-80, 1e-90, 1e-100};
  
  // Compute the chi-squared quantiles
  arma::vec out(temp.n_elem);

    boost::math::non_central_chi_squared dist(l, d);

  for (size_t i = 0; i < temp.n_elem; ++i) {
    //R::pchisq(x, df, ncp, lower_tail, log_p)
    if(temp(i) > 0){
    	out(i) = 1 - boost::math::cdf(dist, temp(i)); // Lower tail = FALSE
    }else{
	out(i) = 1;
    }
  }
  
  // Find the maximum index where the quantile is less than Q_norm1
  arma::uvec IDXvec = arma::find(out < Q_norm1);
  unsigned int IDX = IDXvec.max();
  // Output the corresponding p-value message to std::cout
  std::string result = "Pvalue < " + std::to_string(temp(IDX)) + "e";
  std::cout << result << std::endl;
  return(result);
}


// Function to calculate Q for SKAT meta-analysis
arma::vec SKAT_META_Optimal_Get_Q(arma::vec& Score, arma::vec& r_all) {
  int n_r = r_all.n_elem;
  arma::vec Q_r(n_r, arma::fill::zeros);
 
  double r_corr;
  double scoresqusum = arma::accu(arma::square(Score));
  double scoresumsqu = std::pow(arma::accu(Score), 2);
  for (int i = 0; i < n_r; i++) {
    r_corr = r_all(i);
    Q_r(i) = (1 - r_corr) * scoresqusum + r_corr * scoresumsqu;
  }

  Q_r = Q_r / 2;
  return Q_r;
}

// Function to calculate Q for resampling
arma::mat SKAT_META_Optimal_Get_Q_Res(arma::mat& Score_res, arma::vec& r_all) {
  int n_r = r_all.n_elem;
  int p = Score_res.n_rows;
  arma::mat Q_r(p, n_r, arma::fill::zeros);
 
  double r_corr;
  arma::vec rsvec = arma::sum(arma::square(Score_res), 1);
  arma::vec rsvec_2 = arma::square(arma::sum(Score_res, 1));
  for (int i = 0; i < n_r; i++) {
    r_corr = r_all(i);
    Q_r.col(i) = (1 - r_corr) * rsvec  + r_corr * rsvec_2;
  }
  Q_r =  Q_r / 2; 
  return Q_r;
}





// Davies' p-value computation function
double daviesPValue(arma::vec& eigenvalues, double q, double tol) {
    arma::vec sorted_eigenvalues = arma::sort(eigenvalues, "descend");
    double p_value = 0.0;
    double cumulative = 0.0;

    for (size_t i = 0; i < sorted_eigenvalues.n_elem; ++i) {
        double lambda = sorted_eigenvalues[i];
        if (lambda <= 0) continue;

        boost::math::chi_squared chi2(1.0); 
        double contrib = 1.0 - boost::math::cdf(chi2, q / lambda);
        cumulative += contrib;
        if (cumulative >= 1.0 - tol) {
            break; 
        }
        p_value += contrib;
    }

    return std::min(p_value, 1.0); 
}

// Wrapper for Davies p-value calculation
// [[Rcpp::export]]
Rcpp::List Get_Davies_PVal(arma::mat & Q, arma::mat & W, arma::mat & Q_resampling, bool isFast) {
    // Ensure matrices are in the correct format
    arma::mat K = W / 2;
    
    // Combine Q and Q.resampling if available
    arma::vec Q_all = arma::vectorise(Q);
    arma::vec Q_resampling_vec = arma::vectorise(Q_resampling);

    if (Q_resampling_vec.n_elem > 0) {
        Q_all = arma::join_cols(Q_all, Q_resampling_vec);
    }

    // Calculate the Davies p-value based on the kernel matrix K and the combined Q values
    arma::vec eigvalK = Get_Lambda(K, false, 100);
    arma::ivec dfvec = arma::ones<arma::ivec>(eigvalK.n_elem);
    Rcpp::List re = Get_PValue_Lambda(eigvalK, Q_all, dfvec);

    Rcpp::List param;
    // Adding param information
    arma::vec p_val_liu = re["p_val_liu"];
    arma::ivec is_converge = re["is_converge"];
    param["liu.pval"] = p_val_liu(0);
    param["Is_Converged"] = is_converge(0);  // Assuming convergence is always true for simplicity
    std::string p_val_msg = re["pval_zero_msg"];

    arma::vec p_valuevec = re["p.value"];
    double p_value = p_valuevec(0);
    //std::cout << "Get_Davies_PVal p_value " << p_value << std::endl;

    // Prepare the result list
    Rcpp::List result;
    
    // If resampling is available
    arma::vec p_value_resampling;
    //if (Q_resampling.n_elem > 0) {

    result["p.value"] = p_value;
    result["param"] = param;
    result["p.value.resampling"] = p_value_resampling;
    result["pval.zero.msg"] = p_val_msg;  // Assuming no zero p-values for simplicity

    //std::cout << "Get_Davies_PVal end" << std::endl;
    return result;
}



// The SKAT_davies function using RcppArmadillo
// [[Rcpp::export]]
Rcpp::List SKAT_davies(double q, arma::vec& lambda, arma::ivec& h, arma::vec& delta, double sigma, int lim, double acc) {
  
  // Ensure lambda, h, and delta have the same length
  int r = lambda.n_elem;

  if(h.n_elem == 0){
    h.ones(r); 
  }

  if(delta.n_elem == 0){
    delta.zeros(r); 
  }

  if (h.n_elem != r) Rcpp::stop("lambda and h should have the same length!");
  if (delta.n_elem != r) Rcpp::stop("lambda and delta should have the same length!");
  
  // Initialize result variables
  arma::vec trace(7, arma::fill::zeros); // Trace vector to store intermediate results
  int ifault = 0; // Fault code
  double res = 0.0; // Result of the computation
 
  //std::cout << "Call the qfc_1 function before" << std::endl;
  // Call the qfc_1 function (from the C code)
  //QUADFORM::QuadFormClass* quadForm = new QUADFORM::QuadFormClass(r, lim);
  //std::cout << "Call the qfc_1 function before a" << std::endl;
 /*
  Rcpp::List call_qfc(const arma::vec& lambdas, const arma::vec& noncentral,
                    const arma::ivec& df, int r, double sigma, double q,
		                        int lim, double acc)
*/


  //lambda.print("lambda");
  //delta.print("delta");
  //h.print("h");
  //std::cout << "r sigma q lim acc " << r << " " << sigma << " " << q << " " << lim << " " << acc << std::endl;

/*
  lambda[0] = 5144.902;
  lambda[1] = 4946.847;
  lambda[2] = 4080.661;
  lambda[3] = 3484.483;

  q = 3294.954;
*/
//acc = 1e-5; 
  Rcpp::List results = call_qfc(lambda, delta, h, r, sigma, q, lim, acc);
  //qfc(lambda, delta, h, r, sigma, q, lim, acc, trace, ifault, res);
  //std::cout << "Call the qfc_1 function before b" << std::endl;

  //delete quadForm;

  // Adjust result
  res = results["res"];
  res = 1.0 - res;
  //std::cout << "res " << res << std::endl;
  ifault = results["ifault"];
  //std::cout << "ifault " << ifault << std::endl;
  // Return results as a list
  return Rcpp::List::create(
    Rcpp::Named("trace") = results["trace"],
    Rcpp::Named("ifault") = ifault,
    Rcpp::Named("Qq") = res
  );
}


Rcpp::List Get_PValue_Lambda(arma::vec &lambda, arma::vec &Q, arma::ivec &df1){
  
  int n1 = Q.n_elem;
  
  arma::vec p_val = arma::zeros<arma::vec>(n1);
  arma::vec p_val_liu = arma::zeros<arma::vec>(n1);
  arma::ivec is_converge = arma::zeros<arma::ivec>(n1);
  
  // Get Liu p-values
  //std::cout << "Get_Liu_PVal_MOD_Lambda 1" << std::endl;
  p_val_liu = Get_Liu_PVal_MOD_Lambda(Q, lambda, df1);
  //std::cout << "Get_Liu_PVal_MOD_Lambda 2" << std::endl;
 
  double Qi;
  arma::ivec h;
  arma::vec delta;
  Rcpp::List out;
  //df1.print("df1");
  for (int i = 0; i < n1; i++) {
    Qi = Q(i);
    if (df1.n_elem == 0) {
      // Call SKAT_davies for each Q[i] and lambda
      out = SKAT_davies(Qi, lambda, h, delta, 0, 10000, std::pow(10, -6));
    } else {
      // Call SKAT_davies with h = df1 for each Q[i] and lambda
      out = SKAT_davies(Qi, lambda, df1, delta, 0, 1000, std::pow(10, -6));
    }
    p_val(i) = out["Qq"];
    
    is_converge(i) = 1;
    
    // Check convergence
    if (lambda.n_elem == 1) {
      p_val[i] = p_val_liu[i];
    } else if (out["ifault"] != 0) {
      is_converge[i] = 0;
    }
    
    // Check p-value
    if (p_val[i] > 1 || p_val[i] <= 0) {
      is_converge[i] = 0;
      p_val[i] = p_val_liu[i];
    }
  }
 
    //std::cout << "Get_Liu_PVal_MOD_Lambda 3" << std::endl; 
  // For p-value zero case
  //Rcpp::Nullable<double> p_val_msg = R_NilValue;
  //Rcpp::Nullable<double> p_val_log = R_NilValue;
  std::string p_val_msg;
  arma::vec p_val_log;
  arma::vec  Q0 = Q.subvec(0,0);


    //std::cout << "Get_Liu_PVal_MOD_Lambda 4" << std::endl; 
  if (p_val[0] == 0) {
    // Get Liu Params Mod Lambda
    Rcpp::List param = Get_Liu_Params_Mod_Lambda(lambda, df1);
    //std::string Get_Liu_PVal_MOD_Lambda_Zero(double Q, double muQ, double muX, double sigmaQ, double sigmaX, double l, double d) {
    p_val_msg = Get_Liu_PVal_MOD_Lambda_Zero(Q(0), param["muQ"], param["muX"], param["sigmaQ"], param["sigmaX"], param["l"], param["d"]);
    p_val_log = Get_Liu_PVal_MOD_Lambda(Q0, lambda, df1, true);
  }else{
    p_val_msg="";
  }
    //std::cout << "Get_Liu_PVal_MOD_Lambda 5" << std::endl; 
 
//p_val.print("p_val");
//p_val_liu.print("p_val_liu");
//p_val_log.print("p_val_log");
double p_val_log0;
if(p_val_log.n_elem > 0){
  p_val_log0 = p_val_log[0];
}else{
  p_val_log0 = 0;
}
  // Return the results as a List
  return Rcpp::List::create(
    Rcpp::Named("p.value") = p_val,
    Rcpp::Named("p_val_liu") = p_val_liu,
    Rcpp::Named("is_converge") = is_converge,
    Rcpp::Named("p_val_log") = p_val_log0,
    Rcpp::Named("pval_zero_msg") = p_val_msg
  );
}

Rcpp::List SKAT_Optimal_Each_Q(Rcpp::List &param_m, arma::mat &Q_all, arma::vec &r_all, Rcpp::List & lambda_all, std::string & method) {
  
  int n_r = r_all.n_elem;
  arma::vec c1 = arma::zeros<arma::vec>(4);
  //Q_all.print("Q_all") ;
  int n_q = Q_all.n_rows;
  // Initialize matrices for p-values and minimum p-values
  arma::mat pval = arma::zeros<arma::mat>(n_q, n_r);
  arma::mat pmin_q = arma::zeros<arma::mat>(n_q, n_r);
  arma::mat param_mat(n_r, 3);
  arma::vec Q, Q_norm;
  //arma::vec lambda_temp;
  double r_corr, muQ, varQ, df, sigmaQ;
  arma::ivec empvec;
  //std::cout << "SKAT_Optimal_Each_Q 0" << std::endl;
  for (int i = 0; i < n_r; i++) {
    //std::cout << "SKAT_Optimal_Each_Q 001" << std::endl;
    Q = Q_all.col(i);
    //std::cout << "SKAT_Optimal_Each_Q 002" << std::endl;
    r_corr = r_all[i];
    //std::cout << "SKAT_Optimal_Each_Q 003" << std::endl;
    arma::vec lambda_temp = arma::vec(Rcpp::as<arma::vec>(lambda_all[i]));
    //std::cout << "SKAT_Optimal_Each_Q 004" << std::endl; 
    // Compute c1 components
    c1[0] = arma::accu(lambda_temp);
    c1[1] = arma::accu(arma::square(lambda_temp));
    c1[2] = arma::accu(arma::pow(lambda_temp, 3));
    c1[3] = arma::accu(arma::pow(lambda_temp, 4));
    
    //std::cout << "SKAT_Optimal_Each_Q 00" << std::endl;
    //std::cout << "i " << i << std::endl;
    // Get parameters
    Rcpp::List param_temp = Get_Liu_Params_Mod(c1);
    //std::cout << "Get_Liu_Params_Mod after" << std::endl;
    
    muQ = param_temp["muQ"];
    sigmaQ = param_temp["sigmaQ"];
    varQ = std::pow(sigmaQ, 2);
    df = param_temp["l"];
    boost::math::chi_squared dist(df);  
    // Get p-value
    //std::cout << "SKAT_Optimal_Each_Q 01" << std::endl;
    
    Q_norm = ((Q - muQ) / std::sqrt(varQ)) * std::sqrt(2 * df) + df;
    //Q_norm.print("Q_norm");
    //pval.col(i) = R::pchisq(Q_norm, df, 0, false, false);
    //p_value(i) = 1 - boost::math::cdf(dist, Q_norm1(i));
    //std::cout << "SKAT_Optimal_Each_Q 02" << std::endl;
    for (size_t j = 0; j < Q_norm.n_elem; ++j) {
      if(Q_norm(j) > 0){
        pval(j,i) = 1 - boost::math::cdf(dist, Q_norm(j));
      }else{
        pval(j,i) = 1;
      }
     }
    //std::cout << "SKAT_Optimal_Each_Q 1" << std::endl;
    Rcpp::List pval_lambda;
    std::string method_str;
    // Check for method-specific changes
    if (method != "") {
      method_str = method;
      //std::cout << "method_str " << method_str << std::endl;
      if (method_str == "optimal.mod" || method_str == "optimal.adj" || method_str == "optimal.moment.adj") {
        pval_lambda = Get_PValue_Lambda(lambda_temp, Q, empvec);
	//std::cout << "SKAT_Optimal_Each_Q 1a" << std::endl;
        pval.col(i) = arma::vec(Rcpp::as<arma::vec>(pval_lambda["p.value"]));
      }
    }
    //std::cout << "SKAT_Optimal_Each_Q 2" << std::endl;
    //arma::mat tempmat(1,3);
    //tempmat(0,0) = muQ;
    //tempmat(0,1) = varQ;
    //tempmat(0,2) = df;
    //param_mat.print("param_mat");
    //std::cout << "muQ "<< muQ << std::endl;
    //std::cout << "varQ "<< varQ << std::endl;
    //std::cout << "df "<< df << std::endl;
    param_mat(i,0) = muQ;
    param_mat(i,1) = varQ;
    param_mat(i,2) = df;
    //std::cout << "SKAT_Optimal_Each_Q 2b" << std::endl;
  }
    //std::cout << "SKAT_Optimal_Each_Q 3" << std::endl;
  
  // Calculate the minimum p-values
  arma::vec pmin = arma::min(pval, 1);
  arma::vec q_org, q_q;
  q_org.set_size(pmin.n_elem);
  //param_mat.print("param_mat");
  //pmin.print("pmin");
  for (int i = 0; i < n_r; i++) {
    muQ = param_mat(i, 0);
    varQ = param_mat(i, 1);
    df = param_mat(i, 2);
    boost::math::chi_squared dist(df);    
    for (int j = 0; j < pmin.n_elem; j++) {
        q_org[j] =  boost::math::quantile(dist, 1-pmin[j]);
    }
    //std::cout << "i " << i << std::endl;
    q_q = (q_org - df) / std::sqrt(2 * df) * std::sqrt(varQ) + muQ;
    pmin_q.col(i) = q_q;
  }
  //std::cout << "SKAT_Optimal_Each_Q 4" << std::endl;
  //pval.print("pval");
  // Return the results as a List
  return Rcpp::List::create(
    Rcpp::Named("pmin") = pmin,
    Rcpp::Named("pval") = pval,
    Rcpp::Named("pmin_q") = pmin_q
  );
}


// Davies Integration Function for SKAT
// [[Rcpp::export]]
arma::vec SKAT_Optimal_Integrate_Func_Davies(arma::vec & x,  arma::mat &pmin_q, Rcpp::List &param_m,  arma::vec &r_all) {
  int n_r = r_all.n_elem;
  int n_x = x.n_elem;;

    arma::mat taumat = param_m["tau"];
      arma::mat xt(1, n_x);
        xt.row(0) = x.t();
 //taumat.print("taumat");
 //xt.print("xt");
	  arma::mat temp1 = arma::kron(taumat ,  xt);

 //temp1.print("temp1");
 //std::cout << "temp1.n_rows " << temp1.n_rows << std::endl;
 //std::cout << "temp1.n_cols " << temp1.n_cols << std::endl;
arma::mat temp_sub = arma::repmat(pmin_q, 1, temp1.n_cols) - temp1;  // Replicate pmin_q to 7x2000
arma::mat r_all_expanded = arma::repmat(r_all, 1, temp1.n_cols);  // Replicate r_all to 7x2000
arma::mat tempmat = temp_sub / (1 - r_all_expanded);
//tempmat.print("tempmat");
 arma::rowvec temp_min = arma::min(tempmat, 0);
 //temp_min.print("temp_min");
  arma::vec re = arma::zeros<arma::vec>(n_x);
  arma::vec lambda = param_m["lambda"];
  double min1, min1_temp, sd1, min1_st, xi, VarQ, VarRemain, temp, MuQ;
  arma::ivec h; 
  VarQ = param_m["VarQ"];
  VarRemain = param_m["VarRemain"];
  MuQ = param_m["MuQ"];
  int ifault = 0;
  for (int i = 0; i < n_x; i++) {
    //std::cout << "i " << i << std::endl;
    min1 = temp_min[i];
    if (min1 > arma::accu(lambda) * 1e4) {
      temp = 0;
    } else {
       //std::cout << "min1 " << min1 << std::endl;
       //std::cout << "MuQ " << MuQ << std::endl;
       //std::cout << "VarRemain " << VarRemain << std::endl;
       //std::cout << "VarQ " << VarQ << std::endl;
      min1_temp = min1 - MuQ;
      sd1 = std::sqrt(VarQ - VarRemain) / std::sqrt(VarQ);
      min1_st = min1_temp * sd1 + MuQ;
      arma::vec deltavec = arma::zeros<arma::vec>(lambda.n_elem);
      Rcpp::List dav_re = SKAT_davies(min1_st, lambda, h, deltavec, 0, 10000, 1e-6);
      temp = dav_re["Qq"];
      ifault = dav_re["ifault"];
      if (ifault != 0) {
        Rcpp::stop("dav_re$ifault is not 0");
      }
    }
    if (temp > 1) {
      temp = 1;
    }
    xi = x[i];

    boost::math::chi_squared dist1(1);
    re[i] = (1 - temp) * boost::math::pdf(dist1, x(i)) ;
  }
  return re;
}

// Wrapper class to integrate over x using Boost's Gauss-Kronrod quadrature
class DaviesFunction {
    arma::mat &pmin_q;
    Rcpp::List &param_m;
    arma::vec &r_all;

public:
    DaviesFunction(arma::mat &pmin_q_, Rcpp::List &param_m_, arma::vec &r_all_)
        : pmin_q(pmin_q_), param_m(param_m_), r_all(r_all_) {}

    double operator()(double x) const {
        // Convert scalar x to vector form for SKAT_Optimal_Integrate_Func_Davies
        arma::vec x_vec(1, arma::fill::value(x));
        arma::vec result = SKAT_Optimal_Integrate_Func_Davies(x_vec, pmin_q, param_m, r_all);
        return result[0]; // Assuming scalar output for each input x
    }
};


// Wrapper class to integrate over x using Boost's Gauss-Kronrod quadrature
class LiuFunction {
    arma::mat &pmin_q;
    Rcpp::List &param_m;
    arma::vec &r_all;

public:
    LiuFunction(arma::mat &pmin_q_, Rcpp::List &param_m_, arma::vec &r_all_)
        : pmin_q(pmin_q_), param_m(param_m_), r_all(r_all_) {}

    double operator()(double x) const {
        // Convert scalar x to vector form for SKAT_Optimal_Integrate_Func_Davies
        arma::vec x_vec(1, arma::fill::value(x));
        arma::vec result = SKAT_Optimal_Integrate_Func_Liu(x_vec, pmin_q, param_m, r_all);
        return result[0]; // Assuming scalar output for each input x
    }
};


// Perform adaptive integration over x using Gauss-Kronrod
double integrate_SKAT_Optimal_Davies(arma::vec &pmin_q, Rcpp::List &param_m, arma::vec &r_all, double lower, double upper, int subdivisions, double abs_tol){
    DaviesFunction func(pmin_q, param_m, r_all);

    // Perform integration with adaptive Gauss-Kronrod quadrature (15-point rule)
   try {
    boost::math::quadrature::gauss_kronrod<double, 200> integrator;
    //double tolerance = 1e-25; // Error tolerance
    // Integrate from lower to upper bound
    double result = integrator.integrate(func, lower, upper, abs_tol);
    //std::cout << "integrate_SKAT_Optimal_Davies result line 823 " << result << std::endl;
    return result;
   } catch (const std::exception& e) {
   	//std::cerr << "Integration error: " << e.what() << std::endl;
	return -1.0;                 
    }
}

double integrate_SKAT_Optimal_Davies_v2(arma::vec &pmin_q, Rcpp::List &param_m, arma::vec &r_all, double lower, double upper, int subdivisions, double abs_tol){
    arma::vec x = arma::linspace(lower, upper, subdivisions + 1); // Points in the range
     arma::vec y(subdivisions + 1);
	y = SKAT_Optimal_Integrate_Func_Davies(x, pmin_q, param_m, r_all);
        double result = 0.5 * (y[0] + y[y.n_elem - 1]);
        result = result + arma::accu(y);
        result *= (upper - lower) / subdivisions;

        // Check if the result is below the tolerance
        //if (result < abs_tol) {
        //    result = 0;
        //}else{
	//    result = -1.0;	
	//}		    
	return(result);
}

/*
// Define the trapezoidal rule integration function
double integrate_SKAT_Optimal_Davies(arma::vec &pmin_q, Rcpp::List &param_m, arma::vec &r_all, double lower, double upper, int subdivisions, double abs_tol) {
    double integral = 0.0;  // Initialize the result of the integral

    // Calculate the width of each subdivision
    double step_size = (upper - lower) / subdivisions;

    // Perform the trapezoidal rule integration
    double last_integral = 0.0;  // Store the previous integral value for convergence checking
    bool issuccess = false;
    for (int i = 1; i < subdivisions; ++i) {
        arma::vec x1 = arma::vec(1, arma::fill::value(lower + i * step_size));
        arma::vec x2 = arma::vec(1, arma::fill::value(lower + (i + 1) * step_size));  // Vector for x2

        // Compute the values of the function at x1 and x2
       arma::vec f1vec = SKAT_Optimal_Integrate_Func_Davies(x1, pmin_q, param_m, r_all);
        arma::vec f2vec = SKAT_Optimal_Integrate_Func_Davies(x2, pmin_q, param_m, r_all);
	double f1 = f1vec[0];
	double f2 = f2vec[0];
        // Compute the area of the trapezoid and accumulate the result
        integral += 0.5 * (f1 + f2) * step_size;

        // Check for convergence using absolute tolerance
        if (std::abs(integral - last_integral) < abs_tol) {
	    issuccess = true;
            break;  // Stop if the change is less than the tolerance
        }

        // Store the current integral value for the next iteration
        last_integral = integral;
    }

    if(!issuccess){
	integral = -1.0;
	std::cout << "An unknown error occurred during integration." << std::endl;
	}
    // Return the final integral value
    return integral;
}



// Numerical integration using the trapezoidal rule over a vector
// [[Rcpp::export]]
arma::vec integrate_SKAT_Optimal_Davies(arma::vec &pmin_q, Rcpp::List &param_m, arma::vec &r_all, double lower, double upper, int subdivisions, double abs_tol) {
    int n_r = r_all.n_elem;
    arma::vec results = arma::zeros<arma::vec>(n_r);

    arma::vec x = arma::linspace(lower, upper, subdivisions + 1); // Points in the range

    for (int i = 0; i < n_r; i++) {
        // Use each row of r_all to compute the integral
        arma::vec r_vec = r_all.subvec(i, i);
        arma::vec y = SKAT_Optimal_Integrate_Func_Davies(x, pmin_q, param_m, r_vec);

        // Apply trapezoidal rule for numerical integration
        double result = 0.5 * (y[0] + y[y.n_elem - 1]);
        for (int j = 1; j < y.n_elem - 1; j++) {
            result += y[j];
        }
        result *= (upper - lower) / subdivisions;

        // Check if the result is below the tolerance
        if (result < abs_tol) {
            result = 0;
        }
        results[i] = result;
    }
    return results;
}
*/


// Function for computing the Davies p-value
double SKAT_Optimal_PValue_Davies(arma::vec &pmin_q, Rcpp::List &param_m,
                                     arma::vec &r_all, double pmin) {
    double pvalue;
    //try {
        // Perform the integration
	//
	//std::cout << "integrate_SKAT_Optimal_Davies before" << std::endl;
	//pmin_q.print("pmin_q");
	//param_m.print("param_m");
	//r_all.print("r_all");
        double re = integrate_SKAT_Optimal_Davies(pmin_q, param_m, r_all, 0, 40, 2000, 1e-25);
        //double re = integrate_SKAT_Optimal_Davies_v2(pmin_q, param_m, r_all, 0, 40, 2000, 1e-25);
	//std::cout << "integrate_SKAT_Optimal_Davies after" << std::endl;
	//std::cout << "re" << re << std::endl;
        // Compute the p-value
	if(re > 0.0){
        	pvalue = 1 - re;
        	if (pmin != arma::datum::nan) {
            		if ((pmin * (r_all.n_elem)) < pvalue) {
                		pvalue = pmin * r_all.n_elem;
            	}
		}
        	
	}else{
        // Return the result
    //} catch (const std::exception &e) {
        // Log the error message for debugging
        	//Rcpp::Rcerr << "Error in SKAT_Optimal_PValue_Davies: " << std::endl;

        // Fall back to Liu method
        	pvalue = SKAT_Optimal_PValue_Liu(pmin_q, param_m, r_all, pmin);
	}
	//std::cout << "pmin " << pmin << std::endl;
	//std::cout << "pvalue " << pvalue << std::endl;

	return pvalue;
    //}
}

//arma::vec & x,  arma::mat &pmin_q, Rcpp::List &param_m,  arma::vec &r_all
// Liu Integration Function for SKAT
arma::vec SKAT_Optimal_Integrate_Func_Liu(arma::vec & x,  arma::mat &pmin_q,  Rcpp::List &param_m,  arma::vec &r_all) {
  int n_r = r_all.n_elem;
  int n_x = x.n_elem;

  arma::mat taumat = param_m["tau"];
//  x.print("x");
  arma::mat xt(1, n_x);
  xt.row(0) = x.t(); 
//  xt.print("xt SKAT_Optimal_Integrate_Func_Liu");
//  taumat.print("taumat SKAT_Optimal_Integrate_Func_Liu");
  arma::mat temp1 = arma::kron(taumat , xt);

  //temp1.print("temp1 SKAT_Optimal_Integrate_Func_Liu");
  //pmin_q.print("pmin_q SKAT_Optimal_Integrate_Func_Liu");
  //r_all.print("r_all SKAT_Optimal_Integrate_Func_Liu");
arma::mat temp_sub = arma::repmat(pmin_q, 1, temp1.n_cols) - temp1;  // Replicate pmin_q to 7x2000
arma::mat r_all_expanded = arma::repmat(r_all, 1, temp1.n_cols);  // Replicate r_all to 7x2000
arma::mat temp = temp_sub / (1 - r_all_expanded); 

  //arma::mat temp = (pmin_q - temp1.each_col()) / (1 - r_all);
  //temp.print("temp SKAT_Optimal_Integrate_Func_Liu");
  //std::cout << "temp.n_cols " << temp.n_cols << std::endl;
  //std::cout << "temp.n_rows " << temp.n_rows << std::endl;
 arma::rowvec temp_min = arma::min(temp, 0);

  //temp_min.print("temp_min SKAT_Optimal_Integrate_Func_Liu");
  
  
  double MuQ, VarQ, Df;
  MuQ = param_m["MuQ"];
  VarQ = param_m["VarQ"];
  Df = param_m["Df"];
  
  arma::vec temp_q(n_x);
  arma::vec re(n_x);
  boost::math::chi_squared dist(Df);
  boost::math::chi_squared dist1(1);
  for (int i = 0; i < n_x; i++) {
	//std::cout << "i " << i << std::endl;
	double temp_min_val = temp_min(i);
	double  temp_q_val = (temp_min_val - MuQ)/(std::sqrt(VarQ)) * (std::sqrt(2*Df)) + Df;
	double temp_cdf, temp_pdf;
	if(temp_q_val <= 0){
		temp_cdf = 0;
	}else{
		temp_cdf = boost::math::cdf(dist, temp_q_val);
	}
	if(x(i) <= 0){
		temp_pdf = 0;
	}else{
		temp_pdf =  boost::math::pdf(dist1, x(i));
	}
	re(i) = temp_cdf * temp_pdf;	
  }
  
  return re;
}



// Numerical integration using the trapezoidal rule over a vector
double integrate_SKAT_Optimal_Liu_v2(arma::vec &pmin_q, Rcpp::List &param_m, arma::vec &r_all, double lower, double upper, int subdivisions, double abs_tol){

    arma::vec x = arma::linspace(lower, upper, subdivisions + 1); // Points in the range
     arma::vec y(subdivisions + 1);
     
        y = SKAT_Optimal_Integrate_Func_Liu(x, pmin_q, param_m, r_all);
	//std::cout << "integrate_SKAT_Optimal_Liu_v2 b" << std::endl;

        // Apply trapezoidal rule for numerical integration
        double result = 0.5 * (y[0] + y[y.n_elem - 1]);
	result = result + arma::accu(y);
        result *= (upper - lower) / subdivisions;

        // Check if the result is below the tolerance
        if (result < abs_tol) {
            result = 0;
        }
    return result;
}


double integrate_SKAT_Optimal_Liu_v3(arma::vec &pmin_q, Rcpp::List &param_m, arma::vec &r_all, double lower, double upper, int subdivisions, double abs_tol){
    
    // Initialize variables for integration
    double previous_area = 0;
    double area = 1;
    size_t n_points = 21;  // Initial number of points for integration
    double result = -1.0;   // Default result for failure
    
    // Loop until the area converges within the given tolerance
    while (area ==0.0 || std::abs(area - previous_area) > abs_tol) {
        
	// Update previous area for the next iteration
        previous_area = area;
        
        // Create the vector of x values (points) for integration
        arma::vec x = arma::linspace(lower, upper, n_points);
        arma::vec y(n_points);

        // Compute the y values for the function at each x
        y = SKAT_Optimal_Integrate_Func_Liu(x, pmin_q, param_m, r_all);
	//x.print("x");
	//y.print("y");
        // Perform numerical integration using the trapezoidal rule
        arma::mat areamat = arma::trapz(x, y);
	//areamat.print("areamat");
	area = areamat(0,0);
        // Output debug information
        std::cout << "Integration with " << n_points << " points, area: " << area << std::endl;

        // Check if the number of points exceeds the maximum allowed subdivisions
        if (n_points > subdivisions && (std::abs(area - previous_area) > abs_tol && area != 0.0)) {
            std::cout << "Maximum subdivisions reached, integration failed to converge." << std::endl;
            return -1.0;  // Return -1.0 to indicate failure
        }

	


        // Double the number of points for better precision
        n_points *= 2;
    }

    // Return the successfully computed area
    result = area;
    return result;
}

// Perform adaptive integration over x using Gauss-Kronrod
double integrate_SKAT_Optimal_Liu(arma::vec &pmin_q, Rcpp::List &param_m, arma::vec &r_all, double lower, double upper, int subdivisions, double abs_tol){
    LiuFunction func(pmin_q, param_m, r_all);

    // Perform integration with adaptive Gauss-Kronrod quadrature (15-point rule)
    boost::math::quadrature::gauss_kronrod<double, 200> integrator;
    //double tolerance = 1e-25; // Error tolerance
    // Integrate from lower to upper bound
    double result = integrator.integrate(func, lower, upper, abs_tol);
    
    
    //double result = integrate_SKAT_Optimal_Liu_v2(pmin_q, param_m, r_all, lower, upper, subdivisions, abs_tol);
    //double result = integrate_SKAT_Optimal_Liu_v3(pmin_q, param_m, r_all, lower, upper, subdivisions, abs_tol);
    //std::cout << "result integrate_SKAT_Optimal_Liu " << result << std::endl;

    return result;
}




// Function for computing the Liu p-value
double SKAT_Optimal_PValue_Liu(arma::vec &pmin_q, Rcpp::List &param_m,
                                     arma::vec &r_all, double pmin){
  double re = integrate_SKAT_Optimal_Liu(pmin_q, param_m, r_all, 0, 40, 2000, 1e-25);


  double pvalue = 1 - re;

  if (pmin != arma::datum::nan){
    if ((pmin * r_all.n_elem) < pvalue) {
      pvalue = pmin * r_all.n_elem;
    }
  }

  return pvalue;
}

Rcpp::List SKAT_META_Optimal_Get_Pvalue( arma::mat &Q_all,  arma::mat &Phi,  arma::vec &r_all, std::string &method, bool isFast) {

  int n_r = r_all.n_elem;
  int n_q = Q_all.n_rows;
  int p_m = Phi.n_cols;

  // Initialize lambda.all as a list
  Rcpp::List lambda_all(n_r);

  double r_corr;
  arma::mat R_M, L, Phi_rho;
  bool success;
  arma::uvec pivot;
  //r_all.print("r_all");
  //std::cout << "p_m " << p_m << std::endl;
  //Phi.print("Phi");
  for (int i = 0; i < n_r; i++) {
     //std::cout << "i " << i << std::endl;
     r_corr = r_all[i];
     //std::cout << "r_corr " << r_corr << std::endl; 
     R_M = diagmat(arma::ones<arma::vec>(p_m) * (1 - r_corr)) + arma::ones<arma::mat>(p_m, p_m) * r_corr;
     //R_M.print("R_M");
     success = arma::chol(L, R_M, "upper");
     //L.print("L");
     Phi_rho = L * (Phi * L.t());
     //Phi_rho.print("Phi_rho");
     arma::vec lambda_all_i = Get_Lambda(Phi_rho, isFast);
     arma::vec lambda_all_i_sort = arma::sort(lambda_all_i, "descend");
     lambda_all[i] = lambda_all_i_sort;
     //lambda_all_i_sort.print("lambda_all_i_sort");
  }

  //std::cout << "SKAT_META_Optimal_Param before" << std::endl;
  // Get Mixture parameters
  Rcpp::List param_m = SKAT_META_Optimal_Param(Phi, r_all);
  //std::cout << "SKAT_META_Optimal_Param after" << std::endl;


  // Call SKAT_Optimal_Each_Q
  Rcpp::List Each_Info = SKAT_Optimal_Each_Q(param_m, Q_all, r_all, lambda_all, method);
  //std::cout << "SKAT_Optimal_Each_Q after" << std::endl;
  arma::mat pmin_q = Each_Info["pmin_q"];
  arma::vec pval = arma::zeros<arma::vec>(n_q);

  // Get the pmin values
  arma::vec pmin = Each_Info["pmin"];

  // Calculate p-values based on the chosen method
  if (method == "davies" || method == "optimal" || method == "optimal.adj" || method == "optimal.mod") {
    for (int i = 0; i < n_q; i++) {
      arma::vec pminqvec = pmin_q.row(i).t();
      pval[i] = SKAT_Optimal_PValue_Davies(pminqvec, param_m, r_all, pmin[i]);
      //std::cout << "i " << i << " pval[i] " << pval[i] << std::endl;
      //std::cout << "SKAT_Optimal_PValue_Davies after" << std::endl;
    }
  } else if (method == "liu" || method == "liu.mod") {
    for (int i = 0; i < n_q; i++) {
      arma::vec pminqvec = pmin_q.row(i).t();
      pval[i] = SKAT_Optimal_PValue_Liu(pminqvec, param_m, r_all, pmin[i]);
    }
  } else {
    Rcpp::stop("Invalid Method: " + method);
  }

  // Correct p-values conservatively
  int multi = 3;
  if (r_all.n_elem < 3) {
    multi = 2;
  }

  //pval.print("pval here");
  arma::mat pvalmat = Each_Info["pval"];
  int npval = pvalmat.n_rows;
  arma::vec pval_each(npval);
  arma::uvec IDX;
  double pval1;
  arma::mat pvalmatt = pvalmat.t();
  for (int i = 0; i < n_q; i++) {
    IDX.clear();
    pval_each = pvalmatt.col(i);
    IDX = arma::find(pval_each > 0);
    //pval_each.print("pval_each");
    //std::cout << "i n_q " << std::endl;
    //IDX.print("IDX");
    //r_all.print(" r_all");
    pval1 = arma::min(pval_each) * multi;
    if (pval[i] <= 0 || IDX.n_elem < r_all.n_elem) {
      pval[i] = pval1;
    }

    // If pval == 0, use non-zero minimum each.pval as p-value
    if (pval[i] == 0) {
      if (IDX.n_elem > 0) {
        pval[i] = arma::min(pval_each(IDX));
      }
    }
  }

  //pvalmat.print("pvalmat");
  //pval.print("pval");
  //std::cout << "the end of SKAT_META_Optimal_Get_Pvalue" << std::endl;
  return Rcpp::List::create(Rcpp::Named("p.value") = pval, Rcpp::Named("p.val.each") = pvalmat);
}

Rcpp::List SKAT_META_Optimal( arma::vec &Score,  arma::mat &Phi,  arma::vec &r_all, 
                              std::string method,  arma::mat &Score_Resampling, bool isFast) {

  //r_all.print("r_all herea0");
  // Adjust r_all values greater than or equal to 0.999
  arma::uvec IDX = find(r_all >= 0.999);
  if (IDX.n_elem > 0) {
    r_all(IDX).fill(0.999);
  }

  int p_m = Phi.n_cols;
  int n_r = r_all.n_elem;

  // Compute Q.r and Q.r.res
  arma::vec out_Q = SKAT_META_Optimal_Get_Q(Score, r_all);
  arma::mat Q_res;
  arma::mat  out_Q_res;
  if (!Score_Resampling.is_empty()) {
    out_Q_res = SKAT_META_Optimal_Get_Q_Res(Score_Resampling, r_all);
  }
  arma::mat Q_all = join_rows(out_Q.t(), out_Q_res);

  // Compute P-values
  arma::mat Phihalf = Phi / 2;

  //r_all.print("r_all herea");
  //Phihalf.print("Phihalf");
  Rcpp::List out = SKAT_META_Optimal_Get_Pvalue(Q_all, Phihalf, r_all, method, isFast);
  //r_all.print("r_all hereb");


  //std::cout << "after SKAT_META_Optimal_Get_Pvalue" << std::endl;

  // Initialize result parameters
  Rcpp::List param;
  arma::mat pvaleachmat = out["p.val.each"];
  //pvaleachmat.print("pvaleachmat");
  arma::mat pvaleachmatt = pvaleachmat.t();
  arma::vec pvaleachvec = pvaleachmatt.col(0);
  //arma::vec pvaleachvec = pvaleachmat(0, arma::span::all);
  //pvaleachvec.print("pvaleachvec");
  param["p.val.each"] = pvaleachvec;
  param["q.val.each"] = Q_all.row(0);
  param["rho"] = r_all;
  param["minp"] = arma::min(pvaleachvec);

  // Adjust rho values if necessary
  double minp =  param["minp"];
  arma::vec p_val_each = param["p.val.each"];
  arma::vec rho = param["rho"];
  arma::uvec id_temp = arma::find(p_val_each == minp);
  arma::uvec id_temp1 = arma::find(rho >= 0.999);  // treat rho >= 0.999 as 1
  arma::vec rhosubvec;
  //id_temp1.print("id_temp1");
  if (id_temp1.n_elem > 0) {
    rho(id_temp1).fill(1);
    param["rho"] = rho;
  }
  //rho.print("rho");

  rhosubvec = rho(id_temp);
  param["rho_est"] = rhosubvec;

  // Extract p-value
  arma::vec pvalvec=out["p.value"];
  double p_value = pvalvec[0];
  arma::vec p_value_resampling;
  if (!Q_res.is_empty()) {
    p_value_resampling = pvalvec.subvec(1, pvalvec.n_elem - 1);
  }

  return Rcpp::List::create(
    Rcpp::Named("p.value") = p_value,
    Rcpp::Named("param") = param,
    Rcpp::Named("p.value.resampling") = p_value_resampling
  );
}

//Rcpp::List Met_SKAT_Get_Pvalue( arma::vec &Score,  arma::mat &Phi,  arma::vec &r_corr, std::string &method,  arma::mat &Score_Resampling, bool isFast) {

// [[Rcpp::export]]
Rcpp::List Met_SKAT_Get_Pvalue( arma::vec &Score,  arma::mat &Phi,  arma::vec &r_corr, std::string &method, bool isFast) {

  int p_m = Phi.n_rows;
  arma::mat Q_res;
 
  //std::cout << "Met_SKAT_Get_Pvalue here" << std::endl;
  // If Phi is a zero matrix
  if (arma::accu(arma::abs(Phi)) == 0) {
    Rcpp::warning("No polymorphic SNPs!");
    return Rcpp::List::create(Rcpp::Named("p.value") = 1, 
                              Rcpp::Named("p.value.resampling") = R_NilValue,
                              Rcpp::Named("pval.zero.msg") = R_NilValue);
  }

  // Resampling handling
  arma::mat Score_Res;
  //if(Score_Resampling.n_rows != 0){
  //  Score_Res = Score_Resampling.t();
  //}

  //Phi.print("Phi");

  if(Phi.n_rows <= 1){
     r_corr.set_size(1);
     r_corr(0) = 0;
  }else{
     if(Phi.n_cols <= 10){
	//Phi.print("Phi"); 

        arma::mat Q, R;
        arma::qr(Q, R, Phi);
        int rank = arma::rank(R);
        if(rank <= 1){
		r_corr.set_size(1);
      	        r_corr(0) = 0;
        }
	//std::cout << "rank " << rank << std::endl;
     }
  }
 
  //r_corr.print("r_corr");

  arma::mat Score_Resampling;
  // Check if r_corr has more than 1 value
  if (r_corr.n_elem > 1) {
    Rcpp::List re = SKAT_META_Optimal(Score, Phi, r_corr, method, Score_Resampling, isFast);
    return re;
  }

  arma::vec Q;
  //Score_Resampling.print("Score_Resampling");
  //std::cout << "Score_Resampling.is_empty() " << Score_Resampling.is_empty() << std::endl; 


  // If r_corr == 0
  if (r_corr(0) == 0) {
    Q.set_size(1);
    Q(0)= arma::accu(arma::square(Score)) / 2;

    if (!Score_Resampling.is_empty()) {
      Q_res = arma::sum(arma::square(Score_Res), 1) / 2;
    }
    //return Rcpp::List::create(Rcpp::Named("p.value") = sum(Q), 
    //                          Rcpp::Named("p.value.resampling") = Q_res, 
    //                          Rcpp::Named("pval.zero.msg") = Rcpp::CharacterVector::create("r_corr is 0"));
  } else if (r_corr(0) == 1) {
     Q = SKAT_META_Optimal_Get_Q(Score, r_corr);
    if (!Score_Resampling.is_empty()) {
      arma::mat Q_res = SKAT_META_Optimal_Get_Q_Res(Score_Res, r_corr);
    }
    double a = arma::accu(Phi);
    arma::mat amat(1,1);
    amat(0,0) = a;
    Rcpp::List re = Get_Liu_PVal(Q, amat, Q_res);
    return re;
  }else {
  // General case for 0 < r_corr < 1
    Q = SKAT_META_Optimal_Get_Q(Score, r_corr);
    if (!Score_Resampling.is_empty()) {
      arma::mat Q_res = SKAT_META_Optimal_Get_Q_Res(Score_Res, r_corr);
    }

    // Compute the correlation matrix and its Cholesky decomposition
    arma::vec pmone;
    pmone.ones(p_m);
    bool success;
    arma::mat L, Phi_rho;   
    arma::mat R_M = diagmat(arma::ones<arma::vec>(p_m) * (1 - r_corr)) + arma::ones<arma::mat>(p_m, p_m) * r_corr;
     success = arma::chol(L, R_M, "lower");
     Phi_rho = L * (Phi * L.t());
     Phi.clear();
     Phi = Phi_rho;
  }
    //std::cout << "Get_Davies_PVal here" << std::endl;
    //r_corr.print("r_corr");
    Rcpp::List re = Get_Davies_PVal(Q, Phi, Q_res, isFast);
    if (r_corr.n_elem == 1) {
      re["Q"] = Q;
    }
    return re;
  
}


// Helper function to compute parameters for the null approximation
// [[Rcpp::export]]
Rcpp::List Get_Liu_Params( arma::vec& c1) {
    // Extract moments from c1
    double muQ = c1(0);
    double sigmaQ = std::sqrt(2.0 * c1(1));
    double s1 = c1(2) / std::pow(c1(1), 1.5);
    double s2 = c1(3) / std::pow(c1(1), 2.0);

    // Compute beta coefficients
    double beta1 = std::sqrt(8.0) * s1;
    double beta2 = 12.0 * s2;

    // Parameters for type determination
    double a, d, l;
    int type1 = 0;

    if (std::pow(s1, 2) > s2) {
        a = 1.0 / (s1 - std::sqrt(std::pow(s1, 2) - s2));
        d = s1 * std::pow(a, 3) - std::pow(a, 2);
        l = std::pow(a, 2) - 2.0 * d;
    } else {
        type1 = 1;
        a = 1.0 / s1;
        d = 0.0;
        l = 1.0 / std::pow(s1, 2);
    }

    // Compute final parameters
    double muX = l + d;
    double sigmaX = std::sqrt(2.0) * a;

    // Return results as a list
    return Rcpp::List::create(
        Rcpp::Named("l") = l,
        Rcpp::Named("d") = d,
        Rcpp::Named("muQ") = muQ,
        Rcpp::Named("muX") = muX,
        Rcpp::Named("sigmaQ") = sigmaQ,
        Rcpp::Named("sigmaX") = sigmaX
    );
}

// [[Rcpp::export]]
Rcpp::List Get_Liu_PVal( arma::vec& Q,  arma::mat& W,  arma::mat& Q_resampling) {
    // Convert inputs to appropriate Armadillo structures
    arma::vec Q_all = arma::join_cols(Q, Q_resampling);

    // Calculate intermediate matrices
    arma::mat A1 = W / 2.0;
    arma::mat A2 = A1 * A1;

    // Compute c1 vector
    arma::vec c1(4);
    c1(0) = arma::trace(A1);
    c1(1) = arma::trace(A2);
    c1(2) = arma::accu(A1 % A2);  // Element-wise product and sum
    c1(3) = arma::accu(A2 % A2);

    // Get Liu parameters
    Rcpp::List param = Get_Liu_Params(c1);

    // Extract parameters
    double muQ = param["muQ"];
    double sigmaQ = param["sigmaQ"];
    double muX = param["muX"];
    double sigmaX = param["sigmaX"];
    double l = param["l"];
    double d = param["d"];

    // Normalize Q values
    arma::vec Q_Norm = (Q_all - muQ) / sigmaQ;
    arma::vec Q_Norm1 = Q_Norm * sigmaX + muX;

    // Compute p-values
    arma::vec p_values(Q_Norm1.n_elem);
    boost::math::non_central_chi_squared dist(l, d);


    for (size_t i = 0; i < Q_Norm1.n_elem; i++) {
        p_values(i) = 1 - boost::math::cdf(dist, Q_Norm1(i)); // pchisq with lower.tail=FALSE
    }

    // Extract primary and resampling p-values
    double p_value = p_values(0);
    arma::vec p_value_resampling;
    if (Q_resampling.n_elem > 0) {
        p_value_resampling = p_values.subvec(1, Q_Norm1.n_elem - 1);
    }

    // Return results as a list
    return Rcpp::List::create(
        Rcpp::Named("p.value") = p_value,
        Rcpp::Named("param") = param,
        Rcpp::Named("p.value.resampling") = p_value_resampling // Wrap in Rcpp object
    );
}

// [[Rcpp::export]]
arma::mat forceSymmetric(const arma::mat& K) {
    return 0.5 * (K + K.t());
}





// [[Rcpp::export]]
void get_SKAT_pvalue_cpp(arma::vec& Score, 
                               arma::mat& Phi, 
                               arma::vec& r_corr,
			       double& Pvalue_SKATO,
			       double& Pvalue_Burden, 
			       double& Pvalue_SKAT,
			       double& BETA_Burden,
			       double& SE_Burden,
			       int& error_code
			       ) {

  std::string method = "optimal.adj";
  bool isFast = false;
  //r_corr.print("r_corr here");
  Rcpp::List out_SKAT_List = Met_SKAT_Get_Pvalue(Score, Phi, r_corr, method, isFast);  
  BETA_Burden = arma::accu(Score) / arma::trace(Phi);
  error_code = 0;
  double p_value = out_SKAT_List["p.value"];
  

if(r_corr.n_elem > 1){
  //std::cout << "p_value " << p_value << std::endl;
  //arma::vec Pvalue(3, arma::fill::zeros);
  Rcpp::RObject paramElement = out_SKAT_List["param"];
  Rcpp::List param = Rcpp::as<Rcpp::List>(paramElement); 
  arma::vec rho = param["rho"];
  arma::vec p_val_each = param["p.val.each"];
  Pvalue_SKATO = p_value;


  // Check conditions similar to the R logic
  if (p_val_each.is_empty() || rho.is_empty()) {
    error_code = 2;
    BETA_Burden = NA_REAL;
  } else if (!arma::any(rho == 0) || !arma::any(rho == 1)) {
    SE_Burden = std::abs(BETA_Burden / R::qnorm(p_value/2, 0, 1, true, false));
    Pvalue_SKAT = p_value;
    Pvalue_Burden = p_value;
  } else {
    int pos00 = arma::as_scalar(arma::find(rho == 0, 1, "first"));
    int pos01 = arma::as_scalar(arma::find(rho == 1, 1, "first"));
    Pvalue_SKAT = p_val_each[pos00];
    Pvalue_Burden = p_val_each[pos01];
    SE_Burden = std::abs(BETA_Burden / R::qnorm(p_val_each[pos01] / 2, 0, 1, true, false));
  }
}else if(r_corr.n_elem == 1){
    if(r_corr[0] == 0){	
	Pvalue_SKAT = p_value;
    }else if(r_corr[0] == 1){
	Pvalue_Burden = p_value;
    }
  }
}


// [[Rcpp::export]]
double get_jointScore_pvalue(arma::vec& Score, arma::mat& Phi) {
    arma::mat Phi_inv = arma::inv(Phi); // Inverse of Phi
    double Teststat = arma::as_scalar(Score.t() * Phi_inv * Score);
    
    // Degrees of freedom
    int df = Score.n_elem;
    
    // Debugging print statements (optional, remove in production)
    //Rcpp::Rcout << "Score:" << std::endl << Score << std::endl;
    //Rcpp::Rcout << "Phi:" << std::endl << Phi << std::endl;
    //Rcpp::Rcout << "Teststat: " << Teststat << std::endl;
    //Rcpp::Rcout << "df: " << df << std::endl;
    
    // Create central chi-squared distribution
    boost::math::chi_squared chi2(df);

    // Compute the p-value
    double p_value;
    if (Teststat > 0) {
        p_value = 1 - boost::math::cdf(chi2, Teststat); // 1 - CDF for upper-tail probability
    } else {
        p_value = 1.0;
    }    
    
    return p_value;
}

// [[Rcpp::export]]
void SPA_ER_kernel_related_Phiadj_fast_new_cpp(arma::vec& p_new,
                                                 arma::vec& Score,
                                                 arma::mat& Phi,
                                                 double p_value_burden,
                                                 std::string regionTestType,
						 arma::vec& scaleFactor
						 ) {
    int p_m = Score.n_elem;
    arma::vec zscore_all_0 = Score;
    arma::vec zscore_all_1(p_m, arma::fill::zeros);

    arma::vec VarS_org = (regionTestType != "BURDEN") ? Phi.diag() : arma::vec(Phi);

    arma::vec stat_qtemp = arma::square(Score) / VarS_org;

    arma::uvec idx_0 = arma::find(VarS_org > 0);
    arma::uvec idx_p0 = arma::find_finite(p_new);

    arma::vec VarS = arma::square(zscore_all_0) / 500.0;

    
    if (!idx_p0.is_empty()) {
        double df = 1.0;
    	for (size_t i = 0; i < idx_p0.n_elem; i++) {
		unsigned int idx_p0_i = idx_p0(i);
		double p_new_i = p_new(idx_p0_i);
		VarS(idx_p0_i) = std::pow(zscore_all_0(idx_p0_i),2) / qchisq_log(p_new_i, df);
    	}
    }
    arma::uvec vars_inf = arma::find(VarS == arma::datum::inf);
    if (regionTestType != "BURDEN") {
        if (!vars_inf.is_empty()) {
            VarS(vars_inf).fill(0);
            zscore_all_1(vars_inf).fill(0);
            Phi.rows(vars_inf).zeros();
            Phi.cols(vars_inf).zeros();
        }
    } else {
        if (!vars_inf.is_empty()) {
            VarS(vars_inf).fill(0);
            zscore_all_1(vars_inf).fill(0);
            Phi(vars_inf).zeros();
        }
    }

    scaleFactor = arma::sqrt(VarS / VarS_org);
    zscore_all_1 = zscore_all_0;

    arma::vec VarStoorg = VarS / VarS_org;

    arma::mat G2_adj_n;
    if (regionTestType != "BURDEN") {
        arma::mat PhiVarS = (Phi % arma::sqrt(VarStoorg)).t();
        G2_adj_n = (PhiVarS % arma::sqrt(VarStoorg)).t();
    } else {
        G2_adj_n = Phi % VarStoorg;
    }

    double VarQ = arma::accu(G2_adj_n);
    double Q_b = std::pow(arma::accu(zscore_all_1), 2);

    double VarQ_2 = Q_b / qchisq_log(p_value_burden, 1);
    double r = (VarQ_2 == 0) ? 1.0 : VarQ / VarQ_2;
    r = std::min(r, 1.0);

    //arma::mat Phi_ccadj;
    Phi = G2_adj_n / r;
    scaleFactor /= std::sqrt(r);

}


// Helper function for qchisq with log.p = TRUE
double qchisq_log(double pval, double df) {
    boost::math::chi_squared dist(df);
    double qval;
    qval = boost::math::quantile(boost::math::complement(dist, std::exp(pval)));
}


// [[Rcpp::export]]
void get_newPhi_scaleFactor_cpp(double q_sum,
                                  arma::vec& mu_a,
                                  arma::vec& g_sum,
                                  arma::vec& p_new,
                                  arma::vec& Score,
                                  arma::mat& Phi,
                                  std::string regionTestType,
                                  arma::vec& scaleFactor) {
    double m1 = arma::accu(mu_a % g_sum);
    double var1 = arma::accu(mu_a % (1.0 - mu_a) % arma::square(g_sum));

    double qinv = -std::copysign(1.0, q_sum - m1) * std::abs(q_sum - m1) + m1;
    //        pval_noadj <- pchisq((q.sum - m1)^2/var1, lower.tail = FALSE, df = 1, log.p = TRUE)
    boost::math::chi_squared chi2(1.0);
    double qval = std::pow(q_sum - m1, 2) / var1;
    double pval_noadj = 1.0 - boost::math::cdf(chi2, qval);
    pval_noadj = std::log(pval_noadj);

    bool isSPAConverge = true;

    double p_value_burden;
    if (std::abs(q_sum - m1) / std::sqrt(var1) < 2) {
        p_value_burden = pval_noadj;
    } else {
        double epsilon = std::numeric_limits<double>::epsilon();
	double epsilon_0_25 = std::pow(epsilon, 0.25);
	p_value_burden = SPA_pval(mu_a, g_sum, q_sum, qinv, pval_noadj, epsilon_0_25, true, "binary", isSPAConverge); 
    }

    SPA_ER_kernel_related_Phiadj_fast_new_cpp(p_new, Score, Phi, p_value_burden, regionTestType, scaleFactor); 	
}

/*
//https://github.com/yixuan/RcppNumerical
//double integrate_SKAT_Optimal_Liu_v2(arma::vec &pmin_q, Rcpp::List &param_m, arma::vec &r_all, double lower, double upper, int subdivisions, double abs_tol){
// [[Rcpp::export]]
Rcpp::List integrate_SKATliu(arma::vec &pmin_q, Rcpp::List &param_m, arma::vec &r_all, double lower, double upper, int subdivisions, double abs_tol)
{

    SKATLiu f(pmin_q, param_m, r_all);

    // Variables to store integration error estimate and code
    double err_est;
    int err_code;

    // Perform numerical integration using Numer::integrate
    const double res = Numer::integrate(f, lower, upper, err_est, err_code, subdivisions, abs_tol);

    // Return the results as an Rcpp list
    return Rcpp::List::create(
        Rcpp::Named("res") = res,
        Rcpp::Named("error_estimate") = err_est,
        Rcpp::Named("error_code") = err_code
    );
}
*/
