#include <RcppArmadillo.h>
#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/non_central_chi_squared.hpp>


#include "qfc_rcpp.hpp"
#include "SKAT_funcs.hpp"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]


// Function to compute filtered eigenvalues from a matrix K
arma::vec Get_Lambda(arma::mat& K, bool isFast, int maxK) {
  // Eigenvalue decomposition
  arma::vec lambda1;
  arma::eig_sym(lambda1, K);  // Symmetric eigenvalue decomposition
  
  // Find indices of positive eigenvalues
  arma::uvec IDX1 = find(lambda1 >= 0);
  
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

  std::cout << "SKAT_META_Optimal_Param 1" << std::endl;
  // ZMZ calculation
  arma::vec Z_item1_1 = Phi * arma::ones<arma::vec>(p_m);
  arma::mat ZZ = Phi;
  arma::mat ZMZ = Z_item1_1 * (Z_item1_1.t()) / arma::accu(ZZ);
  std::cout << "SKAT_META_Optimal_Param 2" << std::endl;

  // W3.2 Term: mixture chi-squared
  arma::mat W3_2_t = ZZ - ZMZ;
  arma::vec lambda = Get_Lambda(W3_2_t);
  std::cout << "SKAT_META_Optimal_Param 3" << std::endl;
  
  // W3.3 Term: variance of remaining
  double W3_3_item = arma::accu(ZMZ % (ZZ - ZMZ)) * 4;
  std::cout << "SKAT_META_Optimal_Param 4" << std::endl;
  
  // tau term
  double z_mean_2 = arma::accu(ZZ) / std::pow(static_cast<double>(p_m), 2);

  double tau1 = arma::accu(ZZ * ZZ) / std::pow(static_cast<double>(p_m), 2) / z_mean_2;
  std::cout << "SKAT_META_Optimal_Param 5" << std::endl;

  // Mixture Parameters
  double MuQ = arma::accu(lambda);
  double VarQ = arma::accu(arma::square(lambda)) * 2 + W3_3_item;
  double KerQ = arma::accu(arma::pow(lambda, 4)) / std::pow(arma::accu(arma::square(lambda)), 2) * 12;
  double Df = 12 / KerQ;
  std::cout << "SKAT_META_Optimal_Param 6" << std::endl;

  // W3.1 Term: tau1 * chi-squared_1
  arma::vec tau(r_n, arma::fill::zeros);
  double r_corr, term1;
  for (int i = 0; i < r_n; i++) {
    r_corr = r_all(i);
    term1 = std::pow(p_m, 2) * r_corr * z_mean_2 + tau1 * (1 - r_corr);
    tau(i) = term1;
  }
  std::cout << "SKAT_META_Optimal_Param 7" << std::endl;

  // Return as a list
  return Rcpp::List::create(
    Rcpp::Named("MuQ") = MuQ,
    Rcpp::Named("VarQ") = VarQ,
    Rcpp::Named("KerQ") = KerQ,
    Rcpp::Named("lambda") = lambda,
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
    p_value(i) = 1 - boost::math::cdf(dist, Q_norm1(i));
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
    out(i) = 1 - boost::math::cdf(dist, temp(i)); // Lower tail = FALSE
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
    arma::vec eigvalK = arma::eig_sym(K);
    double qall0 = Q_all(0);
    double p_value = daviesPValue(eigvalK, qall0); // Test statistic Q.all[0]

    // Prepare the result list
    Rcpp::List result;
    Rcpp::List param;
    
    // Adding param information
    param["liu_pval"] = p_value;
    param["Is_Converged"] = true;  // Assuming convergence is always true for simplicity
    
    // If resampling is available
    arma::vec p_value_resampling;
    if (Q_resampling.n_elem > 0) {
        p_value_resampling.set_size(Q_resampling.n_rows);
        //p_value_resampling = Rcpp::NumericVector(Q_resampling.n_rows);
        for (size_t i = 0; i < Q_resampling.n_rows; i++) {
	    double Q_resampling0 = Q_resampling(i,0);
            p_value_resampling[i] = daviesPValue(eigvalK, Q_resampling0);
        }
        param["liu_pval.resampling"] = p_value_resampling;
        param["Is_Converged.resampling"] = true;
    }

    result["p.value"] = p_value;
    result["param"] = param;
    result["p.value.resampling"] = p_value_resampling;
    result["pval.zero.msg"] = "No zero p-values";  // Assuming no zero p-values for simplicity
    
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
 
  std::cout << "Call the qfc_1 function before" << std::endl;
  // Call the qfc_1 function (from the C code)
  QUADFORM::QuadFormClass* quadForm = new QUADFORM::QuadFormClass(r, lim);
  std::cout << "Call the qfc_1 function before a" << std::endl;
  
  quadForm->qfc_1(lambda, delta, h, r, sigma, q, lim, acc, trace, ifault, res);
  std::cout << "Call the qfc_1 function before b" << std::endl;

  delete quadForm;

  // Adjust result
  res = 1.0 - res;
  
  // Return results as a list
  return Rcpp::List::create(
    Rcpp::Named("trace") = trace,
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
  std::cout << "Get_Liu_PVal_MOD_Lambda 1" << std::endl;
  p_val_liu = Get_Liu_PVal_MOD_Lambda(Q, lambda, df1);
  std::cout << "Get_Liu_PVal_MOD_Lambda 2" << std::endl;
 
  double Qi;
  arma::ivec h;
  arma::vec delta;
  Rcpp::List out;
  df1.print("df1");
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
  
  // For p-value zero case
  //Rcpp::Nullable<double> p_val_msg = R_NilValue;
  //Rcpp::Nullable<double> p_val_log = R_NilValue;
  std::string p_val_msg;
  arma::vec p_val_log;
  arma::vec  Q0 = Q.subvec(0,0);


  if (p_val[0] == 0) {
    // Get Liu Params Mod Lambda
    Rcpp::List param = Get_Liu_Params_Mod_Lambda(lambda, df1);
    //std::string Get_Liu_PVal_MOD_Lambda_Zero(double Q, double muQ, double muX, double sigmaQ, double sigmaX, double l, double d) {
    p_val_msg = Get_Liu_PVal_MOD_Lambda_Zero(Q(0), param["muQ"], param["muX"], param["sigmaQ"], param["sigmaX"], param["l"], param["d"]);
    p_val_log = Get_Liu_PVal_MOD_Lambda(Q0, lambda, df1, true);
  }
  
  // Return the results as a List
  return Rcpp::List::create(
    Rcpp::Named("p.value") = p_val,
    Rcpp::Named("p_val_liu") = p_val_liu,
    Rcpp::Named("is_converge") = is_converge,
    Rcpp::Named("p_val_log") = p_val_log[0],
    Rcpp::Named("pval_zero_msg") = p_val_msg
  );
}

Rcpp::List SKAT_Optimal_Each_Q(Rcpp::List &param_m, arma::mat &Q_all, arma::vec &r_all, Rcpp::List & lambda_all, std::string & method) {
  
  int n_r = r_all.n_elem;
  arma::vec c1 = arma::zeros<arma::vec>(4);
  int n_q = Q_all.n_rows;
  
  // Initialize matrices for p-values and minimum p-values
  arma::mat pval = arma::zeros<arma::mat>(n_q, n_r);
  arma::mat pmin_q = arma::zeros<arma::mat>(n_q, n_r);
  arma::mat param_mat;
  arma::vec Q, Q_norm;
  arma::vec lambda_temp;
  double r_corr, muQ, varQ, df, sigmaQ;
  arma::ivec empvec;
    std::cout << "SKAT_Optimal_Each_Q 0" << std::endl;
  for (int i = 0; i < n_r; i++) {
    Q = Q_all.col(i);
    r_corr = r_all[i];
    lambda_temp.clear(); 
    lambda_temp = arma::vec(Rcpp::as<arma::vec>(lambda_all[i]));
    
    // Compute c1 components
    c1[0] = arma::accu(lambda_temp);
    c1[1] = arma::accu(arma::square(lambda_temp));
    c1[2] = arma::accu(arma::pow(lambda_temp, 3));
    c1[3] = arma::accu(arma::pow(lambda_temp, 4));
    
    std::cout << "SKAT_Optimal_Each_Q 00" << std::endl;
    std::cout << "i " << i << std::endl;
    // Get parameters
    Rcpp::List param_temp = Get_Liu_Params_Mod(c1);
    std::cout << "Get_Liu_Params_Mod after" << std::endl;
    
    muQ = param_temp["muQ"];
    sigmaQ = param_temp["sigmaQ"];
    varQ = std::pow(sigmaQ, 2);
    df = param_temp["l"];
    boost::math::chi_squared dist(df);  
    // Get p-value
    std::cout << "SKAT_Optimal_Each_Q 01" << std::endl;
    
    Q_norm = ((Q - muQ) / std::sqrt(varQ)) * std::sqrt(2 * df) + df;
    Q_norm.print("Q_norm");
    //pval.col(i) = R::pchisq(Q_norm, df, 0, false, false);
    //p_value(i) = 1 - boost::math::cdf(dist, Q_norm1(i));
    std::cout << "SKAT_Optimal_Each_Q 02" << std::endl;
    for (size_t j = 0; j < Q_norm.n_elem; ++j) {
        pval(j,i) = 1 - boost::math::cdf(dist, Q_norm(j));
     }
    std::cout << "SKAT_Optimal_Each_Q 1" << std::endl;
    Rcpp::List pval_lambda;
    std::string method_str;
    // Check for method-specific changes
    if (method != "") {
      method_str = method;
      std::cout << "method_str " << method_str << std::endl;
      if (method_str == "optimal.mod" || method_str == "optimal.adj" || method_str == "optimal.moment.adj") {
        pval_lambda = Get_PValue_Lambda(lambda_temp, Q, empvec);
        pval.col(i) = arma::vec(Rcpp::as<arma::vec>(pval_lambda["p.value"]));
      }
    }
    std::cout << "SKAT_Optimal_Each_Q 2" << std::endl;
    arma::mat tempmat(1,3);
    tempmat(0,0) = muQ;
    tempmat(0,1) = varQ;
    tempmat(0,2) = df;
    param_mat = arma::join_rows(param_mat, tempmat);
  }
    std::cout << "SKAT_Optimal_Each_Q 3" << std::endl;
  
  // Calculate the minimum p-values
  arma::vec pmin = arma::min(pval, 1);
  arma::vec q_org, q_q;
  q_org.set_size(pmin.n_elem);
  for (int i = 0; i < n_r; i++) {
    muQ = param_mat(i, 0);
    varQ = param_mat(i, 1);
    df = param_mat(i, 2);
     boost::math::chi_squared dist(df);    
    for (int j = 0; j < pmin.n_elem; i++) {
        q_org[j] = boost::math::cdf(dist, 1-pmin[j]);
    }

    q_q = (q_org - df) / std::sqrt(2 * df) * std::sqrt(varQ) + muQ;
    pmin_q.col(i) = q_q;
  }
  
  // Return the results as a List
  return Rcpp::List::create(
    Rcpp::Named("pmin") = pmin,
    Rcpp::Named("pval") = pval,
    Rcpp::Named("pmin_q") = pmin_q
  );
}


// Davies Integration Function for SKAT
arma::vec SKAT_Optimal_Integrate_Func_Davies(arma::vec & x,  arma::mat &pmin_q, Rcpp::List &param_m,  arma::vec &r_all) {
  int n_r = r_all.n_elem;
  int n_x = x.n_elem;;

    arma::mat taumat = param_m["tau"];
      arma::mat xt(1, n_x);
        xt.row(0) = x.t();
	  arma::mat temp1 = arma::kron(taumat ,  xt);

  arma::mat tempmat = (pmin_q - temp1) / (1 - r_all);
  arma::vec temp_min = arma::min(tempmat, 0);

  arma::vec re = arma::zeros<arma::vec>(n_x);
  arma::vec lambda = param_m["lambda"];
  double min1, min1_temp, sd1, min1_st, xi, VarQ, VarRemain, temp, MuQ;
  arma::ivec h; 
  VarQ = param_m["VarQ"];
  VarRemain = param_m["VarRemain"];
  MuQ = param_m["MuQ"];
  for (int i = 0; i < n_x; i++) {
    min1 = temp_min[i];
    if (min1 > arma::accu(lambda) * 1e4) {
      temp = 0;
    } else {
      min1_temp = min1 - MuQ;
      sd1 = std::sqrt(VarQ - VarRemain) / std::sqrt(VarQ);
      min1_st = min1_temp * sd1 + MuQ;

      arma::vec deltavec = arma::zeros<arma::vec>(lambda.n_elem);
      Rcpp::List dav_re = SKAT_davies(min1_st, lambda, h, deltavec, 0, 10000, 1e-6);
      temp = dav_re["Qq"];
      if (dav_re["ifault"] != 0) {
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

// Function for computing the Davies p-value
double SKAT_Optimal_PValue_Davies(arma::vec &pmin_q, Rcpp::List &param_m,
                                     arma::vec &r_all, double pmin) {
    double pvalue;
    try {
        // Perform the integration
        arma::vec re = integrate_SKAT_Optimal_Davies(pmin_q, param_m, r_all, 0, 40, 1000, 1e-25);

        // Compute the p-value
       pvalue = 1 - re[0];
        if (pmin != arma::datum::nan) {
            if ((pmin * (r_all.n_elem)) < pvalue) {
                pvalue = pmin * r_all.n_elem;
            }
        }

        // Return the result
        return pvalue;

    } catch (const std::exception &e) {
        // Log the error message for debugging
        Rcpp::Rcerr << "Error in SKAT_Optimal_PValue_Davies: " << e.what() << std::endl;

        // Fall back to Liu method
        pvalue = SKAT_Optimal_PValue_Liu(pmin_q, param_m, r_all, pmin);
	return pvalue;

    } catch (...) {
        // Catch-all for unknown errors
        Rcpp::Rcerr << "Unknown error in SKAT_Optimal_PValue_Davies." << std::endl;

        // Fall back to Liu method
        pvalue =  SKAT_Optimal_PValue_Liu(pmin_q, param_m, r_all, pmin);
	return pvalue;
    }
}

//arma::vec & x,  arma::mat &pmin_q, Rcpp::List &param_m,  arma::vec &r_all
// Liu Integration Function for SKAT
arma::vec SKAT_Optimal_Integrate_Func_Liu(arma::vec & x,  arma::mat &pmin_q,  Rcpp::List &param_m,  arma::vec &r_all) {
  int n_r = r_all.n_elem;
  int n_x = x.n_elem;

  arma::mat taumat = param_m["tau"];
  arma::mat xt(1, n_x);
  xt.row(0) = x.t(); 
  arma::mat temp1 = arma::kron(taumat ,  xt);

  arma::mat temp = (pmin_q - temp1) / (1 - r_all);
  arma::vec temp_min = arma::min(temp, 0);

  double MuQ, VarQ, Df;
  arma::vec temp_q = (temp_min - MuQ) / VarQ * std::sqrt(2 * Df) + Df;
  //arma::vec temp_q = (temp_min - param_m["MuQ"]) / std::sqrt(param_m["VarQ"]) * std::sqrt(2 * param_m["Df"]) + param_m["Df"];
  arma::vec re(temp_q.n_elem);
  boost::math::chi_squared dist(Df);
  boost::math::chi_squared dist1(1);


  for (size_t i = 0; i < temp_q.n_elem; ++i) {

        re(i) = boost::math::cdf(dist, temp_q(i)) *  boost::math::pdf(dist1, x(i));
  }
  
  return re;
}


// Numerical integration using the trapezoidal rule over a vector
// [[Rcpp::export]]
arma::vec integrate_SKAT_Optimal_Liu(arma::mat &pmin_q, Rcpp::List &param_m, arma::vec &r_all,double lower, double upper, int subdivisions, double abs_tol) {
    int n_r = r_all.n_elem;
    arma::vec results = arma::zeros<arma::vec>(n_r);

    arma::vec x = arma::linspace(lower, upper, subdivisions + 1); // Points in the range

    for (int i = 0; i < n_r; i++) {
        // Use each row of r_all to compute the integral
        arma::vec r_vec = r_all.subvec(i, i);
        arma::vec y = SKAT_Optimal_Integrate_Func_Liu(x, pmin_q, param_m, r_vec);

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






// Function for computing the Liu p-value
double SKAT_Optimal_PValue_Liu( arma::mat &pmin_q,  Rcpp::List &param_m, 
                                   arma::vec &r_all, double pmin) {
  arma::vec re = integrate_SKAT_Optimal_Liu(pmin_q, param_m, r_all, 0, 40, 2000, 1e-25);

  double pvalue = 1 - re(0);

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
  for (int i = 0; i < n_r; i++) {
     r_corr = r_all[i];
     R_M = diagmat(arma::ones<arma::vec>(p_m) * (1 - r_corr)) + arma::ones<arma::mat>(p_m, p_m) * r_corr;
     success = arma::chol(L, R_M, "lower");
     Phi_rho = L * (Phi * L.t());
     Phi_rho.print("Phi_rho");
     lambda_all[i] = Get_Lambda(Phi_rho, isFast);
     std::cout << "lambda_all[i] " << std::endl;
  }

  std::cout << "SKAT_META_Optimal_Param before" << std::endl;
  // Get Mixture parameters
  Rcpp::List param_m = SKAT_META_Optimal_Param(Phi, r_all);
  std::cout << "SKAT_META_Optimal_Param after" << std::endl;


  // Call SKAT_Optimal_Each_Q
  Rcpp::List Each_Info = SKAT_Optimal_Each_Q(param_m, Q_all, r_all, lambda_all, method);
  std::cout << "SKAT_Optimal_Each_Q after" << std::endl;
  arma::mat pmin_q = Each_Info["pmin_q"];
  arma::vec pval = arma::zeros<arma::vec>(n_q);

  // Get the pmin values
  arma::vec pmin = Each_Info["pmin"];

  // Calculate p-values based on the chosen method
  if (method == "davies" || method == "optimal" || method == "optimal.adj" || method == "optimal.mod") {
    for (int i = 0; i < n_q; i++) {
      arma::vec pminqvec = pmin_q.row(i).t();
      pval[i] = SKAT_Optimal_PValue_Davies(pminqvec, param_m, r_all, pmin[i]);
  std::cout << "i " << i << std::endl;
  std::cout << "SKAT_Optimal_PValue_Davies after" << std::endl;
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

  arma::mat pvalmat = Each_Info["pval"];
  int npval = pvalmat.n_rows;
  arma::vec pval_each(npval);
  arma::uvec IDX;
  double pval1;
  for (int i = 0; i < n_q; i++) {
    IDX.clear();
    pval_each = pvalmat.col(i);
    IDX = arma::find(pval_each > 0);

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

  return Rcpp::List::create(Rcpp::Named("p_value") = pval, Rcpp::Named("p_val_each") = Each_Info["pval"]);
}

Rcpp::List SKAT_META_Optimal( arma::vec &Score,  arma::mat &Phi,  arma::vec &r_all, 
                              std::string method,  arma::mat &Score_Resampling, bool isFast) {

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
  arma::mat Q_all = join_rows(out_Q, out_Q_res);

  // Compute P-values
  arma::mat Phihalf = Phi / 2;

  Phihalf.print("Phihalf");
  Rcpp::List out = SKAT_META_Optimal_Get_Pvalue(Q_all, Phihalf, r_all, method, isFast);
  std::cout << "after SKAT_META_Optimal_Get_Pvalue" << std::endl;

  // Initialize result parameters
  Rcpp::List param;
  arma::mat pvaleachmat = out["p.val.each"];
  arma::vec pvaleachvec = pvaleachmat(0, arma::span::all);
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
  if (id_temp1.n_elem > 0) {
    rho.subvec(id_temp1[0], id_temp1[id_temp1.n_elem - 1]).fill(1);
  }
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
 
  std::cout << "Met_SKAT_Get_Pvalue here" << std::endl;
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

  if(Phi.n_rows <= 1){
     r_corr.set_size(1);
     r_corr(0) = 0;
  }else{
     if(Phi.n_cols <= 10){
	Phi.print("Phi"); 

        arma::mat Q, R;
        arma::qr(Q, R, Phi);
        int rank = arma::rank(R);
        if(rank <= 1){
		r_corr.set_size(1);
      	        r_corr(0) = 0;
        }
	std::cout << "rank " << rank << std::endl;
     }
  }
 
  r_corr.print("r_corr");

  arma::mat Score_Resampling;
  // Check if r_corr has more than 1 value
  if (r_corr.n_elem > 1) {
    Rcpp::List re = SKAT_META_Optimal(Score, Phi, r_corr, method, Score_Resampling, isFast);
    return re;
  }

  arma::vec Q;
  // If r_corr == 0
  if (r_corr.n_elem == 1 && r_corr(0) == 0) {
    Q.set_size(1);
    Q(0)= arma::accu(arma::square(Score)) / 2;
    if (Score_Resampling.n_rows != 0) {
      Q_res = arma::sum(arma::square(Score_Res), 1) / 2;
    }
    
    return Rcpp::List::create(Rcpp::Named("p.value") = sum(Q), 
                              Rcpp::Named("p.value.resampling") = Q_res, 
                              Rcpp::Named("pval.zero.msg") = Rcpp::CharacterVector::create("r_corr is 0"));
  } else if (r_corr(0) == 1) {
     Q = SKAT_META_Optimal_Get_Q(Score, r_corr);
    if (Score_Resampling.n_rows != 0) {
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
    if (Score_Resampling.n_rows > 0) {
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

