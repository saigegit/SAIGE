#ifndef SPA_BINARY_THREADSAFE_HPP
#define SPA_BINARY_THREADSAFE_HPP

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// Result structs are defined externally (in SPA_threadsafe.hpp)
#include "SPA_threadsafe.hpp"  // for RootResult, SaddleResult, SPAResult (ensures consistent types)

// ---------- Core Binomial-style SPA Functions (Standard) ----------

// Log-likelihood component
double Korg_Binom_safe(double t1, const arma::vec& mu, const arma::vec& g);

// First derivative
double K1_adj_Binom_safe(double t1, const arma::vec& mu, const arma::vec& g, double q);

// Second derivative
double K2_Binom_safe(double t1, const arma::vec& mu, const arma::vec& g);

// Newton root solver — thread-safe binary SPA
RootResult getroot_K1_Binom_safe(double init,
                                 const arma::vec& mu,
                                 const arma::vec& g,
                                 double q,
                                 double tol,
                                 int maxiter = 1000);

// Saddle point p-value calculator
SaddleResult Get_Saddle_Prob_Binom_safe(double zeta,
                                        const arma::vec& mu,
                                        const arma::vec& g,
                                        double q,
                                        bool logp = false);

// SPA driver for binary trait
SPAResult SPA_binary_safe(const arma::vec& mu,
                          const arma::vec& g,
                          double q,
                          double qinv,
                          double pval_noadj,
                          double tol,
                          bool logp = false);

// ---------- Fast SPA — NA/Imputation-Aware Models ----------

// K functions for fast version w/ missing data adjustment
double Korg_fast_Binom_safe(double t1,
                            const arma::vec& mu,
                            const arma::vec& g,
                            const arma::vec& gNA,
                            const arma::vec& gNB,
                            const arma::vec& muNA,
                            const arma::vec& muNB,
                            double NAmu,
                            double NAsigma);

double K1_adj_fast_Binom_safe(double t1,
                              const arma::vec& mu,
                              const arma::vec& g,
                              double q,
                              const arma::vec& gNA,
                              const arma::vec& gNB,
                              const arma::vec& muNA,
                              const arma::vec& muNB,
                              double NAmu,
                              double NAsigma);

double K2_fast_Binom_safe(double t1,
                          const arma::vec& mu,
                          const arma::vec& g,
                          const arma::vec& gNA,
                          const arma::vec& gNB,
                          const arma::vec& muNA,
                          const arma::vec& muNB,
                          double NAmu,
                          double NAsigma);

// Root solver for fast model
RootResult getroot_K1_fast_Binom_safe(double init,
                                      const arma::vec& mu,
                                      const arma::vec& g,
                                      double q,
                                      const arma::vec& gNA,
                                      const arma::vec& gNB,
                                      const arma::vec& muNA,
                                      const arma::vec& muNB,
                                      double NAmu,
                                      double NAsigma,
                                      double tol,
                                      int maxiter = 1000);

// Saddle approximation for fast model
SaddleResult Get_Saddle_Prob_fast_Binom_safe(double zeta,
                                             const arma::vec& mu,
                                             const arma::vec& g,
                                             double q,
                                             const arma::vec& gNA,
                                             const arma::vec& gNB,
                                             const arma::vec& muNA,
                                             const arma::vec& muNB,
                                             double NAmu,
                                             double NAsigma,
                                             bool logp = false);

// Top-level SPA function for fast mode
SPAResult SPA_binary_fast_safe(const arma::vec& mu,
                               const arma::vec& g,
                               double q,
                               double qinv,
                               double pval_noadj,
                               bool logp,
                               const arma::vec& gNA,
                               const arma::vec& gNB,
                               const arma::vec& muNA,
                               const arma::vec& muNB,
                               double NAmu,
                               double NAsigma,
                               double tol);

#endif  // SPA_BINARY_THREADSAFE_HPP

