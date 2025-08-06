#ifndef SPA_SURVIVAL_THREADSAFE_HPP
#define SPA_SURVIVAL_THREADSAFE_HPP

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// Pull in result structs and assertions from main SPA header
#include "SPA_threadsafe.hpp"

// ---------- Core Poisson-style SPA Functions (Standard) ----------

double Korg_Poi(double t1, const arma::vec& mu, const arma::vec& g);
double K1_adj_Poi(double t1, const arma::vec& mu, const arma::vec& g, double q);
double K2_Poi(double t1, const arma::vec& mu, const arma::vec& g);

// Newton root solver for survival SPA
RootResult getroot_K1_Poi_safe(double init,
                               const arma::vec& mu,
                               const arma::vec& g,
                               double q,
                               double tol,
                               int maxiter = 1000);

// Saddle-point probability (survival)
SaddleResult Get_Saddle_Prob_Poi_safe(double zeta,
                                      const arma::vec& mu,
                                      const arma::vec& g,
                                      double q,
                                      bool logp = false);

// Top-level SPA for survival
SPAResult SPA_survival_safe(const arma::vec& mu,
                            const arma::vec& g,
                            double q,
                            double qinv,
                            double pval_noadj,
                            double tol,
                            bool logp = false);

// ---------- FAST SPA: Missing/stratified adjustments ----------

double Korg_fast_Poi(double t1,
                     const arma::vec& mu,
                     const arma::vec& g,
                     const arma::vec& gNA,
                     const arma::vec& gNB,
                     const arma::vec& muNA,
                     const arma::vec& muNB,
                     double NAmu,
                     double NAsigma);

double K1_adj_fast_Poi(double t1,
                       const arma::vec& mu,
                       const arma::vec& g,
                       double q,
                       const arma::vec& gNA,
                       const arma::vec& gNB,
                       const arma::vec& muNA,
                       const arma::vec& muNB,
                       double NAmu,
                       double NAsigma);

double K2_fast_Poi(double t1,
                   const arma::vec& mu,
                   const arma::vec& g,
                   const arma::vec& gNA,
                   const arma::vec& gNB,
                   const arma::vec& muNA,
                   const arma::vec& muNB,
                   double NAmu,
                   double NAsigma);

// Root solving with NA adjustment
RootResult getroot_K1_fast_Poi_safe(double init,
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

// Saddle point with NA adjustment
SaddleResult Get_Saddle_Prob_fast_Poi_safe(double zeta,
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

// Final SPA entry point for fast Poisson-based logic
SPAResult SPA_survival_fast_safe(const arma::vec& mu,
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

#endif  // SPA_SURVIVAL_THREADSAFE_HPP

