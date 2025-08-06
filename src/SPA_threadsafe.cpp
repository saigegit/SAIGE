// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppArmadillo)]]

#include "SPA_threadsafe.hpp"
#include "SPA_binary_threadsafe.hpp"
#include "SPA_survival_threadsafe.hpp"
#include "UTIL.hpp"

#if DEBUG_MODE
#include <iostream>
#endif

#include <stdexcept>
#include <sstream>

// Thread-safe wrapper functions

RootResult getroot_K1_safe(const std::string& traitType, double init,
                          const arma::vec& mu, const arma::vec& g,
                          double q, double tol) {
    SPA_ASSERT(mu.n_elem == g.n_elem, "mu and g length mismatch in getroot_K1_safe");

    if (traitType == "binary") {
        return getroot_K1_Binom_safe(init, mu, g, q, tol);
    } else if (traitType == "timeToEvent") {
        return getroot_K1_Poi_safe(init, mu, g, q, tol);
    }

    return RootResult();
}

RootResult getroot_K1_fast_safe(const std::string& traitType, double init,
                               const arma::vec& mu, const arma::vec& g, double q,
                               const arma::vec& gNA, const arma::vec& gNB,
                               const arma::vec& muNA, const arma::vec& muNB,
                               double NAmu, double NAsigma, double tol) {
    SPA_ASSERT(mu.n_elem == g.n_elem, "mu and g length mismatch in getroot_K1_fast_safe");

    if (traitType == "binary") {
        return getroot_K1_fast_Binom_safe(init, mu, g, q, gNA, gNB, muNA, muNB, NAmu, NAsigma, tol);
    } else if (traitType == "timeToEvent") {
        return getroot_K1_fast_Poi_safe(init, mu, g, q, gNA, gNB, muNA, muNB, NAmu, NAsigma, tol);
    }

    return RootResult();
}

SaddleResult get_saddle_prob_safe(const std::string& traitType, double root,
                                 const arma::vec& mu, const arma::vec& g,
                                 double q, bool logp) {
    SPA_ASSERT(mu.n_elem == g.n_elem, "mu and g length mismatch in get_saddle_prob_safe");

    if (traitType == "binary") {
        return Get_Saddle_Prob_Binom_safe(root, mu, g, q, logp);
    } else if (traitType == "timeToEvent") {
        return Get_Saddle_Prob_Poi_safe(root, mu, g, q, logp);
    }

    return SaddleResult();
}

SaddleResult get_saddle_prob_fast_safe(const std::string& traitType, double root,
                                      const arma::vec& mu, const arma::vec& g, double q,
                                      const arma::vec& gNA, const arma::vec& gNB,
                                      const arma::vec& muNA, const arma::vec& muNB,
                                      double NAmu, double NAsigma, bool logp) {
    SPA_ASSERT(mu.n_elem == g.n_elem, "mu and g length mismatch in get_saddle_prob_fast_safe");

    if (traitType == "binary") {
        return Get_Saddle_Prob_fast_Binom_safe(root, mu, g, q, gNA, gNB, muNA, muNB, NAmu, NAsigma, logp);
    } else if (traitType == "timeToEvent") {
        return Get_Saddle_Prob_fast_Poi_safe(root, mu, g, q, gNA, gNB, muNA, muNB, NAmu, NAsigma, logp);
    }

    return SaddleResult();
}

// [[Rcpp::export]]
SPAResult SPA_threadsafe(const arma::vec& mu, const arma::vec& g, double q, double qinv,
                         double pval_noadj, double tol, bool logp, const std::string& traitType) {

    SPA_ASSERT(mu.n_elem == g.n_elem, "mu and g must be the same length in SPA_threadsafe");

    RootResult outuni1 = getroot_K1_safe(traitType, 0, mu, g, q, tol);
    RootResult outuni2 = getroot_K1_safe(traitType, 0, mu, g, qinv, tol);
    double p1 = 0.0, p2 = 0.0;
    bool isConverge = true;

    if (outuni1.isConverge && outuni2.isConverge) {
        SaddleResult getSaddle = get_saddle_prob_safe(traitType, outuni1.root, mu, g, q, logp);
        SaddleResult getSaddle2 = get_saddle_prob_safe(traitType, outuni2.root, mu, g, qinv, logp);

        p1 = getSaddle.isSaddle  ? getSaddle.pval  : (logp ? pval_noadj - std::log(2) : pval_noadj / 2);
        p2 = getSaddle2.isSaddle ? getSaddle2.pval : (logp ? pval_noadj - std::log(2) : pval_noadj / 2);
        isConverge = getSaddle.isSaddle && getSaddle2.isSaddle;

        return SPAResult(logp ? add_logp(p1, p2) : std::abs(p1) + std::abs(p2), isConverge);
    }

    return SPAResult(pval_noadj, false);
}

// [[Rcpp::export]]
SPAResult SPA_fast_threadsafe(const arma::vec& mu, const arma::vec& g, double q, double qinv,
                              double pval_noadj, bool logp,
                              const arma::vec& gNA, const arma::vec& gNB,
                              const arma::vec& muNA, const arma::vec& muNB,
                              double NAmu, double NAsigma, double tol,
                              const std::string& traitType) {
/*
	std::cout << "mu.n_elem " << mu.n_elem << std::endl;
        std::cout << "g.n_elem " << g.n_elem << std::endl;
        std::cout << "gNA.n_elem " << gNA.n_elem << std::endl;
        std::cout << "gNB.n_elem " << gNB.n_elem << std::endl;
        std::cout << "muNA.n_elem " << muNA.n_elem << std::endl;
        std::cout << "muNB.n_elem " << muNB.n_elem << std::endl;
*/
    SPA_ASSERT(mu.n_elem == g.n_elem, "mu and g mismatch in SPA_fast_threadsafe");
    SPA_ASSERT(gNA.n_elem == muNA.n_elem, "gNA and muNA must match");
    SPA_ASSERT(gNB.n_elem == muNB.n_elem, "gNB and muNB must match");

    RootResult outuni1 = getroot_K1_fast_safe(traitType, 0, mu, g, q, gNA, gNB, muNA, muNB, NAmu, NAsigma, tol);
    RootResult outuni2 = getroot_K1_fast_safe(traitType, 0, mu, g, qinv, gNA, gNB, muNA, muNB, NAmu, NAsigma, tol);
    double p1 = 0.0, p2 = 0.0;
    bool isConverge = true;

    if (outuni1.isConverge && outuni2.isConverge) {
        SaddleResult getSaddle = get_saddle_prob_fast_safe(traitType, outuni1.root, mu, g, q,
                                                          gNA, gNB, muNA, muNB, NAmu, NAsigma, logp);
        SaddleResult getSaddle2 = get_saddle_prob_fast_safe(traitType, outuni2.root, mu, g, qinv,
                                                           gNA, gNB, muNA, muNB, NAmu, NAsigma, logp);

        p1 = getSaddle.isSaddle  ? getSaddle.pval  : (logp ? pval_noadj - std::log(2) : pval_noadj / 2);
        p2 = getSaddle2.isSaddle ? getSaddle2.pval : (logp ? pval_noadj - std::log(2) : pval_noadj / 2);
        isConverge = getSaddle.isSaddle && getSaddle2.isSaddle;

        return SPAResult(logp ? add_logp(p1, p2) : std::abs(p1) + std::abs(p2), isConverge);
    }

    return SPAResult(pval_noadj, false);
}

// [[Rcpp::export]]
double SPA_pval_threadsafe(const arma::vec& mu, const arma::vec& g, double q, double qinv,
                           double pval_noadj, double tol, bool logp,
                           const std::string& traitType, bool& isSPAConverge) {
    SPAResult res = SPA_threadsafe(mu, g, q, qinv, pval_noadj, tol, logp, traitType);
    isSPAConverge = res.isConverge;
    return res.pvalue;
}
