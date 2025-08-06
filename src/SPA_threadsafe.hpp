#ifndef SPA_THREADSAFE_HPP
#define SPA_THREADSAFE_HPP

#include <RcppArmadillo.h>
#include <string>
#include <sstream>
#include <stdexcept>

// Enable strict runtime checks for OpenMP/SPA routines
#ifndef DEBUG_MODE
#define DEBUG_MODE 1
#endif

#if DEBUG_MODE
#include <iostream>
#define SPA_ASSERT(cond, msg) \
    do { \
        if (!(cond)) { \
            std::ostringstream _oss; \
            _oss << "DEBUG ASSERTION FAILED: " << msg \
                 << " (" << __FILE__ << ":" << __LINE__ << ")"; \
            throw std::runtime_error(_oss.str()); \
        } \
    } while (0)
#else
#define SPA_ASSERT(cond, msg)
#endif

struct RootResult {
    double root;
    bool isConverge;
    RootResult() : root(0.0), isConverge(false) {}
    RootResult(double r, bool conv) : root(r), isConverge(conv) {}
};

struct SaddleResult {
    double pval;
    bool isSaddle;
    SaddleResult() : pval(0.0), isSaddle(false) {}
    SaddleResult(double p, bool saddle) : pval(p), isSaddle(saddle) {}
};

struct SPAResult {
    double pvalue;
    bool isConverge;
    SPAResult() : pvalue(0.0), isConverge(false) {}
    SPAResult(double p, bool conv) : pvalue(p), isConverge(conv) {}
};

// Wrappers for root and saddle computation (you should implement or link these)
RootResult getroot_K1_safe(const std::string& traitType, double init,
                           const arma::vec& mu, const arma::vec& g, double q, double tol);

RootResult getroot_K1_fast_safe(const std::string& traitType, double init,
                                const arma::vec& mu, const arma::vec& g, double q,
                                const arma::vec& gNA, const arma::vec& gNB, const arma::vec& muNA, const arma::vec& muNB,
                                double NAmu, double NAsigma, double tol);

SaddleResult get_saddle_prob_safe(const std::string& traitType, double root,
                                  const arma::vec& mu, const arma::vec& g, double q, bool logp);

SaddleResult get_saddle_prob_fast_safe(const std::string& traitType, double root,
                                       const arma::vec& mu, const arma::vec& g, double q,
                                       const arma::vec& gNA, const arma::vec& gNB, const arma::vec& muNA, const arma::vec& muNB,
                                       double NAmu, double NAsigma, bool logp);

// Thread-safe SPA main functions
SPAResult SPA_threadsafe(const arma::vec& mu, const arma::vec& g, double q, double qinv,
                         double pval_noadj, double tol, bool logp, const std::string& traitType);

SPAResult SPA_fast_threadsafe(const arma::vec& mu, const arma::vec& g, double q, double qinv,
                              double pval_noadj, bool logp,
                              const arma::vec& gNA, const arma::vec& gNB, const arma::vec& muNA, const arma::vec& muNB,
                              double NAmu, double NAsigma, double tol, const std::string& traitType);

#endif // SPA_THREADSAFE_HPP
