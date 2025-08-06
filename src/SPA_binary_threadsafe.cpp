// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppArmadillo)]]

#include "SPA_threadsafe.hpp"       // Structs: RootResult, SaddleResult, SPAResult, SPA_ASSERT
#include "SPA_binary_threadsafe.hpp" // Function declarations
#include "UTIL.hpp"                  // For add_logp helper or others

#include <boost/math/distributions/normal.hpp>
#include <cmath>
#include <stdexcept>

#if DEBUG_MODE
#include <iostream>
#endif

// Thread-local standard normal distribution for pnorm
static thread_local boost::math::normal normal_dist(0.0, 1.0);

// Thread-safe normal pnorm with fallback
inline double thread_safe_pnorm(double x, bool lower_tail = true, bool log_p = false) {
    try {
        double p = lower_tail
            ? boost::math::cdf(normal_dist, x)
            : boost::math::cdf(boost::math::complement(normal_dist, x));
        return log_p ? std::log(p) : p;
    } catch (...) {
        // Fallback using std::erf for numerical stability
        double result = 0.5 * (1.0 + std::erf(x / std::sqrt(2.0)));
        if (!lower_tail) result = 1.0 - result;
        return log_p ? std::log(result) : result;
    }
}

// ---- Standard Binomial SPA Helper Functions ----

double Korg_Binom_safe(double t1, const arma::vec& mu, const arma::vec& g) {
    SPA_ASSERT(mu.n_elem == g.n_elem, "Korg_Binom_safe: mu/g length mismatch");
    return arma::sum(arma::log(1 - mu + mu % arma::exp(g * t1)));
}

double K1_adj_Binom_safe(double t1, const arma::vec& mu, const arma::vec& g, double q) {
    SPA_ASSERT(mu.n_elem == g.n_elem, "K1_adj_Binom_safe: mu/g length mismatch");
    arma::vec temp1 = (1 - mu) % arma::exp(-g * t1) + mu;
    arma::vec temp2 = mu % g;
    return arma::sum(temp2 / temp1) - q;
}

double K2_Binom_safe(double t1, const arma::vec& mu, const arma::vec& g) {
    SPA_ASSERT(mu.n_elem == g.n_elem, "K2_Binom_safe: mu/g length mismatch");
    arma::vec temp0 = arma::exp(-g * t1);
    arma::vec temp1 = arma::pow((1 - mu) % temp0 + mu, 2);
    arma::vec temp2 = (1 - mu) % mu % arma::pow(g, 2) % temp0;
    return arma::sum(temp2 / temp1);
}

RootResult getroot_K1_Binom_safe(double init, const arma::vec& mu, const arma::vec& g,
                                 double q, double tol, int maxiter) {
    SPA_ASSERT(mu.n_elem == g.n_elem, "getroot_K1_Binom_safe: mu/g length mismatch");

    double gpos = arma::accu(g.elem(arma::find(g > 0)));
    double gneg = arma::accu(g.elem(arma::find(g < 0)));

    if (q >= gpos || q <= gneg)
        return RootResult(std::numeric_limits<double>::infinity(), true);

    double t = init;
    double prevJump = std::numeric_limits<double>::infinity();
    double K1_eval = K1_adj_Binom_safe(t, mu, g, q);
    int rep = 1;
    bool conv = true;

    while (rep <= maxiter) {
        double K2_eval = K2_Binom_safe(t, mu, g);
        double tnew = t - K1_eval / K2_eval;

        if (!std::isfinite(tnew)) { conv = false; break; }
        if (std::abs(tnew - t) < tol) break;
        if (rep == maxiter) { conv = false; break; }

        double newK1 = K1_adj_Binom_safe(tnew, mu, g, q);
        if (arma::sign(K1_eval) != arma::sign(newK1)) {
            if (std::abs(tnew - t) > (prevJump - tol)) {
                tnew = t + arma::sign(newK1 - K1_eval) * prevJump / 2;
                newK1 = K1_adj_Binom_safe(tnew, mu, g, q);
                prevJump /= 2;
            } else {
                prevJump = std::abs(tnew - t);
            }
        }
        rep++; t = tnew; K1_eval = newK1;
    }
    return RootResult(t, conv);
}

SaddleResult Get_Saddle_Prob_Binom_safe(double zeta, const arma::vec& mu, const arma::vec& g,
                                        double q, bool logp) {
    SPA_ASSERT(mu.n_elem == g.n_elem, "Get_Saddle_Prob_Binom_safe: mu/g length mismatch");

    double k1 = Korg_Binom_safe(zeta, mu, g);
    double k2 = K2_Binom_safe(zeta, mu, g);
    double temp1 = zeta * q - k1;

    if (std::isfinite(k1) && std::isfinite(k2) && temp1 >= 0 && k2 >= 0) {
        double w = arma::sign(zeta) * std::sqrt(2 * temp1);
        double v = zeta * std::sqrt(k2);
        if (w != 0) {
            double Ztest = w + (1 / w) * std::log(v / w);
            double pval = (Ztest > 0) ? thread_safe_pnorm(Ztest, false, logp)
                                      : -thread_safe_pnorm(Ztest, true, logp);
            return SaddleResult(pval, true);
        }
    }
    return SaddleResult(logp ? -std::numeric_limits<double>::infinity() : 0.0, false);
}

SPAResult SPA_binary_safe(const arma::vec& mu, const arma::vec& g,
                          double q, double qinv, double pval_noadj,
                          double tol, bool logp) {
    SPA_ASSERT(mu.n_elem == g.n_elem, "SPA_binary_safe: mu/g length mismatch");

    RootResult outuni1 = getroot_K1_Binom_safe(0, mu, g, q, tol);
    RootResult outuni2 = getroot_K1_Binom_safe(0, mu, g, qinv, tol);

    if (!(outuni1.isConverge && outuni2.isConverge))
        return SPAResult(pval_noadj, false);

    SaddleResult s1 = Get_Saddle_Prob_Binom_safe(outuni1.root, mu, g, q, logp);
    SaddleResult s2 = Get_Saddle_Prob_Binom_safe(outuni2.root, mu, g, qinv, logp);

    double p1 = s1.isSaddle ? s1.pval : (logp ? pval_noadj - std::log(2) : pval_noadj / 2);
    double p2 = s2.isSaddle ? s2.pval : (logp ? pval_noadj - std::log(2) : pval_noadj / 2);
    bool isConverge = s1.isSaddle && s2.isSaddle;

    double pval = logp ? add_logp(p1, p2) : std::abs(p1) + std::abs(p2);
    return SPAResult(pval, isConverge);
}

// ---- FAST BINOMIAL FUNCTIONS: NA-AWARE ----

double Korg_fast_Binom_safe(double t1, const arma::vec& mu,
    const arma::vec& g,
    const arma::vec& gNA,
    const arma::vec& gNB,
    const arma::vec& muNA,
    const arma::vec& muNB,
    double NAmu, double NAsigma) {

    SPA_ASSERT(muNB.n_elem == gNB.n_elem, "Korg_fast_Binom_safe: muNB and gNB length mismatch");

    arma::vec temp = arma::log(1 - muNB + muNB % arma::exp(gNB * t1));
    return arma::sum(temp) + NAmu * t1 + 0.5 * NAsigma * std::pow(t1, 2);
}

double K1_adj_fast_Binom_safe(double t1, const arma::vec& mu,
    const arma::vec& g, double q,
    const arma::vec& gNA,
    const arma::vec& gNB,
    const arma::vec& muNA,
    const arma::vec& muNB,
    double NAmu, double NAsigma) {

    SPA_ASSERT(muNB.n_elem == gNB.n_elem, "K1_adj_fast_Binom_safe: muNB and gNB length mismatch");

    arma::vec temp1 = (1 - muNB) % arma::exp(-gNB * t1) + muNB;
    arma::vec temp2 = muNB % gNB;
    double temp3 = NAmu + NAsigma * t1;
    return arma::sum(temp2 / temp1) + temp3 - q;
}

double K2_fast_Binom_safe(double t1, const arma::vec& mu,
    const arma::vec& g,
    const arma::vec& gNA,
    const arma::vec& gNB,
    const arma::vec& muNA,
    const arma::vec& muNB,
    double NAmu, double NAsigma) {

    SPA_ASSERT(muNB.n_elem == gNB.n_elem, "K2_fast_Binom_safe: muNB and gNB length mismatch");

    arma::vec temp0 = arma::exp(-gNB * t1);
    arma::vec temp1 = (1 - muNB) % temp0 + muNB;
    temp1 = arma::pow(temp1, 2);
    arma::vec temp2 = arma::pow(gNB, 2) % temp0;
    temp2 = (1 - muNB) % muNB % temp2;
    return arma::sum(temp2 / temp1) + NAsigma;
}

RootResult getroot_K1_fast_Binom_safe(double init,
    const arma::vec& mu, const arma::vec& g, double q,
    const arma::vec& gNA,
    const arma::vec& gNB,
    const arma::vec& muNA,
    const arma::vec& muNB,
    double NAmu, double NAsigma,
    double tol, int maxiter) {

    SPA_ASSERT(muNB.n_elem == gNB.n_elem, "getroot_K1_fast_Binom_safe: muNB and gNB length mismatch");

    double gpos = arma::accu(g.elem(arma::find(g > 0)));
    double gneg = arma::accu(g.elem(arma::find(g < 0)));

    if (q >= gpos || q <= gneg) {
        return RootResult(std::numeric_limits<double>::infinity(), true);
    }

    double t = init;
    double K1_eval = K1_adj_fast_Binom_safe(t, mu, g, q, gNA, gNB, muNA, muNB, NAmu, NAsigma);
    double prevJump = std::numeric_limits<double>::infinity();
    int rep = 1;
    bool conv = true;

    while (rep <= maxiter) {
        double K2_eval = K2_fast_Binom_safe(t, mu, g, gNA, gNB, muNA, muNB, NAmu, NAsigma);
        double tnew = t - K1_eval / K2_eval;

        if (!std::isfinite(tnew)) { conv = false; break; }
        if (std::abs(tnew - t) < tol) break;
        if (rep == maxiter) { conv = false; break; }

        double newK1 = K1_adj_fast_Binom_safe(tnew, mu, g, q, gNA, gNB, muNA, muNB, NAmu, NAsigma);
        if ((K1_eval * newK1) < 0) {
            if (std::abs(tnew - t) > (prevJump - tol)) {
                tnew = t + arma::sign(newK1 - K1_eval) * prevJump / 2;
                newK1 = K1_adj_fast_Binom_safe(tnew, mu, g, q, gNA, gNB, muNA, muNB, NAmu, NAsigma);
                prevJump /= 2;
            } else {
                prevJump = std::abs(tnew - t);
            }
        }

        rep++; t = tnew; K1_eval = newK1;
    }

    //return RootResult(t, rep, conv);
    return RootResult(t, conv);
}

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
                                             bool logp) {

    SPA_ASSERT(muNB.n_elem == gNB.n_elem, "Get_Saddle_Prob_fast_Binom_safe: muNB and gNB length mismatch");

    double k1 = Korg_fast_Binom_safe(zeta, mu, g, gNA, gNB, muNA, muNB, NAmu, NAsigma);
    double k2 = K2_fast_Binom_safe(zeta, mu, g, gNA, gNB, muNA, muNB, NAmu, NAsigma);
    double temp1 = zeta * q - k1;

    if (std::isfinite(k1) && std::isfinite(k2) && temp1 >= 0 && k2 >= 0) {
        double w = arma::sign(zeta) * std::sqrt(2 * temp1);
        double v = zeta * std::sqrt(k2);
        if (w != 0) {
            double Ztest = w + (1 / w) * std::log(v / w);
            double pval = (Ztest > 0)
                ? thread_safe_pnorm(Ztest, false, logp)
                : -thread_safe_pnorm(Ztest, true, logp);
            return SaddleResult(pval, true);
        }
    }
    return SaddleResult(logp ? -std::numeric_limits<double>::infinity() : 0.0, false);
}

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
                               double tol) {

    SPA_ASSERT(mu.n_elem == g.n_elem, "SPA_binary_fast_safe: mu and g length mismatch");
    SPA_ASSERT(gNA.n_elem == gNB.n_elem, "gNA and gNB length mismatch");
    SPA_ASSERT(muNA.n_elem == muNB.n_elem, "muNA and muNB length mismatch");

    RootResult outuni1 = getroot_K1_fast_Binom_safe(0, mu, g, q, gNA, gNB, muNA, muNB, NAmu, NAsigma, tol);
    RootResult outuni2 = getroot_K1_fast_Binom_safe(0, mu, g, qinv, gNA, gNB, muNA, muNB, NAmu, NAsigma, tol);

    if (!outuni1.isConverge || !outuni2.isConverge)
        return SPAResult(pval_noadj, false);

    SaddleResult s1 = Get_Saddle_Prob_fast_Binom_safe(outuni1.root, mu, g, q, gNA, gNB, muNA, muNB, NAmu, NAsigma, logp);
    SaddleResult s2 = Get_Saddle_Prob_fast_Binom_safe(outuni2.root, mu, g, qinv, gNA, gNB, muNA, muNB, NAmu, NAsigma, logp);

    double p1 = s1.isSaddle ? s1.pval : (logp ? pval_noadj - std::log(2) : pval_noadj / 2);
    double p2 = s2.isSaddle ? s2.pval : (logp ? pval_noadj - std::log(2) : pval_noadj / 2);
    bool isConverge = s1.isSaddle && s2.isSaddle;

    double pval = logp ? add_logp(p1, p2) : std::abs(p1) + std::abs(p2);
    return SPAResult(pval, isConverge);
}

