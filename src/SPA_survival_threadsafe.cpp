// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppArmadillo)]]

#include "SPA_threadsafe.hpp"
#include "SPA_survival_threadsafe.hpp"
#include "UTIL.hpp"

#include <boost/math/distributions/normal.hpp>
#include <cmath>

#if DEBUG_MODE
#include <iostream>
#endif

// Thread-local normal distribution
static thread_local boost::math::normal normal_dist(0.0, 1.0);

double thread_safe_pnorm(double x, bool lower_tail = true, bool log_p = false) {
    try {
        double p = lower_tail
            ? boost::math::cdf(normal_dist, x)
            : boost::math::cdf(boost::math::complement(normal_dist, x));
        return log_p ? std::log(p) : p;
    } catch (...) {
        double r = 0.5 * (1 + std::erf(x / std::sqrt(2.0)));
        if (!lower_tail) r = 1.0 - r;
        return log_p ? std::log(r) : r;
    }
}

// -----------------------------------
// Core SPA functions (Poisson case)
// -----------------------------------

double Korg_Poi(double t1, const arma::vec& mu, const arma::vec& g) {
    SPA_ASSERT(mu.n_elem == g.n_elem, "Korg_Poi: mu and g size mismatch");
    return arma::dot(mu, arma::exp(g * t1) - g * t1 - 1);
}

double K1_adj_Poi(double t1, const arma::vec& mu, const arma::vec& g, double q) {
    SPA_ASSERT(mu.n_elem == g.n_elem, "K1_adj_Poi: mu and g size mismatch");
    return arma::sum(mu % g % (arma::exp(g * t1) - 1)) - q;
}

double K2_Poi(double t1, const arma::vec& mu, const arma::vec& g) {
    SPA_ASSERT(mu.n_elem == g.n_elem, "K2_Poi: mu and g size mismatch");
    return arma::sum(mu % arma::pow(g, 2) % arma::exp(g * t1));
}

RootResult getroot_K1_Poi_safe(double init, const arma::vec& mu, const arma::vec& g,
                                double q, double tol, int maxiter) {
    double t = init;
    double prevJump = std::numeric_limits<double>::infinity();
    double K1_eval = K1_adj_Poi(t, mu, g, q);

    int rep = 1;
    bool conv = true;

    while (rep <= maxiter) {
        double K2_eval = K2_Poi(t, mu, g);
        double tnew = t - K1_eval / K2_eval;

        if (!std::isfinite(tnew)) { conv = false; break; }
        if (std::abs(tnew - t) < tol) break;
        if (rep == maxiter) { conv = false; break; }

        double newK1 = K1_adj_Poi(tnew, mu, g, q);
        if (arma::sign(K1_eval) != arma::sign(newK1)) {
            if (std::abs(tnew - t) > (prevJump - tol)) {
                tnew = t + arma::sign(newK1 - K1_eval) * prevJump / 2;
                newK1 = K1_adj_Poi(tnew, mu, g, q);
                prevJump /= 2;
            } else {
                prevJump = std::abs(tnew - t);
            }
        }

        t = tnew;
        K1_eval = newK1;
        ++rep;
    }

    return RootResult(t, conv);
}

SaddleResult Get_Saddle_Prob_Poi_safe(double zeta, const arma::vec& mu, const arma::vec& g,
                                      double q, bool logp) {
    SPA_ASSERT(mu.n_elem == g.n_elem, "Get_Saddle_Prob_Poi_safe: mu and g size mismatch");

    double k1 = Korg_Poi(zeta, mu, g);
    double k2 = K2_Poi(zeta, mu, g);
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

SPAResult SPA_survival_safe(const arma::vec& mu, const arma::vec& g,
                            double q, double qinv,
                            double pval_noadj, double tol, bool logp) {
    SPA_ASSERT(mu.n_elem == g.n_elem, "SPA_survival_safe: mu and g size mismatch");

    RootResult out1 = getroot_K1_Poi_safe(0, mu, g, q, tol);
    RootResult out2 = getroot_K1_Poi_safe(0, mu, g, qinv, tol);

    if (!(out1.isConverge && out2.isConverge))
        return SPAResult(pval_noadj, false);

    SaddleResult s1 = Get_Saddle_Prob_Poi_safe(out1.root, mu, g, q, logp);
    SaddleResult s2 = Get_Saddle_Prob_Poi_safe(out2.root, mu, g, qinv, logp);

    double p1 = s1.isSaddle ? s1.pval : (logp ? pval_noadj - std::log(2) : pval_noadj / 2);
    double p2 = s2.isSaddle ? s2.pval : (logp ? pval_noadj - std::log(2) : pval_noadj / 2);

    bool converge = s1.isSaddle && s2.isSaddle;
    double finalP = logp ? add_logp(p1, p2) : std::abs(p1) + std::abs(p2);

    return SPAResult(finalP, converge);
}

// -----------------------------------
// Fast survival SPA (NA-aware case)
// -----------------------------------

double Korg_fast_Poi(double t1, const arma::vec& mu, const arma::vec& g,
                     const arma::vec& gNA, const arma::vec& gNB,
                     const arma::vec& muNA, const arma::vec& muNB,
                     double NAmu, double NAsigma) {
    SPA_ASSERT(gNB.n_elem == muNB.n_elem, "Korg_fast_Poi: gNB and muNB size mismatch");
    arma::vec temp = muNB % (arma::exp(gNB * t1) - gNB * t1 - 1);
    return arma::sum(temp) + 0.5 * NAsigma * std::pow(t1, 2);
}

double K1_adj_fast_Poi(double t1, const arma::vec& mu, const arma::vec& g, double q,
                       const arma::vec& gNA, const arma::vec& gNB,
                       const arma::vec& muNA, const arma::vec& muNB,
                       double NAmu, double NAsigma) {
    SPA_ASSERT(gNB.n_elem == muNB.n_elem, "K1_adj_fast_Poi: gNB and muNB size mismatch");

    arma::vec temp = muNB % gNB % (arma::exp(gNB * t1) - 1);
    return arma::sum(temp) + NAsigma * t1 - q;
}

double K2_fast_Poi(double t1, const arma::vec& mu, const arma::vec& g,
                   const arma::vec& gNA, const arma::vec& gNB,
                   const arma::vec& muNA, const arma::vec& muNB,
                   double NAmu, double NAsigma) {
    SPA_ASSERT(gNB.n_elem == muNB.n_elem, "K2_fast_Poi: gNB and muNB size mismatch");

    arma::vec temp = muNB % arma::pow(gNB, 2) % arma::exp(gNB * t1);
    return arma::sum(temp) + NAsigma;
}

RootResult getroot_K1_fast_Poi_safe(double init, const arma::vec& mu, const arma::vec& g, double q,
                                    const arma::vec& gNA, const arma::vec& gNB,
                                    const arma::vec& muNA, const arma::vec& muNB,
                                    double NAmu, double NAsigma,
                                    double tol, int maxiter) {
    double t = init;
    double K1_eval = K1_adj_fast_Poi(t, mu, g, q, gNA, gNB, muNA, muNB, NAmu, NAsigma);
    double prevJump = std::numeric_limits<double>::infinity();
    int rep = 1;
    bool conv = true;

    while (rep <= maxiter) {
        double K2_eval = K2_fast_Poi(t, mu, g, gNA, gNB, muNA, muNB, NAmu, NAsigma);
        double tnew = t - K1_eval / K2_eval;

        if (!std::isfinite(tnew)) { conv = false; break; }
        if (std::abs(tnew - t) < tol) break;
        if (rep == maxiter) { conv = false; break; }

        double newK1 = K1_adj_fast_Poi(tnew, mu, g, q, gNA, gNB, muNA, muNB, NAmu, NAsigma);
        if ((K1_eval * newK1) < 0) {
            if (std::abs(tnew - t) > (prevJump - tol)) {
                tnew = t + arma::sign(newK1 - K1_eval) * prevJump / 2;
                newK1 = K1_adj_fast_Poi(tnew, mu, g, q, gNA, gNB, muNA, muNB, NAmu, NAsigma);
                prevJump /= 2;
            } else {
                prevJump = std::abs(tnew - t);
            }
        }

        t = tnew;
        K1_eval = newK1;
        ++rep;
    }

    return RootResult(t, conv);
}

SaddleResult Get_Saddle_Prob_fast_Poi_safe(double zeta, const arma::vec& mu, const arma::vec& g, double q,
                                           const arma::vec& gNA, const arma::vec& gNB,
                                           const arma::vec& muNA, const arma::vec& muNB,
                                           double NAmu, double NAsigma, bool logp) {
    double k1 = Korg_fast_Poi(zeta, mu, g, gNA, gNB, muNA, muNB, NAmu, NAsigma);
    double k2 = K2_fast_Poi(zeta, mu, g, gNA, gNB, muNA, muNB, NAmu, NAsigma);
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

SPAResult SPA_survival_fast_safe(const arma::vec& mu, const arma::vec& g, double q, double qinv,
                                 double pval_noadj, bool logp,
                                 const arma::vec& gNA, const arma::vec& gNB,
                                 const arma::vec& muNA, const arma::vec& muNB,
                                 double NAmu, double NAsigma, double tol) {
    RootResult out1 = getroot_K1_fast_Poi_safe(0, mu, g, q, gNA, gNB, muNA, muNB, NAmu, NAsigma, tol);
    RootResult out2 = getroot_K1_fast_Poi_safe(0, mu, g, qinv, gNA, gNB, muNA, muNB, NAmu, NAsigma, tol);

    if (!(out1.isConverge && out2.isConverge))
        return SPAResult(pval_noadj, false);

    SaddleResult s1 = Get_Saddle_Prob_fast_Poi_safe(out1.root, mu, g, q, gNA, gNB, muNA, muNB, NAmu, NAsigma, logp);
    SaddleResult s2 = Get_Saddle_Prob_fast_Poi_safe(out2.root, mu, g, qinv, gNA, gNB, muNA, muNB, NAmu, NAsigma, logp);

    double p1 = s1.isSaddle ? s1.pval : (logp ? pval_noadj - std::log(2) : pval_noadj / 2);
    double p2 = s2.isSaddle ? s2.pval : (logp ? pval_noadj - std::log(2) : pval_noadj / 2);

    bool isConverge = s1.isSaddle && s2.isSaddle;
    double finalP = logp ? add_logp(p1, p2) : std::abs(p1) + std::abs(p2);

    return SPAResult(finalP, isConverge);
}

