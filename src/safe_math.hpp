#ifndef SAFE_MATH_HPP
#define SAFE_MATH_HPP

#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <cmath>
#include <stdexcept>

// Thread-safe wrapper for qnorm (inverse normal CDF)
inline double qnorm_cpp(double p, bool lower_tail = true, bool log_p = false) {
    if (log_p) {
        // Value is on log scale: need to exp before passing to quantile
        p = std::exp(p);
    }

    boost::math::normal norm_dist(0.0, 1.0);

    if (!lower_tail) {
        return boost::math::quantile(boost::math::complement(norm_dist, p));
    } else {
        return boost::math::quantile(norm_dist, p);
    }
}

// Thread-safe wrapper for pnorm (normal CDF)
inline double pnorm_cpp(double x, bool lower_tail = true, bool log_p = false) {
    boost::math::normal norm_dist(0.0, 1.0);

    double p = lower_tail
        ? boost::math::cdf(norm_dist, x)
        : boost::math::cdf(boost::math::complement(norm_dist, x));

    if (log_p) {
        return std::log(p);
    } else {
        return p;
    }
}

// Thread-safe wrapper for pchisq (chi-squared CDF)
/*
 * @param x: test statistic
 * @param df: degrees of freedom
 * @param lower_tail: whether to return P[X <= x] (otherwise, return P[X > x])
 * @param log_p: whether to return on log scale
 */
inline double pchisq_cpp(double x, double df, bool lower_tail = true, bool log_p = false) {
    boost::math::chi_squared chi_dist(df);

    double p = lower_tail
        ? boost::math::cdf(chi_dist, x)
        : boost::math::cdf(boost::math::complement(chi_dist, x));

    if (log_p) {
        return std::log(p);
    } else {
        return p;
    }
}

#endif // SAFE_MATH_HPP
