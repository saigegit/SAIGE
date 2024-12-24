// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>

#include "qfc_rcpp.hpp"


#define TRUE  1
#define FALSE 0

//typedef int BOOL;

//#define pi 3.14159265358979
//#define log28 .0866  /* log(2.0) / 8.0 */

// Class to encapsulate the functionality of the original C code

namespace QUADFORM {

//class QuadForm {
    QuadFormClass::QuadFormClass(int r, int lim) : r(r), lim(lim) {
        count = 0;
        ndtsrt = false;
        fail = false;
        lb.resize(r);
        nc.resize(r);
        th.resize(r);
        n.resize(r);
    }

    // Helper functions
    double QuadFormClass::exp1(double x) {
        return x < -50.0 ? 0.0 : exp(x);
    }

    void QuadFormClass::counter() {
        count++;
        if (count > lim) {
            throw std::runtime_error("Exceeded count limit");
        }
    }

    double QuadFormClass::log1(double x, bool first) {
        if (fabs(x) > 0.1) {
            return (first ? log(1.0 + x) : (log(1.0 + x) - x));
        } else {
            double s, s1, term, y, k;
            y = x / (2.0 + x);
            term = 2.0 * pow(y,3);
            k = 3.0;
            s = (first ? 2.0 : -x) * y;
            y = pow(y,2);
            for (s1 = s + term / k; s1 != s; s1 = s + term / k) {
                k = k + 2.0;
                term = term * y;
                s = s1;
            }
            return s;
        }
    }

    void QuadFormClass::order() {
        for (int j = 0; j < r; j++) {
            double lj = fabs(lb[j]);
            int k;
            for (k = j - 1; k >= 0; k--) {
                if (lj > fabs(lb[th[k]])) {
                    th[k + 1] = th[k];
                } else {
                    break;
                }
            }
            th[k + 1] = j;
        }
        ndtsrt = false;
    }

    double QuadFormClass::errbd(double u, double& cx) {
        double sum1 = 0.0, lj, ncj, x, y, xconst;
        counter();
        xconst = u * sigsq;
        sum1 = u * xconst;
        u = 2.0 * u;

        for (int j = r - 1; j >= 0; j--) {
            int nj = n[j];
            lj = lb[j];
            ncj = nc[j];
            x = u * lj;
            y = 1.0 - x;
            xconst += lj * (ncj / y + nj) / y;
            sum1 += ncj * pow(x / y, 2) + nj * (pow(x, 2) / y + log1(-x, FALSE));
        }

        cx = xconst;
        return exp1(-0.5 * sum1);
    }

    double QuadFormClass::ctff(double accx, double& upn) {
        double u1 = 0.0, u2 = upn, c1 = mean, c2, rb, xconst, u;

        rb = 2.0 * ((u2 > 0.0) ? lmax : lmin);
        while (errbd(u2 / (1.0 + u2 * rb), c2) > accx) {
            u1 = u2;
            c1 = c2;
            u2 = 2.0 * u2;
        }

        while ((c1 - mean) / (c2 - mean) < 0.9) {
            u = (u1 + u2) / 2.0;
            if (errbd(u / (1.0 + u * rb), xconst) > accx) {
                u1 = u;
                c1 = xconst;
            } else {
                u2 = u;
                c2 = xconst;
            }
        }
        upn = u2;
        return c2;
    }

    double QuadFormClass::truncation(double u, double tausq) {
        double sum1 = 0.0, sum2 = 0.0, prod1 = 0.0, prod2 = 0.0, prod3 = 0.0;
        int s = 0;
        counter();
        sum2 = (sigsq + tausq) * pow(u, 2);
        prod1 = 2.0 * sum2;
        u = 2.0 * u;

        for (int j = 0; j < r; j++) {
            double lj = lb[j];
            double ncj = nc[j];
            int nj = n[j];
            double x = pow(u * lj, 2);
            sum1 += ncj * x / (1.0 + x);
            if (x > 1.0) {
                prod2 += nj * log(x);
                prod3 += nj * log1(x, TRUE);
                s += nj;
            } else {
                prod1 += nj * log1(x, TRUE);
            }
        }

        sum1 = 0.5 * sum1;
        prod2 = prod1 + prod2;
        prod3 = prod1 + prod3;
        double x = exp1(-sum1 - 0.25 * prod2) / M_PI;
        double y = exp1(-sum1 - 0.25 * prod3) / M_PI;
        double err1 = (s == 0) ? 1.0 : x * 2.0 / s;
        double err2 = (prod3 > 1.0) ? 2.5 * y : 1.0;
        if (err2 < err1) err1 = err2;
        x = 0.5 * sum2;
        err2 = (x <= y) ? 1.0 : y / x;
        return (err1 < err2) ? err1 : err2;
    }

    void QuadFormClass::findu(double& utx, double accx) {
        double u = utx / 4.0;
        if (truncation(u, 0.0) > accx) {
            for (u = utx; truncation(u, 0.0) > accx; u = utx) utx *= 4.0;
        } else {
            utx = u;
            for (u = u / 4.0; truncation(u, 0.0) <= accx; u = u / 4.0) utx = u;
        }
        static double divis[] = {2.0, 1.4, 1.2, 1.1};
        for (int i = 0; i < 4; i++) {
            u = utx / divis[i];
            if (truncation(u, 0.0) <= accx) utx = u;
        }
    }

    void QuadFormClass::integrate(int nterm, double interv, double tausq, bool mainx) {
        double inpi = interv / M_PI;
        for (int k = nterm; k >= 0; k--) {
            double u = (k + 0.5) * interv;
            double sum1 = -2.0 * u * c;
            double sum2 = fabs(sum1);
            double sum3 = -0.5 * sigsq * pow(u, 2);
            for (int j = r - 1; j >= 0; j--) {
                int nj = n[j];
                double x = 2.0 * lb[j] * u;
                double y = pow(x, 2);
                sum3 -= 0.25 * nj * log1(y, TRUE);
                y = nc[j] * x / (1.0 + y);
                double z = nj * atan(x) + y;
                sum1 += z;
                sum2 += fabs(z);
                sum3 -= 0.5 * x * y;
            }
            double x = inpi * exp1(sum3) / u;
            if (!mainx) x *= (1.0 - exp1(-0.5 * tausq * pow(u, 2)));
            sum1 = sin(0.5 * sum1) * x;
            sum2 = 0.5 * sum2 * x;
            intl += sum1;
            ersm += sum2;
        }
    }

    double QuadFormClass::cfe(double x) {
    /* Coefficient of tausq in error when convergence factor of
       exp1(-0.5*tausq*u^2) is used when df is evaluated at x */
    double axl, axl1, axl2, sxl, sum1, lj;
    int j, k, t;
    axl = fabs(x);
    sxl = (x > 0.0) ? 1.0 : -1.0;
    sum1 = 0.0;
    for (j = r - 1; j >= 0; j--) {
        t = th[j];
        if (lb[t] * sxl > 0.0) {
            lj = fabs(lb[t]);
            axl1 = axl - lj * (n[t] + nc[t]);
            axl2 = lj / log(28.0);
            if (axl1 > axl2) {
                axl = axl1;
            } else {
                if (axl > axl2) axl = axl2;
                sum1 = (axl - axl1) / lj;
                for (k = j - 1; k >= 0; k--)
                    sum1 = sum1 + (n[th[k]] + nc[th[k]]);
                goto l;
            }
        }
    }
l:
    if (sum1 > 100.0) {
        fail = true;
        return 1.0;
    } else {
        return pow(2.0, (sum1 / 4.0)) / (M_PI * pow(axl,2));
    }
}

  void QuadFormClass::qfc_1(
	arma::vec& lb1, arma::vec& nc1, arma::ivec& n1, int r,
           double  sigma, double c1, int lim,
           double acc, arma::vec& trace, int ifault,
           double & res
) {
    int count = 0;
    //int r = r1[0], lim = lim1[0], count = 0;
    double qfval = -1.0, sd, lmax = 0.0, lmin = 0.0, mean = 0.0;
    double xlim = (double)lim, utx, tausq, almx, xnt, xntm, intv, intv1, x, up, un;
    double d1, d2, lj, ncj;

    // Initialize trace
    trace.fill(0.0);
    //ifault.fill(0);
    ifault = 0;

     arma::ivec rats={1,2,4,8};

    // Convert input arma::vec to std::vector for compatibility
    //std::vector<int> n(n1.begin(), n1.end());
    //std::vector<double> lb(lb1.begin(), lb1.end());
    //std::vector<double> nc(nc1.begin(), nc1.end());

    // Parameter validation
    sigsq = std::pow(sigma,2);
    sd = sigsq;

    int nj;
    //double lj, ncj;
    for (int j = 0; j < r; ++j) {
        nj = n1[j];
        lj = lb1[j];
        ncj = nc1[j];
        if (nj < 0 || ncj < 0.0) {
            ifault = 3; // Invalid parameters
            return;
        }
        sd += std::pow(lj,2) * (2 * nj + 4.0 * ncj);
        mean += lj * (nj + ncj);
        lmax = std::max(lmax, lj);
        lmin = std::min(lmin, lj);
    }

    if (sd == 0.0) {
        qfval = (c > 0.0) ? 1.0 : 0.0;
        res = qfval;
        return;
    }

    if (lmin == 0.0 && lmax == 0.0 && sigma == 0.0) {
        ifault = 3; // Invalid parameters
        return;
    }

    sd = std::sqrt(sd);
    almx = (lmax < -lmin) ? -lmin : lmax;

    // Starting values for findu, ctff
    utx = 16.0 / sd;
    up = 4.5 / sd;
    un = -up;

    // Truncation point with no convergence factor
    findu(utx, 0.5 * acc);

    if (c != 0.0 && (almx > 0.07 * sd)) {
        tausq = 0.25 * acc / cfe(c);
        if (false) { /* Fail condition placeholder */ }
        else if (truncation(utx, tausq) < 0.2 * acc) {
            sigsq += tausq;
            findu(utx, 0.25 * acc);
            trace[5] = std::sqrt(tausq);
        }
    }
    trace[4] = utx;
    acc = 0.5 * acc;
    // Find range of distribution
    while (true) {
        d1 = ctff(acc, up) - c;
        if (d1 < 0.0) { qfval = 1.0; break; }
        d2 = c - ctff(acc, un);
        if (d2 < 0.0) { qfval = 0.0; break; }

        intv = 2.0 * M_PI / std::max(d1, d2);
        xnt = utx / intv;
        xntm = 3.0 / std::sqrt(acc);

        if (xnt > xntm * 1.5) {
            if (xntm > xlim) {
                ifault = 1; // Accuracy not achieved
                break;
            }
            int ntm = std::floor(xntm + 0.5);
            intv1 = utx / ntm;
            x = 2.0 * M_PI / intv1;
            if (x <= std::abs(c)) break;

            tausq = 0.33 * acc / (1.1 * (cfe(c - x) + cfe(c + x)));
            if (false) break; // Placeholder fail condition

            acc = 0.67 * acc;
            integrate(ntm, intv1, tausq, false);
            xlim -= xntm;
            sigsq += tausq;
            trace[2] += 1;
            trace[1] += ntm + 1;

            findu(utx, 0.25 * acc);
            acc = 0.75 * acc;
            continue;
        }

        trace[3] = intv;
        if (xnt > xlim) { ifault = 1; break; }

        int nt = std::floor(xnt + 0.5);
        integrate(nt, intv, 0.0, true);
        trace[2] += 1;
        trace[1] += nt + 1;

        qfval = 0.5 - intl;
        trace[0] = ersm;

        up = ersm;
        x = up + acc / 10.0;
        for (int j = 0; j < 4; ++j) {
            if (rats[j] * x == rats[j] * up) ifault = 2; // Round-off error
        }
        break;
    }

    res = qfval;


}    



//};

}
