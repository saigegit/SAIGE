#ifndef QFC_HPP
#define QFC_HPP


#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <cmath>



#define UseDouble 0             /* all floating point double */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/*#include <malloc.h>*/
#include <setjmp.h>
#include "qfc.hpp"


#define TRUE  1
#define FALSE 0
typedef int BOOL;


#define pinew 3.14159265358979
#define log28 .0866  /*  log(2.0) / 8.0  */

//extern "C" {

static double sigsq, lmax, lmin, mean, c;
static double intl, ersm;
static int count0, r, lim;  static BOOL ndtsrt, fail;
static int *n,*th; static double *lb,*nc;
static jmp_buf env;




static double exp1(double x);
static void counter(void);
static double square(double x);
static double cube(double x);
static double  log1(double x, BOOL first);
static void order(void);
static double   errbd(double u, double* cx);
static double  ctff(double accx, double* upn);
static double truncation(double u, double tausq);
static void findu(double* utx, double accx);
static void integrate(int nterm, double interv, double tausq, BOOL mainx);
static double cfe(double x);

void  qfc_1(double* lb1, double* nc1, int* n1, int *r1, double *sigma, double *c1, int *lim1, double *acc, double* trace, int* ifault, double *res);

#endif
