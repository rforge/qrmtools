/* C header for worst_VaR.c ***************************************************/

#ifndef worst_VaR_H
#define worst_VaR_H

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

/* Auxiliary tools */
double min(double *x, int n); /* vectorized min(x) */
double max(double *x, int n); /* vectorized max(x) */
double *row_sums(double **x, int n, int d, int j); /* rowSums(x) */

/* For RA and extensions specifically */
int num_opp_order_cols(double **X, int N, int d); /* count number of oppositely ordered cols in X */
double *opp_order(double *x, double *y, int n); /* oppositely order x w.r.t. y */
double **opp_order_mat(double **X, int N, int d); /* iterate opp_order() over all columns of a matrix */
SEXP RA_aux_C(SEXP X, SEXP method, SEXP err, SEXP maxiter, SEXP eps); /* Steps 4 and 5 for the RA */

#endif

