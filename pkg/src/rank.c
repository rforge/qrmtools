/* C function for computing the ranks faster than rank() **********************/

#include "rank.h"

SEXP rank_aux_(SEXP x)
{
	/* Setup */
	double *x_ = REAL(x); /* pointer to x */
	int *N = LENGTH(x); /* define N */
	SEXP res=PROTECT(allocVector(REALSXP, N)); /* result containing the ranks */

	/* Compute rank() */
	y = (double *) R_alloc(N, sizeof(double)); /* N-vector */
	R_orderVector(y, N, Rf_lang1(x), TRUE, /* nalast (use same default as order()) */
		      FALSE); /* decreasing */
	R_orderVector(res, N, Rf_lang1(y), TRUE, /* nalast (use same default as order()) */
		      FALSE); /* decreasing */

	/* Return */
	UNPROTECT(1);
	return(res);
}