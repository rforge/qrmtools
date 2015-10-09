/* C functions for computing VaR bounds ***************************************/

#include "VaR_bounds.h"

/**
 * @title Determine the indices which order any increasing (!) vector y
 *        oppositely to the given vector x
 * @param x N-vector (one column in the (A)RA()/rearrange() input matrix X)
 * @return order(order(, decreasing=TRUE))
 * @author Marius Hofert, Kurt Hornik
 * @note y <- c(5, 7, 6, 8) => sort(y, decreasing=TRUE)[rank(y)]
 *       = sort(y, decreasing=TRUE)[order(order(y))]
 *       = sort(y)[order(order(y, decreasing=TRUE))]
 */
SEXP indices_opp_ordered_to(SEXP x)
{
	/* Setup */
	int N = LENGTH(x); /* define N */
	SEXP res = PROTECT(allocVector(INTSXP, N));
	SEXP ind = PROTECT(allocVector(INTSXP, N));
	int i, *y;

	/* Compute rank() */
	y = INTEGER(ind);
	R_orderVector(y, /* result */
		      N,
		      Rf_lang1(x), /* argument */
		      TRUE, /* nalast (use same default as order()) */
		      TRUE); /* decreasing TRUE */
	y = INTEGER(res);
	R_orderVector(y, N, Rf_lang1(ind),
		      TRUE, /* nalast (use same default as order()) */
		      FALSE); /* decreasing FALSE */
	for(i = 0; i < N; i++) y[i] += 1; /* increase all by 1 */

	/* Return */
	UNPROTECT(2);
	return(res);
}

/**
 * @title Fast C version of split(x, col(x))
 * @param x (N,d)-matrix ((A)RA()/rearrange() input matrix X)
 * @return split(x, col(x))
 * @author Marius Hofert, Kurt Hornik
 */
SEXP col_split(SEXP x)
{
	/* Setup */
	int *dims = INTEGER(getAttrib(x, R_DimSymbol));
	int N = dims[0], d = dims[1];
	SEXP res = PROTECT(allocVector(VECSXP, d));
	int i = 0, j, k; /* i runs globally, j runs over all cols, k runs over all rows */

	/* Distinguish int/real matrices */
	switch (TYPEOF(x)) {
	case INTSXP:
		for(j = 0; j < d; j++) {
			SET_VECTOR_ELT(res, j, allocVector(INTSXP, N));
			int *e = INTEGER(VECTOR_ELT(res, j));
			for(k = 0 ; k < N ; i++, k++) {
				e[k] = INTEGER(x)[i];
			}
		}
		break;
	case REALSXP:
		for(j = 0; j < d; j++) {
			SET_VECTOR_ELT(res, j, allocVector(REALSXP, N));
			double *e = REAL(VECTOR_ELT(res, j));
			for(k = 0 ; k < N ; i++, k++) {
				e[k] = REAL(x)[i];
			}
		}
		break;
	}

	/* Return */
	UNPROTECT(1);
	return(res);
}
