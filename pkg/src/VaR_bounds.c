/* C functions for computing VaR bounds ***************************************/

#include "VaR_bounds.h"

SEXP rank_(SEXP x)
{
    /* Setup */
    // double *x_ = REAL(x); /* pointer to x */
    int N = LENGTH(x);				/* define N */
    SEXP res = PROTECT(allocVector(INTSXP, N)); /* result containing the ranks */
    SEXP ind = PROTECT(allocVector(INTSXP, N));
    int i, *y;

    /* compute rank() */
    y = INTEGER(ind);
    R_orderVector(y, N, Rf_lang1(x),
		  TRUE,	      /* nalast (use same default as order()) */
		  FALSE);     /* decreasing */
    y = INTEGER(res);
    R_orderVector(y, N, Rf_lang1(ind),
		  TRUE,	      /* nalast (use same default as order()) */
		  FALSE);     /* decreasing */

    for(i = 0; i < N; i++) {
	y[i] += 1;
    }

    /* Return */
    UNPROTECT(2);
    return(res);
}

SEXP colsplit_(SEXP x) 
{
    int *dims = INTEGER(getAttrib(x, R_DimSymbol));
    int N = dims[0], d = dims[1];
    SEXP res = PROTECT(allocVector(VECSXP, d));
    int i = 0, j, k;
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

    UNPROTECT(1);
    return(res);
}
