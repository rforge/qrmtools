/* C functions for computing VaR bounds ***************************************/

#include "VaR_bounds.h"

/**
 * @title Determine the indices which order any increasing (!) vector y
 *        oppositely to the given vector x
 * @param x N-vector (one column in the (A)RA()/rearrange() input matrix X)
 * @return order(order(x, decreasing=TRUE)) = N+1-rank(x) (since
 *         order(order(x)) = rank(x))
 * @author Marius Hofert, Kurt Hornik
 */
SEXP indices_opp_ordered_to(SEXP x)
{
    int N = LENGTH(x); /* define N */
    int i; /* running index */

    double *x_; /* define pointer to x */
    x_ = REAL(x); /* set pointer to x */

    SEXP ind = PROTECT(allocVector(INTSXP, N)); /* vector of indices */
    int *ind_; /* define pointer to ind */
    ind_ = INTEGER(ind); /* set pointer to ind */
    for(i=0; i<N; i++) ind_[i] = i; /* init */

    SEXP res = PROTECT(allocVector(INTSXP, N)); /* result vector */
    int *res_; /* define pointer to res */
    res_ = INTEGER(res); /* set pointer to res */

    SEXP xx = PROTECT(allocVector(REALSXP, N)); /* copy of x */
    double *xx_; /* define pointer to xx */
    xx_ = REAL(xx); /* set pointer to xx */

    for(i=0; i<N; i++) xx_[i] = x_[i]; /* as x is changed otherwise */
    rsort_with_index(xx_, ind_, N); /* => ind_ contains the perm. which sorts x_, i.e., order(x_) */

    for(i=0; i<N; i++) res_[ind_[i]] = i; /* res = rank(x)-1 */
    for(i=0; i<N; i++) res_[i] = N-res_[i]; /* compute N+1-rank(x) */
    UNPROTECT(3);
    return(res); /* return */

    /* Via R_orderVector() (too slow, should use R_orderVector1(), then maybe faster): */

    /* /\* Setup *\/ */
    /* int N = LENGTH(x); /\* define N *\/ */
    /* SEXP res = PROTECT(allocVector(INTSXP, N)); */
    /* SEXP ind = PROTECT(allocVector(INTSXP, N)); */
    /* int i, *y; */

    /* /\* Compute order(order(, decreasing=TRUE)) *\/ */
    /* y = INTEGER(ind); */
    /* R_orderVector(y, /\* result *\/ */
    /* 	      N, /\* length *\/ */
    /* 	      Rf_lang1(x), /\* argument *\/ */
    /* 	      TRUE, /\* nalast (use same default as order()) *\/ */
    /* 	      TRUE); /\* decreasing TRUE *\/ */
    /* y = INTEGER(res); */
    /* R_orderVector(y, /\* result *\/ */
    /* 	      N, /\* length *\/ */
    /* 	      Rf_lang1(ind), /\* argument *\/ */
    /* 	      TRUE, /\* nalast (use same default as order()) *\/ */
    /* 	      FALSE); /\* decreasing FALSE *\/ */
    /* for(i = 0; i < N; i++) y[i] += 1; /\* increase all by 1 *\/ */

    /* /\* Return *\/ */
    /* UNPROTECT(2); */
    /* return(res); */
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
