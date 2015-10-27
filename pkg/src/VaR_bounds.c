/* C functions for computing VaR bounds ***************************************/

#include "VaR_bounds.h"

/* Rank pairs struct */
struct rank_pair {
    double val;
    size_t ind;
};

/* Comparison function for rank pairs */
int cmp_rank_pair(const void* a, const void* b) {
    struct rank_pair *lhs = (struct rank_pair*)a;
    struct rank_pair *rhs = (struct rank_pair*)b;
    return lhs->val < rhs->val ? -1 : (lhs->val > rhs->val ? 1 : 0);
}

/* Compute rank(a) and store result in r */
/* See http://stackoverflow.com/questions/33358521/how-to-compute-the-rank-of-a-vector-in-c/33358729#33358729 */
void rank(double a[], int r[], size_t n) {
    struct rank_pair *tmp = R_alloc(n, sizeof(struct rank_pair));
    for (int i = 0 ; i != n ; i++) {
    	tmp[i].val = a[i];
    	tmp[i].ind = i;
    }
    qsort(tmp, n, sizeof(struct rank_pair), cmp_rank_pair);
    for (int i = 0 ; i != n ; i++) {
    	r[tmp[i].ind] = i+1;
    }
}

/**
 * @title Determine the indices which order any increasing (!) vector y
 *        oppositely to the given vector x
 * @param x N-vector (one column in the (A)RA()/rearrange() input matrix X)
 * @return order(order(, decreasing=TRUE))
 * @author Marius Hofert, Kurt Hornik
 * @note What we are looking for is N+1-rank(x) which is also equal to
 *       or order(order(x, decreasing=TRUE)) since rank(x) = order(order(x))
 */
SEXP indices_opp_ordered_to(SEXP x)
{
    int N = LENGTH(x); /* define N */
    double *x_; /* define pointer to x */
    x_ = REAL(x); /* set pointer to x */
    int i; /* running index */
    SEXP res = PROTECT(allocVector(INTSXP, N)); /* result vector */
    int *res_; /* define pointer to res */
    res_ = INTEGER(res); /* set pointer to res */
    rank(x_, res_, N); /* compute rank(x) */
    for(i=0; i<N; i++) res_[i] = N+1-res_[i]; /* compute N+1-rank(x) */
    UNPROTECT(1);
    return(res); /* return */

    /* Via R_orderVector() (too slow, should use R_orderVector1()): */

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
