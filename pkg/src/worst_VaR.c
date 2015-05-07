/* C functions for computing Steps 4 and 5 of the RA **************************/

#include "worst_VaR.h"


/**
 * @title Minimum
 * @param x vector x
 * @param n length of x
 * @return min(x)
 * @author Marius Hofert
 */
double min(double *x, int n)
{
	int i;
	double res = x[0];
	for(i=1; i<n; i++) if(x[i] < res) res = x[i];
	return res;
}

/**
 * @title Maximum
 * @param x vector x
 * @param n length of x
 * @return max(x)
 * @author Marius Hofert
 */
double max(double *x, int n)
{
	int i;
	double res = x[0];
	for(i=1; i<n; i++) if(x[i] > res) res = x[i];
	return res;
}

/**
 * @title Number of Columns of an (N,d)-Matrix which are Oppositely
 *        Ordered to the Sum of All Other Columns
 * @param X (N,d)-matrix
 * @param N number of rows of X
 * @param d number of columns of X
 * @return number of columns oppositely orderd to the sum of all other columns
 * @author Marius Hofert
 */
int num_opp_order_cols(double **X, int N, int d)
{
    double *rs; /* for *r*ow *s*ums (without jth column) */
    rs = (double *) R_alloc(N, sizeof(double)); /* (double *) R_alloc(n, sizeof(double)) */
    Rboolean col_opp_ordered; /* current col opp. ordered to the sum of all others? */
    int num_cols_opp_ordered = 0;
    int j, k, l, i;
    for(int j=0; j<d; j++) { /* walk over all columns */
	/* Compute row sums over all other columns except jth */
        for(k=0; k<N; k++) {
            rs[k] = 0.0;
            for(l=0; l<d; l++) { if(l!=j) rs[k] += X[k][l]; }
        }
        /* Determine number of columns which are oppositely order to the sum of all others */
        col_opp_ordered = TRUE;
        for(i=0; i<N; i++) {
            for(k=0; k<N; k++) {
                if((X[i][j]-X[k][j])*(rs[i]-rs[k]) > 0) {
                    col_opp_ordered = FALSE;
                    goto end; /* appropriate use of goto */
                }
            }
        }
        end:
        if(col_opp_ordered == TRUE) num_cols_opp_ordered += 1;
    }
    return num_cols_opp_ordered;
}

/**
 * @title Auxiliary Function for Computing Steps 4 and 5 of the RA
 * @param X (N,d)-matrix (either \underline{X}^\alpha or \overline{X}^\alpha)
 * @param N nrow(X)
 * @param d ncol(X)
 * @param method character indicating which VaR is approximated (worst/best)
 *        ("worst" or "best")
 * @param err character string indicating the error function used
 *        ("absolute" or "relative")
 * @param maxiter maximal number of iterations; if < 0, then the iteration
 *        is done until convergence determined by eps
 * @param eps epsilon error to determine convergence; if < 0, then the
 *        iteration is done until the matrix doesn't change anymore
 * @param m_row_sums_size
 * @param individual_err (individual) error reached
 * @param m_row_sums minimal [for worst VaR] or maximal [for best VaR] row sums
 *        for each iteration
 * @param num_opp_ordered number of oppositely ordered columns
 * @param count number of iterations through the matrix columns
 * @return void
 * @author Marius Hofert
 */
void RA_aux(double **X, int N, int d, const char *method, const char *err,
            int maxiter, double eps, int m_row_sums_size,
            double *individual_err, double *m_row_sums,
            int *num_opp_ordered, int *count)
{
    /* Running variables */
    int i, j, l;

    /* Allocate memory */
    int cnt;
    double mrs_old, mrs_new, ierr;
    SEXP rs = PROTECT(allocVector(REALSXP, N)); /* vector of row sums; needs to be SEXP for R_orderVector() */
    double *rs_ = REAL(rs); /* pointer to the value of rs */
    double *Y_j;
    Y_j = (double *) R_alloc(N, sizeof(double)); /* Y[,j] */
    int *ind; /* for order (permutation of 0:(N-1)) computed by R_orderVector() */
    ind = (int *) R_alloc(N, sizeof(int));
    double **Y; /* matrix; for new iteration of oppositely ordering */
    Y = (double **) R_alloc(N, sizeof(double));
    for(i=0; i<N; i++) Y[i] = (double *) R_alloc(d, sizeof(double));
    Rboolean stp, change;

    /* Loop */
    cnt = 0;
    for(;;) {

	/* Counter related quantities */
        cnt++; /* increase counter */
        if(cnt == 1) {
            for(i=0; i<N; i++) {
                rs_[i] = 0.0;
                for(j=0; j<d; j++) rs_[i] += X[i][j];
            }
	    /* compute min/max row sum */
	    if(strcmp(method, "worst") == 0) {
		    mrs_old = min(rs_, N);
	    } else {
		    mrs_old = max(rs_, N);
	    }
        } else { mrs_old = mrs_new; } /* old min/max row sum */

        /* Go through all columns of X and oppositely order them w.r.t. to the sum of all others */
        for(j=0; j<d; j++) {
            /* Compute the row sum over all columns except jth */
            for(i=0; i<N; i++) {
                rs_[i] = 0.0;
                for(l=0; l<d; l++) if(l!=j) rs_[i] += Y[i][l];
            }
            /* Oppositely order Y[,j] with respect to rs_; so sort(Y[,j])[rev(rank(rs_))] */
            for(i=0; i<N; i++) Y_j[i] = Y[i][j]; /* pick out jth column of Y */
	    printf("Hello 0\n");
            R_orderVector(ind, N, rs, TRUE, /* nalast (use same default as order()) */
			  TRUE); /* decreasing */
	    /* => ind == rev(rank(rs_)) */
	    printf("Hello 1\n");
	    R_rsort(Y_j, N); /* R's sort() for real arguments */
	    printf("Hello 2\n");
	    for(i=0; i<N; i++) Y[i][j] = Y_j[ind[i]];  /* update jth column of Y */
        }
	    printf("Hello 3\n");
        /* Check whether m_row_sums has space left */
        if(cnt >= m_row_sums_size) {
		m_row_sums = S_realloc(m_row_sums, m_row_sums_size,
				       m_row_sums_size + 64,
				       m_row_sums_size * sizeof(double));
		m_row_sums_size += 64;
        }
	    printf("Hello 4\n");
        /* Compute and store minimal/maximal row sums */
	for(i=0; i<N; i++) {
            rs_[i] = 0.0;
            for(l=0; l<d; l++) rs_[i] += Y[i][l];
	}
	/* compute new min/max row sum */
        if(strcmp(method, "worst") == 0) {
		mrs_new = min(rs_, N);
	} else {
		mrs_new = max(rs_, N);
	}
	/* append it to m_row_sums */
        m_row_sums[cnt-1] = mrs_new;
	    printf("Hello 5\n");
        /* Check convergence (we use "<= eps" as it entails eps=0) */
        if(cnt == maxiter) stp = TRUE; else { /* check number of iterations */
            if(eps < 0) { /* check whether there was no change in the matrix */
                change = FALSE;
                for(i=0; i<N; i++) {
                    for(j=0; j<d; j++) {
                        if(Y[i][j]!=X[i][j]) {
                            change = TRUE;
                            goto end; /* appropriate use of goto */
                        }
                    }
                }
                end:
                if(change) { stp = FALSE; } else { stp = TRUE; }
             }  else { /* check individual convergence criterion */
		    if(strcmp(err, "absolute") == 0) {
			    ierr = abs(mrs_new - mrs_old);
		    } else {
			    ierr = abs((mrs_new - mrs_old) / mrs_old);
		    }
		    if(ierr <= eps) { stp = TRUE; } else { stp = FALSE; }
	    }
        }
		    printf("Hello 6\n");
        if(stp) {
	    /* Set counter */
	    (*count) = cnt;
            /* Count number of oppositely ordered columns */
	    (*num_opp_ordered) = num_opp_order_cols(Y, N, d);
            /* Compute/set individual error */
	    (*individual_err) = ierr;
	    break;
  	} else { memcpy(X, Y, N*d * sizeof(double)); } /* destination, source */
    }

    /* clean-up */
    UNPROTECT(1); /* clean-up; number = # of PROTECT() calls */
}

/**
 * @title R Interface to C for Computing Steps 4 and 5 of the RA
 * @param X (N,d)-matrix (either \underline{X}^\alpha or \overline{X}^\alpha)
 * @param method character indicating which VaR is approximated (worst/best)
 *        ("worst" or "best")
 * @param err character string indicating the error function used
 *        ("absolute" or "relative")
 * @param maxiter maximal number of iterations; if < 0, then the iteration
 *        is done until convergence determined by eps
 * @param eps epsilon error to determine convergence; if < 0, then the
 *        iteration is done until the matrix doesn't change anymore
 * @return 5-list containing the
 *         1) computed (lower or upper [depending on X]) bound for (worst or
 *            best [depending on method]) VaR
 *         2) (individual) error reached
 *         3) minimal [for worst VaR] or maximal [for best VaR] row sums
 *            for each iteration
 *         4) number of oppositely ordered columns
 *         5) number of iterations through the matrix columns
 * @author Marius Hofert
 */
SEXP RA_aux_(SEXP X, SEXP method, SEXP err, SEXP maxiter, SEXP eps)
{
    /* Input parameters */
    double *X_          = REAL(X); /* (N,d)-matrix encoded as N*d vector; returns a pointer */
    const char *method_ = CHAR(STRING_ELT(method, 0)); /* character(1); CHAR() returns a pointer, has to be const */
    const char *err_    = CHAR(STRING_ELT(err, 0)); /* character(1) */
    int maxiter_        = asInteger(maxiter); /* numeric(1); convert double to integer */
    double *eps_        = REAL(eps); /* numeric(1); returns a pointer */
    int *dim            = INTEGER(getAttrib(X, R_DimSymbol)); /* dim(X) */
    int N = dim[0]; /* nrow(X) */
    int d = dim[1]; /* ncol(X) */

    /* Define auxiliary variables */
    int i, j;
    int m_row_sums_size = 64; /* keep track of the length of m_row_sums */

    /* Allocate matrix and copy in X */
    double **X__;
    X__ = (double **) R_alloc(N, sizeof(double));
    for(i=0; i<N; i++) X__[i] = (double *) R_alloc(d, sizeof(double));
    for(i=0; i<N; i++) for(j=0; j<d; j++) X__[i][j] = X_[j + i*d];

    /* Allocate memory for output objects and construct pointers to them
       Kurt: - For scalars, use "double *foo;" and then "ScalarReal(*foo)" (= SEXP)
             - For vectors, use "SEXP foo = PROTECT(allocVector(REALSXP, <size>))"
               and then define a pointer to it via "double *foo_ = REAL(foo)" */
    double *individual_err;
    SEXP m_row_sums = PROTECT(allocVector(REALSXP, m_row_sums_size)); /* length 64 (expanded if required) */
    double *m_row_sums_ = REAL(m_row_sums); /* pointer to the value of m_row_sums */
    int *num_opp_ordered;
    int *count;

    /* Main */
    RA_aux(X__, N, d, method_, err_, maxiter_, (*eps_), /* inputs */
	   m_row_sums_size, /* auxiliary */
	   individual_err, m_row_sums_, num_opp_ordered, count); /* outputs */
    double bound = m_row_sums_[(*count)-1]; /* extract bound */

    /* Create result object */
    SEXP res = PROTECT(allocVector(VECSXP, 5)); /* list of length 5 */
    SET_VECTOR_ELT(res, 0, ScalarReal(bound)); /* computed bound (min/max row sum); numeric(1) */
    SET_VECTOR_ELT(res, 1, ScalarReal(*individual_err)); /* individual error; numeric(1) */
    SET_VECTOR_ELT(res, 2, m_row_sums); /* min/max row sums; numeric(<some length>) */
    SET_VECTOR_ELT(res, 3, ScalarInteger(*num_opp_ordered)); /* number of oppositely ordered cols; integer(1) */
    SET_VECTOR_ELT(res, 4, ScalarInteger(*count)); /* number of iterations; integer(1) */

    /* Name sublists */
    char *names[5] = {"bound", "individual.err", "m.row.sums",
		      "num.opp.ordered", "num.iter"};
    for(i = 0; i < 5; i++) SET_STRING_ELT(res, i, mkString(names[i]));
    SEXP nms = PROTECT(allocVector(STRSXP, 5)); /* names as SEXP */
    setAttrib(res, R_NamesSymbol, nms);

    /* Return */
    UNPROTECT(3); /* clean-up; number = # of PROTECT() calls */
    return res;
}
