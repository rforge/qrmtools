/* C functions for computing the worst VaR_alpha for given margins ************/

#include "worst_VaR.h"


/**
 * Determine the row sums of a given matrix, excluding the j-th column
 * note: if j < 0, then this is rowSums()
 *
 * @param x pointer to an (n,d)-matrix
 * @param n number of rows of x
 * @param d number of columns of x
 * @param j index of the column not used for building row sums
 * @return rowSums(x)
 * @author Marius Hofert
 */
double *row_sums(double **x, int n, int d, int j)
{
    /* static double *rs[n]; => not possible as storage size not constant */
    double *rs; /* for *r*ow *s*ums (without jth column) */
    rs = (double *) malloc(n * sizeof(double)); /* (double *) R_alloc(n, sizeof(double)) */
    for(int i=0; i<n; i++) {
        rs[i] = .0;
        for(int k=0; k<d; k++) if(k!=j) rs[i] += x[i][k];
    }
    return rs;
}

/**
 * Determine the number of columns of an (N,d)-matrix which are oppositely
 * ordered to the sum of all other columns
 *
 * @param X (N,d)-matrix
 * @param N number of rows of X
 * @param d number of columns of X
 * @return number of columns oppositely orderd to the sum of all other columns
 * @author Marius Hofert
 */
int num_opp_order_cols(double **X, int N, int d)
{
    double *rs; /* for *r*ow *s*ums (without jth column) */
    rs = (double *) malloc(N * sizeof(double)); /* (double *) R_alloc(n, sizeof(double)) */
    Rboolean col_opp_ordered; /* current col opp. ordered to the sum of all others? */
    int num_cols_opp_ordered = 0;
    for(int j=0; j<d; j++) { /* walk over all columns */
        rs = row_sums(X, N, d, j); /* row sums over all other columns except jth */
        /* Determine number of columns which are oppositely order to the sum of all others */
        col_opp_ordered = TRUE;
        for(int i=0; i<N; i++) {
            for(int k=0; k<N; k++) {
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
 * Oppositely order an n-dimensional vector x w.r.t. y, i.e., compute
 *
 *     sort(x)[rev(rank(y))]
 *
 * (use case: x = jth col of Y; y = row sum of Y for all but the jth col)
 *
 * @param x n-vector (to be ordered oppositely to y)
 * @param y n-vector (determining order)
 * @param n length of x (and y)
 * @return n-vector oppositely ordered with respect to y
 * @author Marius Hofert
 * Note: - For the C functions used, see Writing R Extensions (2014, Section 6.10)
 *         and ./src/main/sort.c
 *       - If you need your own order(), see http://stackoverflow.com/questions/28790745/how-to-determine-the-ordering-of-elements-in-a-vector/28791220#28791220
 *       - If you need your own sort(), use that x[order(x)] is sorted
 */
double *opp_order(double *x, double *y, int n)
{
    int *indx, *res;
    indx   = (int *) malloc(n * sizeof(int)); /* R's order(); contains the order (permutation of 0:(n-1)) after R_orderVector() */
    res = (double *) malloc(n * sizeof(double));
    R_orderVector(indx, n, Rf_lang1(*x), /* convert x to SEXP */
                  TRUE, /* nalast (use same default as order()) */
                  TRUE); /* decreasing */
    /* => indx == rev(rank(y)) */
    R_rsort(x, n); /* R's sort(x) for real x */
    for(int i=0; i<n; i++) res[i] = x[indx[i]];
    return res;
}

/* TODO: comment */
void RA_aux_c(const double **X_, int N, int d, const char *method_, const char *err_,
              const int maxiter_, const double eps_,
	      double *mrs_new, double *individual_err, double *m_row_sums,
              int *m_row_sums_size, int *num_opp_ordered, int *count)
{
    /* Define auxiliary functions (declare before due to 'conditional definition' */
    double optim_fun(double *x, int n); /* declaration */
    double error_fun(double x, double y); /* declaration */
    if(method_ == "worst") { /* optim_fun = min */
        double optim_fun(double *x, int n) {
            double min = x[0];
	    for(int i=1; i<n; i++) if(x[i] < min) min = x[i];
	    return min;
	}
    } else { /* optim_fun = max */
        double optim_fun(double *x, int n) {
            double max = x[0];
	    for(int i=1; i<n; i++) if(x[i] > max) max = x[i];
	    return max;
	}
    }
    if(err_=="absolute") {
        double err_fun(double x, double y) { abs(x-y); }
    } else {
        double err_fun(double x, double y) { abs((x-y)/y); }
    }

    /* Allocate memory */
    double mrs_old;
    double *rs, *Y_j, *Y_j_oo;
    rs       = (double *) malloc(N * sizeof(double)); /* or (double *) R_alloc(n, sizeof(double)) */
    Y_j      = (double *) malloc(N * sizeof(double)); /* Y[,j] */
    Y_j_oo   = (double *) malloc(N * sizeof(double)); /* oppositely ordered Y[,j] */
    double **Y; /* matrix; for new iteration of oppositely ordering */
    Y        = (double **) malloc(N * sizeof(double));
    for(int i=0; i<N; i++) { Y[i] = (double *) malloc(d * sizeof(double)); }

    Rboolean stp, change;

    /* Loop */
    *count = 0;
    for(;;) {

        /* Counter related quantities */
	(*count)++; /* increase counter */
        if((*count) == 1) {
	    rs = row_sums(X_, N, d, -1);
            mrs_old = optim_fun(rs, N);
        } else { mrs_old = mrs_new; } /* old min/max row sum */

        /* Go through all columns of X and oppositely order them w.r.t. to the sum of all others */
        for(int j=0; j<d; j++) {
            rs = row_sums(Y, N, d, j); /* row sum over all columns except jth */
            for(int i=0; i<N; i++) Y_j[i] = Y[i][j]; /* pick out jth column of Y */
            Y_j_oo = opp_order(Y_j, rs, N); /* oppositely order Y[,j] w.r.t. rs */
            for(int i=0; i<N; i++) Y[i][j] = Y_j_oo[i]; /* update jth column of Y */
        }

        /* Compute minimal/maximal row sums */
        rs = row_sums(Y, N, d, -1);
        mrs_new = optim_fun(rs, N); /* new min/max row sum */

        /* Store minimal/maximal row sums */
        if((*count) >= m_row_sums_size) {
	    (*m_row_sums_size) += 64;
            m_row_sums = realloc(m_row_sums, m_row_sums_size * sizeof(double));
        }
        m_row_sums[(*count)-1] = mrs_new; /* append min/max row sum */

        /* Check convergence (we use "<= eps" as it entails eps=0) */
        if(*count == maxiter_) stp = TRUE; else { /* check number of iterations */
            if(eps_ < 0) { /* check whether there was no change in the matrix */
                change = FALSE;
                for(int i=0; i<N; i++) {
                    for(int j=0; j<d; j++) {
                        if(Y[i][j]!=X[i][j]) {
                            change = TRUE;
                            goto end; /* appropriate use of goto */
                        }
                    }
                }
                end:
                if(change) { stp = FALSE; } else { stp = TRUE; }
             }  else { /* check convergence criterion */
	         if(err_fun(mrs_new, mrs_old) <= eps_) { stp = TRUE; } else { stp = FALSE; }
             }
        }
        if(stp) {
            /* Count number of oppositely ordered columns */
            num_opp_ordered = num_opp_order_cols(Y, N, d);
            /* Compute the (individual) error */
            individual_err = err_fun(mrs_new, mrs_old);
                break;
  	} else { memcpy(X, Y, N*d * sizeof(double)); } /* destination, source */
    }
}

/**
 * Auxiliary Function for Computing Steps 4 and 5 for the RA
 *
 * @param X (N, d)-matrix (either \underline{X}^\alpha or \overline{X}^\alpha)
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
 *         3) minimal [for worst VaR] or maximal [for best VaR] row sum
 *         4) number of oppositely ordered columns
 *         5) number of iterations through the matrix columns used for the last N
 * @author Marius Hofert
 */
SEXP RA_aux_(SEXP X, SEXP method, SEXP err, SEXP maxiter, SEXP eps)
{
    /* Input parameters */
    const double **X_   = PROTECT(REAL(X)); /* (N, d)-matrix */
    const char *method_ = PROTECT(CHAR(STRING_ELT(method, 0))); /* character(1) */
    const char *err_    = PROTECT(CHAR(STRING_ELT(err, 0))); /* character(1) */
    const int maxiter_  = PROTECT(INTEGER(maxiter)); /* integer(1) */
    const double eps_   = PROTECT(REAL(eps)); /* numeric(1) */
    int *dim = INTEGER(GET_DIM(X));
    int N = dim[0]; /* nrow(X) */
    int d = dim[1]; /* ncol(X) */

    /* Allocate memory for output objects */
    double mrs_new, individual_err;
    double m_row_sums = malloc(64 * sizeof(double)); /* length 64 (expanded if required) */
    int m_row_sums_size = 64; /* keep track of its length */
    int *num_opp_ordered, *count;
    num_opp_ordered = malloc(sizeof(int));
    count           = malloc(sizeof(int));

    /* Main */
    RA_aux_c(X_, N, d, method_, err_, maxiter_, eps_, /* inputs */
	     mrs_new, individual_err, m_row_sums, m_row_sums_size, /* outputs */
             num_opp_ordered, count);

    /* Create result object */
    SEXP res = PROTECT(allocVector(VECSXP, 5)); /* list of length 5 */
    SET_VECTOR_ELT(res, 0, ScalarReal(mrs_new)); /* computed bound (min/max row sum); numeric(1) */
    SET_VECTOR_ELT(res, 1, ScalarReal(individual_err)); /* individual error; numeric(1) */
    SET_VECTOR_ELT(res, 2, REAL(m_row_sums)); /* min/max row sums; numeric(???) */
    SET_VECTOR_ELT(res, 3, ScalarInteger(num_opp_ordered)); /* number of oppositely ordered cols; integer(1) */
    SET_VECTOR_ELT(res, 4, ScalarInteger(count)); /* number of iterations; integer(1) */
    const char *names[5] = {"bound", "individual.err", "m.row.sums",
                            "num.opp.ordered", "num.iter"}; /* name the sublists */
    for(int i = 0; i < 5; i++) SET_STRING_ELT(res, i, mkString(names[i]));
    SEXP nms = PROTECT(allocVector(STRSXP, 5)); /* names as SEXP */
    setAttrib(res, R_NamesSymbol, nms);

    /* Return */
    UNPROTECT(9); /* clean-up */
    /* TODO: free memory */
    return res;
}
