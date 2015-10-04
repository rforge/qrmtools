// C++ functions for computing Steps 4 and 5 of the RA /////////////////////////

#include <Rcpp.h>

using namespace Rcpp;

// order()
// See http://stackoverflow.com/questions/21609934/ordering-permutation-in-rcpp-i-e-baseorder
IntegerVector order_(NumericVector x)
{
  if(is_true(any(duplicated(x)))) {
    Rf_warning("There are duplicates in 'x'; order not guaranteed to match that of R's base::order");
  }
  NumericVector sorted = clone(x).sort();
  return match(sorted, x);
}

// rev()
// See http://stackoverflow.com/questions/23073270/how-to-reverse-an-integer-vector-efficiently
IntegerVector rev_(NumericVector x)
{
	return rev(x); // Rcpp sugar; alternatively, use: std::reverse(x.begin(), x.end())
}

// rowSums() without jth column
NumericVector rowSumsMinus_(NumericMatrix x, int j)
{
	int i, k;
	NumericVector res;
	for(j=0; j<d; j++) { /* walk over all columns */

	}
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
    for(j=0; j<d; j++) { /* walk over all columns */
	/* Compute row sums over all other columns except jth */
        for(k=0; k<N; k++) {
            rs[k] = 0.0;
            for(l=0; l<d; l++) { if(l!=j) rs[k] += X[k][l]; }
        }
        /* Determine number of columns which are oppositely order to the sum of all others */
        col_opp_ordered = TRUE;
        for(i=0; i<N; i++) {
            for(k=0; k<N; k++) {
		    /* printf("imr: %f\n", X[i][j]-X[k][j]); TODO */
		    /* printf("imr2: %f\n", rs[i]-rs[k]); TODO */
	        if((X[i][j]-X[k][j])*(rs[i]-rs[k]) > 0) {
		    col_opp_ordered = FALSE;
		    goto end; /* appropriate use of goto */
	        }
            }
        }
        end:
        if(col_opp_ordered == TRUE) num_cols_opp_ordered++;
    }
    return num_cols_opp_ordered;
}

/**
 * @title Auxiliary Function for Computing Steps 4 and 5 of the RA
 * @param X (N, d)-matrix (either \underline{X}^\alpha or \overline{X}^\alpha)
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
 * @param m_row_sums_xtnd_size
 * @param individual_err (individual) error reached
 * @param m_row_sums_xtnd minimal [for worst VaR] or maximal [for best VaR] row sums
 *        for each iteration
 * @param num_opp_ordered number of oppositely ordered columns
 * @param count number of iterations through the matrix columns
 * @return void
 * @author Marius Hofert
 */
void RA_aux(double **X, int N, int d, const char *method, const char *err,
            int maxiter, double eps, int m_row_sums_xtnd_size,
            double *individual_err, double *m_row_sums_xtnd,
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
    Y = (double **) R_alloc(N, sizeof(double)); /* (N, d) matrix */
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
	/* Initialize Y (to be X) */
        for(i=0; i<N; i++) {
            for(j=0; j<d; j++) Y[i][j] = X[i][j];
        }
        /* For all j in {1,..,d}, oppositely reorder Y[,j] w.r.t. to the sum of all others */
        for(j=0; j<d; j++) {
            /* Compute the row sum over all columns except jth */
            for(i=0; i<N; i++) {
                rs_[i] = 0.0;
                for(l=0; l<d; l++) if(l!=j) rs_[i] += Y[i][l];
            }
            /* Oppositely order Y[,j] with respect to rs_ */
            for(i=0; i<N; i++) Y_j[i] = Y[i][j]; /* pick out jth column of Y */
            /* Determine rank(rs) [= order(rs)] */
	    R_orderVector(ind, N, Rf_lang1(rs), TRUE, /* nalast (use same default as order()) */
			  TRUE); /* decreasing */
            /* Would like to order Y_j in decreasing order but go with increasing order
               and then appropriately subset in the following line */
	    R_rsort(Y_j, N); /* R's sort() for real arguments */
	    for(i=0; i<N; i++) Y[i][j] = Y_j[N-1-ind[i]]; /* update jth column of Y */
        }
	/* printf("num opp ordered = %d", num_opp_order_cols(Y, N, d)); TODO: should be at least 1 at this stage! */
        /* Check whether m_row_sums_xtnd has space left */
        if(cnt >= m_row_sums_xtnd_size) {
		m_row_sums_xtnd_size += 64; /* enlarge size */
		m_row_sums_xtnd = (double *) S_realloc((char *) m_row_sums_xtnd,
                                                       m_row_sums_xtnd_size, /* new size */
                                                       m_row_sums_xtnd_size - 64, /* old size */
                                                       m_row_sums_xtnd_size * sizeof(double));
		/* An alternative would be (does not require Free() here): */
		/* m_row_sums_xtnd = (double *) Realloc(m_row_sums_xtnd, */
                /*                                      m_row_sums_xtnd_size, /\* new size *\/ */
	        /*                                      double); */
        }
        /* Compute and store minimal/maximal row sums */
	for(i=0; i<N; i++) {
            rs_[i] = 0.0;
            for(l=0; l<d; l++) rs_[i] += Y[i][l];
	}
	/* Compute new min/max row sum */
        if(strcmp(method, "worst") == 0) {
		mrs_new = min(rs_, N);
	} else {
		mrs_new = max(rs_, N);
	}
	/* Append it to m_row_sums_xtnd */
        m_row_sums_xtnd[cnt-1] = mrs_new;
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
			    ierr = fabs(mrs_new - mrs_old);
		    } else {
			    ierr = fabs((mrs_new - mrs_old) / mrs_old);
		    }
		    if(ierr <= eps) { stp = TRUE; } else { stp = FALSE; }
	    }
        }
        if(stp) {
	    /* Set counter */
	    *count = cnt;
            /* Count number of oppositely ordered columns */
	    *num_opp_ordered = num_opp_order_cols(Y, N, d);
            /* Compute/set individual error */
	    *individual_err = ierr;
	    break;
  	} else { memcpy(X, Y, N*d * sizeof(double)); } /* destination, source */
    }

    /* clean-up */
    UNPROTECT(1); /* clean-up; number = # of PROTECT() calls */
}


// @title R Interface to C++ for Computing Steps 4 and 5 of the RA
// @param X An (N,d)-matrix (either \underline{X}^\alpha or \overline{X}^\alpha)
// @param method A character string indicating which VaR is approximated (worst/best)
//        ("worst" or "best")
// @param err A character string indicating the error function used
//        ("absolute" or "relative")
// @param maxiter The maximal number of iterations; if < 0, then the iteration
//        is done until convergence determined by eps
// @param eps The (individual) epsilon error to determine convergence; if < 0,
//        the iteration is done until the matrix doesn't change anymore
// @return 5-list containing the
//         1) computed (lower or upper [depending on X]) bound for (worst or
//            best [depending on method]) VaR
//         2) (individual) error reached
//         3) minimal [for worst VaR] or maximal [for best VaR] row sums
//            for each iteration
//         4) number of oppositely ordered columns
//         5) number of iterations through the matrix columns
// @author Marius Hofert
RcppExport SEXP rearrange_(SEXP X, SEXP method, SEXP err, SEXP maxiter, SEXP eps)
{
	// Convert the SEXPs to C++ objects
	NumericMatrix X.(X);
	CharacterVector method.(method);
	CharacterVector err.(err);
	int maxiter. = as<int>(maxiter);
	double eps. = as<double>(eps);

	// Auxiliary variables
	int N = X.nrow();
	int d = X.ncol();
	int j;
	NumericVector Y.j;
	NumericMatrix row.sums(N, 0); // (N, 0)-matrix filled with 0s
	NumericVector mrs.old(N); // N-vector filled with 0s
	NumericMatrix Y.(N,d);

	// Loop through the colums
	while(true) {
		memcpy(Y., X., sizeof(Y.)); // destination, source
		for(j=0; j<d; j++) {

			ord.Y.j = order_(Y.(_,j)); // order(jth column)

			// TODO
		}
	}

	double individual_err;
	int num_opp_ordered, count;



	// Return
	List res;
	res["bound"] = ;
	res["err"] = ;
	res["num.iter"] = ;
	res["row.sums"] = ;
	res["m.row.sums"] = ;
	res["num.opp.ordered"] = ;
	return res;
}



	if(method_ == "worst") {
		double optim.fun()
{
	int i;
	double res = x[0];
	for(i=1; i<n; i++) if(x[i] < res) res = x[i];
	return res;
}

	} else {

	}
	optim.fun <- if(method=="worst") min else max
							  err.fun <- if(err.type=="absolute") {
              function(x, y) abs(x-y)
          } else {
              function(x, y) abs((x-y)/y)
          }
 	// Loop through the columns

	row.sums <- matrix(, nrow=nrow(X), ncol=0) # (N, 0)-matrix of row sums
		mrs.old <- optim.fun(rowSums(X)) # old minimal row sum (currently: the one of X)
		while (TRUE) {
			// Oppositely order X (=> Y)
			Y <- X
				for(j in 1:d)
					Y[,j] <- sort(Y[,j], decreasing=TRUE)[rank(rowSums(Y[,-j, drop=FALSE]))]
						// Compute row sums and minimal/maximal row sums
						Y.rs <- rowSums(Y)
						row.sums <- cbind(row.sums, Y.rs) # append the new row sums
						mrs.new <- optim.fun(Y.rs) # compute new minimal row sum
						// Check convergence (we use "<= eps" as it entails eps=0)
						stp <- (ncol(row.sums) == maxiter) || if(is.null(eps)) { all(Y == X) } else {
			err.fun(mrs.new, mrs.old) <= eps
				}
			if(stp) {
				num.opp.ordered <- sum(sapply(seq_len(d), function(j)
							      all(sort(Y[,j], decreasing=TRUE)[rank(rowSums(Y[,-j, drop=FALSE]))] == Y[,j]))) # count number of oppositely ordered columns
					err <- err.fun(mrs.new, mrs.old) # compute the (individual) error
					break
					} else {
				mrs.old <- mrs.new # update mrs.old
					X <- Y # update X
					}
}
	// Return
	colnames(row.sums) <- NULL # remove column names so that they don't appear in output
          list(bound=mrs.new, err=err, num.iter=ncol(row.sums),
               row.sums=row.sums, m.row.sums=apply(row.sums, 2, optim.fun),
               num.opp.ordered=num.opp.ordered)


}

