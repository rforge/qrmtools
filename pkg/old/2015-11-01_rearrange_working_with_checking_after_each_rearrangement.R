##' @title Basic rearrangement function for (A)RA
##' @param X (N, d)-matrix \underline{X}^\alpha or \overline{X}^\alpha
##' @param tol Tolerance to determine (the individual) convergence;
##'        if NULL, column rearrangements are done until the matrix doesn't
##'        change anymore d consecutive times
##' @param tol.type Character string indicating the tolerance function used
##'        ("relative" or "absolute")
##' @param max.num.col.ra Maximal number of column rearrangements
##' @param method Character indicating which VaR is approximated (worst/best)
##'        determines optimizing function (min for worst VaR; max
##'        for best VaR)
##' @param sample A logical indicating whether each column of the working
##'        matrix is randomly permuted before the rearrangements begin
##' @param is.sorted A logical indicating whether X is columnwise sorted in
##'        increasing order
##' @param trace A logical indicating whether the underlying matrix is
##'        printed after each rearrangement step
##' @return List containing the
##'         1) Computed (lower or upper [depending on X]) bound for (worst or
##'            best [depending on method]) VaR
##'         2) (Individual) tolerance reached
##'         3) Logical indicating whether the algorithm has converged
##'         4) Vector of minimal [for worst VaR] or maximal [for best VaR]
##'            row sums after each considered column rearrangement
##'         5) The (optimally) rearranged (N, d)-matrix
##' @author Marius Hofert and Kurt Hornik
##' @note - We use "<= tol" to determine convergence instead of "< tol" as
##'         this then also nicely works with "= 0" (if tol=0) which stops in
##'         case the matrices are identical (no change at all).
##'       - We conduct checks of convergence after rearranging each column
##'         (not only after rearranging all d columns)
##'       - The columns of X have to be given in increasing order if !is.sorted!
##'       - No checking here due to speed!
rearrange <- function(X, tol=0, tol.type=c("relative", "absolute"),
                      max.num.col.ra=Inf, method=c("worst", "best"),
                      sample=TRUE, is.sorted=FALSE, trace=FALSE)
{
    N <- nrow(X)
    d <- ncol(X)
    tol.type <- match.arg(tol.type)
    method <- match.arg(method)

    ## Define helper functions
    optim.fun <- if(method=="worst") min else max
    tol.fun <- if(tol.type=="absolute") {
        function(x, y) abs(x-y)
    } else {
        function(x, y) abs((x-y)/y)
    }

    ## Tracing
    if(trace) {
        B <- X
        colnames(B) <- rep("", d)
        print(B)
    }

    ## Keep the sorted X
    X.lst.sorted <- if(is.sorted) {
        .Call(C_col_split, X)
    } else {
        .Call(C_col_split, apply(X, 2, sort)) # need to sort first
    }

    ## Sample the columns (if chosen), compute the initial row sum
    ## and the corresponding min/max row sum
    if(sample) {
        X.lst <- lapply(X.lst.sorted, sample) # list of (resampled) columns of X
        X.rs <- .rowSums(do.call(cbind, X.lst), N, d) # row sums of X
    } else {
        X.lst <- X.lst.sorted # list of columns of X
        X.rs <- .rowSums(X, m=N, n=d) # initial row sum
    }
    m.rs.last.col <- optim.fun(X.rs) # initial minimal row sum

    ## Go through the columns and rearrange one at a time
    iter <- 0 # current iteration number
    j <- 0 # current column number
    num.cols.no.change <- 0 # number of consecutively rearranged columns with no change
    m.row.sums <- c() # vector of minimal/maximal row sums after each considered column
    while (TRUE) {

        ## Update the running indices
        iter <- iter+1 # current iteration number (in IN)
        j <- if(j >= d) 1 else j+1 # current column

        ## Update the working 'matrix'
        Y.lst <- X.lst # define 'matrix' Y (former 'matrix' X) to work with
        Y.rs <- X.rs # row sum of Y (= row sum of X)

        ## Oppositely order the jth column to the sum of all others
        yj <- Y.lst[[j]] # pick out jth column
        rs <- Y.rs - yj # sum over all other columns (but the jth)
        ## Note: The elements of X.lst.sorted are sorted in increasing order
        ##       which is required for oppositely reordering them
        yj. <- X.lst.sorted[[j]][indices_opp_ordered_to(rs)] # oppositely reorder Y_j

        ## Update the working 'matrix'
        Y.lst[[j]] <- yj. # update with rearranged jth column
        Y.rs <- rs + yj. # update row sum of Y

        ## Tracing
        if(trace) {
            B <- do.call(cbind, Y.lst)
            colnames(B) <- rep("", d)
            no.change <- identical(yj, yj.)
            colnames(B)[j] <- if(no.change) "=" else "|"
            B <- cbind(B, rs, sum=.rowSums(B, m=N, n=d))
            colnames(B)[d+1] <- paste0("-",j)
            print(B)
        }

        ## Update the vector of computed minimal/maximal row sums
        m.rs.cur.col <- optim.fun(Y.rs) # compute new minimal/maximal row sum
        m.row.sums <- c(m.row.sums, m.rs.cur.col) # append it

        ## Check convergence
        ## Idea: After we ran through all columns once (so iter > d), we compare the
        ##       current tol with the one from d iterations earlier ('last round')
        ##       and check whether the change was small according to our 'convergence'
        ##       criterion.
        ## Note: - This is a bit more efficient than just checking after running through
        ##         all columns (which we would get by checking only if j==d)
        ##       - Checking only two consecutive columns (maybe already before even
        ##         iterating through all columns once) led to bad behavior for ARA()
        ##         in some cases (e.g., real OpRisk data): Both the individual and the joint
        ##         relative error were satisfied but far off with reltol[1]=0.001. Of course
        ##         one could check d consecutive columns for *all* fulfilling the
        ##         'convergence' criterion, but then what's the reached tolerance tol
        ##         if more than two columns are involved? Maybe the maximum computed over the
        ##         previous d many rearranged columns? Probably no gain...
        if(iter > d) {
            m.rs.d.cols.ago <- m.row.sums[iter-d]
            if(is.null(tol)) { # tol = NULL
                num.cols.no.change <- if(identical(yj, yj.)) num.cols.no.change + 1 else 0
                if(num.cols.no.change == d) { # => matrix has not changed in d consecutive col rearrangements
                    tol. <- 0 # as there was no change
                    tol.reached <- TRUE # as we reached 'no change' in d consecutive steps (we don't care wheter max.num.col.ra has been reached)
                    break
                } else { # check whether we have to stop due to max.num.col.ra
                    if(iter == max.num.col.ra) {
                        tol. <- tol.fun(m.rs.cur.col, m.rs.d.cols.ago) # compute the attained tolerance (over last d cols)
                        tol.reached <- FALSE # as num.cols.no.change < d
                        break
                    }
                }
            } else { # tol >= 0
                tol. <- tol.fun(m.rs.cur.col, m.rs.d.cols.ago) # compute the attained tolerance
                tol.reached <- tol. <= tol
                if(iter == max.num.col.ra || tol.reached) break
            }
        }

        ## Updates for the next column rearrangement
        m.rs.last.col <- m.rs.cur.col # update m.rs.last.col
        X.rs <- Y.rs # update the row sums
        X.lst <- Y.lst # update the working 'matrix'

    }

    ## Return
    list(bound=m.rs.cur.col, # computed bound (\underline{s}_N or \overline{s}_N)
         tol=tol., # tolerance for the computed bound
         converged=tol.reached, # indicating whether converged
         m.row.sums=m.row.sums, # the computed row sums after each column rearrangement
         X.rearranged=do.call(cbind, Y.lst)) # the rearranged matrix X
}
