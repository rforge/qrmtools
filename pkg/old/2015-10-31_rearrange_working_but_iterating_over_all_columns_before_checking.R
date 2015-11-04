##' @title Basic rearrangement function for (A)RA
##' @param X (N, d)-matrix \underline{X}^\alpha or \overline{X}^\alpha
##' @param tol Tolerance to determine (the individual) convergence;
##'        if NULL, the iteration is done until the matrix doesn't change
##' @param tol.type Character string indicating the tolerance function used
##'        ("relative" or "absolute")
##' @param maxiter Maximal number of iterations
##' @param method Character indicating which VaR is approximated (worst/best)
##'        determines optimizing function (min for worst VaR; max
##'        for best VaR)
##' @param sample A logical indicating whether each column of the working
##'        matrix is sampled before the iteration begins
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
##'            row sums after each iteration over all cols
##'         5) The (optimally) rearranged (N, d)-matrix
##' @author Marius Hofert and Kurt Hornik
##' @note - We use "<= tol" to determine convergence instead of "< tol" as
##'         this then also nicely works with "= 0" (if tol=0) which stops in
##'         case the matrices are identical (no change at all).
##'       - No checking here due to speed!
##'       - The columns of X have to be given in increasing order if !is.sorted!
rearrange <- function(X, tol=0, tol.type=c("relative", "absolute"),
                      maxiter=Inf, method=c("worst", "best"),
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

    ## Output initial matrix
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
    m.rs.old <- optim.fun(X.rs) # initial minimal row sum

    ## Loop through the columns
    m.row.sums <- c() # vector of minimal/maximal row sums after each iteration over all cols
    while (TRUE) {

        ## Oppositely order X (=> Y)
        ## Note: - The elements of X.lst.sorted are in increasing order
        ##         => required for oppositely reordering them
        ##       - One could check whether d consecutive column-rearrangements
        ##         did not lead to a change and then stop (as all columns are
        ##         oppositely ordered to the sum of all others in this case).
        ##         This is doable for smaller matrices, but typically neither
        ##         the case nor efficient to do for larger matrices.
        Y.lst <- X.lst
        Y.rs <- X.rs # row sum of Y
        for(j in 1:d) { # one iteration over all columns of the matrix
            yj <- Y.lst[[j]] # jth column of Y
            rs <- Y.rs - yj # sum over all other columns (but the jth)
            yj <- X.lst.sorted[[j]][indices_opp_ordered_to(rs)] # oppositely reorder
            Y.lst[[j]] <- yj # update list with rearranged jth column
            Y.rs <- rs + yj # update row sum of Y
            if(trace) { # for debugging
                B <- do.call(cbind, Y.lst)
                colnames(B) <- rep("", d)
                colnames(B)[j] <- "|"
                B <- cbind(B, rs)
                colnames(B)[d+1] <- paste0("-",j)
                print(B)
            }
        }

        ## Compute the minimal/maximal row sums
        m.rs.new <- optim.fun(Y.rs) # compute new minimal/maximal row sum
        m.row.sums <- c(m.row.sums, m.rs.new) # append the new minimal/maximal row sum

        ## Check convergence (we use "<= tol" as it allows for tol=0)
        tol. <- tol.fun(m.rs.new, m.rs.old) # attained tolerance
        tol.reached <- if(is.null(tol)) { # tol reached?
            identical(Y.lst, X.lst)
        } else {
            tol. <= tol
        }
        if(length(m.row.sums) == maxiter || tol.reached) {
            break
        } else {
            m.rs.old <- m.rs.new # update m.rs.old
            X.rs <- Y.rs
            X.lst <- Y.lst
        }

    }

    ## Return
    list(bound=m.rs.new, # computed bound (\underline{s}_N or \overline{s}_N)
         tol=tol., # tolerance for the computed bound
         converged=tol.reached, # indicating whether converged
         m.row.sums=m.row.sums, # the computed row sums after each iteration through all cols
         X.rearranged=do.call(cbind, Y.lst)) # the rearranged matrix X
}
