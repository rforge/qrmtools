##' @title Compute Log-Returns or the Inverse Transformation
##' @param x matrix of values to be turned into log-returns (if inverse = FALSE)
##'        or log-returns to be turned into data (if inverse = TRUE)
##' @param inverse logical indicating whether the inverse transformation
##'        (data from given log-returns) is to be computed (if TRUE, this
##'        requires start to be specified)
##' @param start if inverse = TRUE, the last available value of the time
##'        series
##' @param include.start if inverse = TRUE, a logical indicating whether the
##'        last available value is included
##' @param drop logical indicating whether 1-column matrices which are not xts
##'        objects are returned as vectors
##' @return matrix containing the log-returns or their 'inverses'
##' @author Marius Hofert
##' @note For *negative* log-returns, use -log_returns(x) or
##'       log_returns(-x, inverse = TRUE, start = ...)
log_returns <- function(x, inverse = FALSE, start, include.start = FALSE,
                        drop = TRUE)
{
    if(!is.matrix(x)) x <- cbind(x)
    if(inverse) {
        ## Note:
        ## X_t = log(S_t/S_{t-1})
        ## => S_t = S_{t-1} * exp(X_t) = ... = S_{last index} * exp(X_1 + X_2 + .. + X_t)
        d <- ncol(x)
        stopifnot(!missing(start), length(start) == d)
        x.csum <- apply(x, 2, cumsum) # 'xts' lost here
        if(include.start) x.csum <- rbind(rep(0, d), x.csum) # include S_t
        start.factors <- matrix(rep(start, each = nrow(x.csum)), ncol = d)
        r <- start.factors * exp(x.csum)
        if(drop && ncol(r) == 1) as.vector(r) else r
        ## Note: We didn't incorporate inherits(..., "zoo") here as a zoo object
        ##       would get lost due to the apply() anyways
    } else {
        if(inherits(x, "zoo")) {
            diff(log(x))[-1,] # diff() works componentwise here; leaves the first row in there as NA (remove it here to be consistent)
        } else {
            r <- apply(x, 2, function(x.) diff(log(x.)))
            if(drop && ncol(r) == 1) as.vector(r) else r
        }
    }
}

