### Working with returns #######################################################

##' @title Compute Log-Returns or the Inverse Transformation
##' @param x A matrix of values to be turned into log-returns (if inverse = FALSE)
##'        or log-returns to be turned into data (if inverse = TRUE)
##' @param inverse A logical indicating whether the inverse transformation
##'        (data from given log-returns) is to be computed (if TRUE, this
##'        requires start.value to be specified)
##' @param start.value If inverse = TRUE, the last available value of the time
##'        series
##' @param drop A logical indicating whether 1-column matrices which are not xts
##'        objects are returned as vectors
##' @return A matrix containing the log-returns or their 'inverses'
##' @author Marius Hofert
##' @note For *negative* log-returns, use -log_returns(x) or
##'       log_returns(-x, inverse = TRUE, start.value = ...)
log_returns <- function(x, inverse = FALSE, start.value, drop = TRUE)
{
    if(!is.matrix(x)) x <- cbind(x)
    if(inverse) {
        ## Note:
        ## X_t = log(S_t/S_{t-1})
        ## => S_t = S_{t-1} * exp(X_t) = ... = S_{last index} * exp(X_1 + X_2 + .. + X_t)
        d <- ncol(x)
        stopifnot(!missing(start.value), length(start.value) == d)
        x.csum <- rbind(rep(0, d),
                        apply(x, 2, cumsum)) # 'xts' lost here
        start.value.factors <- matrix(rep(start.value, each = nrow(x.csum)), ncol = d)
        r <- start.value.factors * exp(x.csum)
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