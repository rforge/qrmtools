### Allocation #################################################################

##' @title Euler Allocations for Elliptical Distributions
##' @param total total to be allocated (typically the risk measure of the
##'        sum of the underlying loss random variables)
##' @param loc location vector of the elliptical distribution of the loss
##'        random vector
##' @param Sigma scale (covariance) matrix of the elliptical distribution of
##'        the loss random vector
##' @return allocation according to the Euler principle
##' @author Marius Hofert and Takaaki Koike
##' @note 1) See MFE (2015, Corollary 8.43) for bm(L) ~ E_d(bm{0}, Sigma, psi)
##'          and positive-homogeneous and law-invariant risk measures.
##'          After summing over all i, you get AC (total) / AC_j = sum(Sigma) /
##'          rowSums(Sigma), or, equivalently, rowSums(Sigma) / sum(Sigma)
##'          = AC_j / AC. Since AC = total, we obtain AC_j = (rowSums(Sigma) /
##'          sum(Sigma)) * total.
##'       2) If, additionally, the risk measure is also translation invariant,
##'          MFE (2015, Theorem 8.28 (1)) implies that r_varrho(bm(lambda))
##'          gets an additional summand bm(lambda)^T bm(mu). When differentiated
##'          w.r.t. lambda_i, this implies a summand of mu_i, so tilde(AC)_j
##'          := AC_j - mu_j satisfies 1) and thus tilde(AC)_j = (rowSums(Sigma) /
##'          sum(Sigma)) * total and thus AC_j = mu_j + (rowSums(Sigma) /
##'          sum(Sigma)) * total
alloc_ellip <- function(total, loc, scale)
{
    d <- ncol(scale)
    if(length(loc) == 1) loc <- rep(loc, d)
    stopifnot(nrow(scale) == d, length(loc) == d)
    rs <- rowSums(scale)
    loc + (rs / sum(rs)) * total # Euler allocation
}

##' @title Sub-sampling based on a Risk Measure of the Sum
##' @param x (n, d)-data matrix
##' @param level confidence level(s)
##' @param risk.measure character string or function specifying the risk measure
##'        computed from the row sums of x based on the given level(s) in order
##'        to determine the conditioning region.
##' @param ... additional arguments passed to 'risk.measure'
##' @return sub-sample of x that satisfies that its row sums are in the
##'         conditioning region specified by 'risk.measure' computed from the row
##'         sums of x and the given 'level'
##' @author Marius Hofert
##' @note This function could at some point get an argument 'type' for using,
##'       say, the row maximum instead of the row sum to determine the
##'       conditioning region.
conditioning <- function(x, level, risk.measure = "VaR_np", ...)
{
    ## Basics
    if(!is.matrix(x))
        x <- as.matrix(x)
    if(length(level) == 1) level <- c(level, 1)
    stopifnot(length(level) == 2, 0 <= level, level <= 1)
    is.function.rm <- is.function(risk.measure)
    if(!is.function.rm) {
        stopifnot(is.character(risk.measure), existsFunction(risk.measure)) # check
        ## Create an expression of the function call and evaluate that below
        expr <- as.call(c(as.name(risk.measure), # 'unquote' string
                          quote(x), # includes 'x' as first argument (which exists in here)
                          quote(level.), # includes placeholder to be substituted below
                          as.expression(list(...)))) # includes '...' (which exists in here)
    }

    ## Estimate risk.measure(S) for the two confidence levels
    ## Note: We could also call risk.measure(S, level = level, ...) directly
    ##       but that assumes the provided risk.measure() is vectorized in 'level'
    S <- rowSums(x) # row sums
    rm.level.S <- if(is.function.rm) {
                      sapply(level, function(l) risk.measure(S, level = l, ...))
                  } else {
                      sapply(level, function(l) eval(expr, list(level. = l)))
                  }
    if(length(rm.level.S) != 2) # sanity check
        stop("The output of the evaluated 'risk.measure' does not have length 2.")

    ## Return sub-sample
    x[(rm.level.S[1] < S) & (S <= rm.level.S[2]),]
}
