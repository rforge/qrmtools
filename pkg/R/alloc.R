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
