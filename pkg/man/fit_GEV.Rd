\name{fit_GEV}
\alias{logLik_GEV}
\alias{fit_GEV}
\title{Maximum Likelihood Fitting of Generalized Extreme Value Distributions}
\description{
  Log-likelihood and maximum-likelihood fitting for the generalized
  extreme value distribution (GEV).
}
\usage{
logLik_GEV(param, x)
fit_GEV(x, init = NULL, estimate.cov = TRUE, control = list(), ...)
}
\arguments{
  \item{param}{vector containing the value of the shape \eqn{\xi}{xi} (a
    real), location \eqn{\mu}{mu} (a real) and scale \eqn{\sigma}{sigma}
    (positive real) parameters of the GEV distribution in this order.}
  \item{x}{numeric vector of data. In the block maxima method, these are the
    block maxima.}
  \item{init}{vector of length three specifying the initial values
    for \eqn{\xi}{xi}, \eqn{\mu}{mu}, \eqn{\sigma}{sigma};
    if \code{NULL}, determined by a procedure to find initial values.}
  \item{estimate.cov}{\code{\link{logical}} indicating whether the
    asymptotic covariance matrix of the parameter estimators is to be
    estimated (inverse of observed Fisher information (negative Hessian
    of log-likelihood evaluated at MLE)).}
  \item{control}{\code{\link{list}}; passed to the underlying
    \code{\link{optim}()}.}
  \item{\dots}{additional arguments passed to the underlying
    \code{\link{optim}()}.}
}
\value{
  \code{logLik_GEV()} computes the log-likelihood of the GEV
  (\code{-Inf} if not admissible).

  \code{fit_GEV()} returns the return object of \code{\link{optim}()}
  and, appended, the estimated asymptotic covariance matrix of the
  parameter estimators, if \code{estimate.cov}.
}
\details{
  Caution: See Coles (2001, p. 55) for how to interpret \eqn{\xi\le
    -0.5}{xi <= -0.5}; in particular, the standard asymptotic properties
  of the MLE do not apply.
}
\author{Marius Hofert}
\references{
   McNeil, A. J., Frey, R., and Embrechts, P. (2015).
   \emph{Quantitative Risk Management: Concepts, Techniques, Tools}.
   Princeton University Press.

   Coles, S. (2001).
   \emph{An Introduction to Statistical Modeling of Extreme Values}.
   Springer-Verlag.
}
\examples{
## Simulate some data
xi <- 0.5
mu <- -2
sig <- 3
n <- 1000
set.seed(271)
X <- rGEV(n, xi = xi, mu = mu, sigma = sig)

## Plot the log-likelihood for fixed mu and sigma (as generated)
xi. <- seq(-0.1, 0.8, length.out = 65)
LL <- sapply(xi., function(xi..) logLik_GEV(c(xi.., mu, sig), x = X))
plot(xi., LL, type = "l", xlab = expression(xi),
     ylab = expression("Log-likelihood for fixed"~mu~"and"~sigma))
## => Numerically non-trivial (for this seed!)

## Plot the profile likelihood for given xi's
## Note: optim() fails if logLik_GEV is -Inf (which can happen for the specified
##       xi's). We could simply choose a smaller xi range to avoid this problem
##       but simply not call optim() in these cases (these cases are not
##       interesting anyways)
pLL <- sapply(xi., function(xi..) {
    nuisance.init <- qrmtools:::fit_GEV_init(X)[2:3] # initialize mu, sigma
    if(!is.finite(logLik_GEV(c(xi.., nuisance.init), x = X))) {
        NA
    } else {
        optim(nuisance.init, fn = function(nuisance)
            logLik_GEV(c(xi.., nuisance), x = X),
            control = list(fnscale = -1))$value
    }
})
plot(xi., pLL, type = "l", xlab = expression(xi), ylab = "Profile log-likelihood")

## Fit
(fit <- fit_GEV(X))
(est <- fit$par) # estimates of xi, mu, sigma
stopifnot(all.equal(est[1], xi, tol = 0.07),
          all.equal(est[2], mu, tol = 0.12),
          all.equal(est[3], sig, tol = 0.06))
sqrt(diag(fit$Cov)) # estimated asymp. variances of MLEs = std. errors of MLEs
}
\keyword{distribution}
