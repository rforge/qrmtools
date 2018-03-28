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
    if \code{NULL}, determined internally.}
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
  Initial values are computed from \eqn{\xi = 0}{xi = 0} (but actually
  with a small positive value as Nelder--Mead could fail otherwise).
  Furthermore, \eqn{\sigma}{sigma} (obtained from the case
  \eqn{\xi = 0}{xi = 0} is doubled in order to guarantee a finite
  log-likelihood for the first evaluation (at the initial values).

  After several experiments (see the source code), one can safely
  say that finding initial values for fitting GEVs via MLE is
  non-trivial; see also the block maxima method script about the
  Black Monday event on \url{qrmtutorial.org}.

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

## Fit
(fit <- fit_GEV(X))
(est <- fit$par) # estimates of xi, mu, sigma
stopifnot(all.equal(est[1], xi, tol = 0.07),
          all.equal(est[2], mu, tol = 0.12),
          all.equal(est[3], sig, tol = 0.06))
sqrt(diag(fit$Cov)) # estimated asymp. variances of MLEs = std. errors of MLEs

## Plot the log-likelihood in xi for fixed mu and sigma (fixed as generated)
xi. <- seq(-0.1, 0.8, length.out = 65)
logLik <- sapply(xi., function(xi..) logLik_GEV(c(xi.., mu, sig), x = X))
plot(xi., logLik, type = "l", xlab = expression(xi),
     ylab = expression("Log-likelihood for fixed"~mu~"and"~sigma))
## => Numerically quite challenging (for this seed!)

## Plot the profile likelihood for these xi's
## Note: As initial values for the nuisance parameters mu, sigma, we
##       use their values in the case xi = 0 (for all fixed xi = xi. as above
##       in particular those xi != 0). Furthermore, for the given data X
##       and xi = xi., we make sure the initial value for sigma is so large
##       that the density is not 0 and thus the log-likelihood is finite.
pLL <- sapply(xi., function(xi..) {
    sig.init <- sqrt(6 * var(X)) / pi
    mu.init <- mean(X) - sig.init * 0.5772157
    while(!is.finite(logLik_GEV(c(xi.., mu.init, sig.init), x = X)))
        sig.init <- sig.init * 2
    optim(c(mu.init, sig.init), fn = function(nuis)
                logLik_GEV(c(xi.., nuis), x = X),
    		        control = list(fnscale = -1))$value
})
plot(xi., pLL, type = "l", xlab = expression(xi),
     ylab = "Profile log-likelihood")
}
\keyword{distribution}
