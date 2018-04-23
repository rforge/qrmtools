\name{fit_GPD}
\alias{logLik_GPD}
\alias{fit_GPD_MLE}
\title{Maximum Likelihood Fitting of the Generalized Pareto Distribution}
\description{
  Log-likelihood and maximum-likelihood fitting for the generalized
  Pareto distribution (GPD).
}
\usage{
logLik_GPD(param, x)
fit_GPD_MLE(x, init = NULL, estimate.cov = TRUE, control = list(), ...)
}
\arguments{
  \item{param}{vector containing the value of the shape \eqn{\xi}{xi} (a
    real) and scale \eqn{\beta}{beta}
    (positive real) parameters of the GPD in this order.}
  \item{x}{numeric vector of data. In the peaks-over-threshold method,
    these are the excesses (exceedances minus threshold).}
  \item{init}{vector of length two specifying the initial values
    for \eqn{\xi}{xi} and \eqn{\beta}{beta};
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
  \code{logLik_GPD()} computes the log-likelihood of the GPD
  (\code{-Inf} if not admissible).

  \code{fit_GPD_MLE()} returns the return object of \code{\link{optim}()}
  and, appended, the estimated asymptotic covariance matrix of the
  parameter estimators, if \code{estimate.cov}.
}
\details{
  Initial values are computed via probability weighted moments.
  Note that, in contrast to the GEV, there is no guarantee that
  the log-likelihood at the initial value is finite and so the
  underlying
  \code{\link{optim}()} may fail.
}
\author{Marius Hofert}
\references{
   McNeil, A. J., Frey, R. and Embrechts, P. (2015).
   \emph{Quantitative Risk Management: Concepts, Techniques, Tools}.
   Princeton University Press.

   Hosking, J. R. M. and Wallis, J. R. (1987).
   Parameter and Quantile Estimation for the Generalized Pareto Distribution.
   \emph{Technometrics} \bold{29}(3), 339--349.

   Landwehr, J. M., Matalas, N. C. and Wallis, J. R. (1979).
   Estimation of Parameters and Quantiles of Wakeby Distributions.
   \emph{Water Resourches Research} \bold{15}(6), 1361--1379.
}
\examples{
## Simulate some data
xi <- 0.5
beta <- 3
n <- 1000
set.seed(271)
X <- rGPD(n, xi = xi, beta = beta)

## Fit
(fit <- fit_GPD_MLE(X))
(est <- fit$par) # estimates of xi, mu, sigma
stopifnot(all.equal(est[["xi"]],   xi,   tol = 0.12),
          all.equal(est[["beta"]], beta, tol = 0.11))
fit$SE # estimated asymp. variances of MLEs = std. errors of MLEs

## Plot the log-likelihood in xi for fixed beta (fixed as generated)
xi. <- seq(-0.1, 0.8, length.out = 65)
logLik <- sapply(xi., function(xi..) logLik_GPD(c(xi.., beta), x = X))
plot(xi., logLik, type = "l", xlab = expression(xi),
     ylab = expression("GPD log-likelihood for fixed"~beta))

## Plot the profile likelihood for these xi's
## (with an initial interval for the nuisance parameter beta such that
##  logLik_GPD() is finite)
pLL <- sapply(xi., function(xi..) {
    ## Choose beta interval for optimize()
    int <- if(xi.. >= 0) {
               ## Method-of-Moment estimator
               mu.hat <- mean(X)
               sig2.hat <- var(X)
               xi.hat <- (1-mu.hat^2/sig2.hat)/2
               beta.hat <- mu.hat*(1-xi.hat)
               ## log-likelihood always fine for xi.. >= 0 for all beta
               c(1e-8, 2 * beta.hat)
           } else { # xi.. < 0
               ## Make sure logLik_GPD() is finite at endpoints of int
               mx <- max(X)
               -xi.. * mx * c(1.01, 100) # -xi * max(X) * scaling
               ## Note: for xi.. closer to 0, the upper scaling factor
               ##       needs to be chosen sufficiently large in order
               ##       for optimize() to find an optimum (not just the
               ##       upper end point). Try it with '2' instead of '100'.
           }
    ## Optimization
    optimize(function(nuis) logLik_GPD(c(xi.., nuis), x = X),
             interval = int, maximum = TRUE)$maximum
})
plot(xi., pLL, type = "l", xlab = expression(xi),
     ylab = "GPD profile log-likelihood")
}
\keyword{distribution}
