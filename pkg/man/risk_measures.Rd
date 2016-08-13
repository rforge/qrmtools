\name{risk_measures}
\alias{VaR_np}
\alias{VaR_t}
\alias{VaR_Par}
\alias{ES_np}
\alias{ES_t}
\alias{ES_Par}
\title{Risk Measures}
\description{
  Computing risk measures for various distributions.
}
\usage{
VaR_np(x, alpha, names = FALSE, type = 1, ...)
VaR_t(alpha, mu = 0, sigma = 1, df = Inf)
VaR_Par(alpha, theta, kappa = 1)

ES_np(x, alpha, method = c(">", ">="), verbose = FALSE, ...)
ES_t(alpha, mu = 0, sigma = 1, df = Inf)
ES_Par(alpha, theta, kappa = 1)
}
\arguments{
  \item{x}{The vector of losses.}
  \item{alpha}{The confidence level \eqn{\alpha\in[0,1]}{alpha in
      [0,1]}.}
  \item{names}{See \code{?\link{quantile}}.}
  \item{type}{See \code{?\link{quantile}}.}
  \item{mu}{The location parameter.}
  \item{sigma}{The scale parameter, a positive number.}
  \item{df}{The degrees of freedom, a positive number; choose \code{df = Inf}
    for the normal distribution.}
  \item{theta}{The Pareto shape parameter, a positive number.}
  \item{kappa}{The Pareto scale parameter, a positive number.}
  \item{method}{A \code{\link{character}} string indicating the method
    for computing expected shortfall.}
  \item{verbose}{A \code{\link{logical}} indicating whether verbose
    output is given (in case the mean is computed over (too) few
    observations).}
  \item{\dots}{For
    \describe{
      \item{\code{VaR_np()}:}{Additional arguments passed to the
	underlying \code{\link{quantile}()}.}
      \item{\code{ES_np()}:}{Additional arguments passed to \code{VaR_np()}.}
    }
  }
}
\value{
  \code{VaR_np()}, \code{ES_np()} estimate Value-at-Risk and expected
  shortfall non-parametrically. For the latter, the mean over all losses
  (strictly) beyond Value-at-Risk is computed. If \code{method = ">="},
  there is always at least one such loss, whereas if \code{method =
    ">"}, there might be no such loss, in which case \code{NaN} is returned.

  \code{VaR_t()}, \code{ES_t()} compute Value-at-Risk and expected
  shortfall for the \eqn{t} (or normal) distribution.

  \code{VaR_Par()}, \code{ES_Par()} compute Value-at-Risk and expected
  shortfall for the Pareto distribution.
}
\details{
  The distribution function of the Pareto distribution is given by
  \deqn{F(x) = 1-(\kappa/(\kappa+x))^{\theta},\ x\ge 0,}{F(x) = 1-(kappa/(kappa+x))^{theta}, x
  >= 0,} where \eqn{\theta > 0}{theta > 0}, \eqn{\kappa > 0}{kappa > 0}.
}
\author{Marius Hofert}
\references{
   McNeil, A. J., Frey, R., and Embrechts, P. (2015).
   \emph{Quantitative Risk Management: Concepts, Techniques, Tools}.
   Princeton University Press.
}
\examples{
## Generate some losses and (non-parametrically) estimate VaR_alpha and ES_alpha
set.seed(271)
L <- rlnorm(1000, meanlog = -1, sdlog = 2) # L ~ LN(mu, sig^2)
## Note: - meanlog = mean(log(L)) = mu, sdlog = sd(log(L)) = sig
##       - E(L) = exp(mu + (sig^2)/2), var(L) = (exp(sig^2)-1)*exp(2*mu + sig^2)
##         To obtain a sample with E(L) = a and var(L) = b, use:
##         mu = log(a)-log(1+b/a^2)/2 and sig = sqrt(log(1+b/a^2))
VaR_np(L, alpha = 0.99)
ES_np(L,  alpha = 0.99)

## Example 2.16 in McNeil, Frey, Embrechts (2015)
V <- 10000 # value of the portfolio today
sig <- 0.2/sqrt(250) # daily volatility (annualized volatility of 20%)
nu <- 4 # degrees of freedom for the t distribution
alpha <- seq(0.001, 0.999, length.out = 256) # confidence levels
VaRnorm <- VaR_t(alpha, sigma = V*sig, df = Inf)
VaRt4 <- VaR_t(alpha, sigma = V*sig*sqrt((nu-2)/nu), df = nu)
ESnorm <- ES_t(alpha, sigma = V*sig, df = Inf)
ESt4 <- ES_t(alpha, sigma = V*sig*sqrt((nu-2)/nu), df = nu)
ran <- range(VaRnorm, VaRt4, ESnorm, ESt4)
plot(alpha, VaRnorm, type = "l", ylim = ran, xlab = expression(alpha), ylab = "")
lines(alpha, VaRt4, col = "royalblue3")
lines(alpha, ESnorm, col = "darkorange2")
lines(alpha, ESt4, col = "maroon3")
legend("bottomright", bty = "n", lty = rep(1,4), col = c("black",
       "royalblue3", "darkorange3", "maroon3"),
       legend = c(expression(VaR[alpha]~~"for normal model"),
                  expression(VaR[alpha]~~"for "*t[4]*" model"),
                  expression(ES[alpha]~~"for normal model"),
                  expression(ES[alpha]~~"for "*t[4]*" model")))
}
\keyword{models}