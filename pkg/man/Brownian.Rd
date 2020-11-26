\name{Brownian}
\alias{rBrownian}
\title{Brownian and Related Motions}
\description{
  Simulate paths of dependent, multivariate Brownian motions, geometric
  Brownian motions and Brownian bridges.
}
\usage{
rBrownian(N, t, d = 1, U = matrix(runif(N * n * d), ncol = d),
          drift = rep(0, d), vola = rep(1, d),
          type = c("BM", "GBM", "BB"), init = rep(1, d))
}
\arguments{
  \item{N}{number \eqn{N} of paths to simulate (positive integer).}
  \item{t}{\eqn{n}-vector of the form
    \eqn{(t_1,\dots,t_n)}{(t_1,...,t_n)}
    with \eqn{0 < t_1 < \dots < t_n}{0 < t_1 < ... < t_n}
    containing the time points where to simulate the stochastic
    process. Note that the return object of \code{rBrownian()}
    contains the values at the time points
    \eqn{0, t_1,\dots, t_n}{0, t_1, ..., t_n} then.}
  \item{d}{dimension \eqn{d} of the stochastic process, so the number
    of components (positive integer).}
  \item{U}{\eqn{(N\cdot n, d)}{(N*n, d)}-\code{\link{matrix}} of copula
    realizations to be converted to the joint increments of
    the stochastic process.}
  \item{drift}{\eqn{d}-\code{\link{vector}} of drifts
    (typically denoted by \eqn{\mu}{mu}). Note that risk-neutral
    drifts are \eqn{r - \sigma^2/2}{r - sigma^2/2}, where \eqn{r}
    is the risk-free interest rate and \eqn{\sigma}{sigma}
    the volatility.}
  \item{vola}{\eqn{d}-\code{\link{vector}} of volatilities
    (typically denoted by \eqn{\sigma}{sigma}).}
  \item{type}{\code{\link{character}} string indicating whether
    a Brownian motion (\code{"BM"}), geometric Brownian motion
    (\code{"GBM"}) or Brownian bridge (\code{"BB"}) is to be
    simulated.}
  \item{init}{initial values (typically stock prices at time 0)
    for \code{type = "GBM"}.}
}
\value{
  \code{rBrownian()} returns an \eqn{(N, n+1, d)}-array containing
  the \eqn{N} paths of the specified \eqn{d}-dimensional
  stochastic process simulated at the \eqn{n+1} time points
  (\eqn{0, t_1,\dots, t_n}{0, t_1, ..., t_n}).
}
\author{Marius Hofert}
\examples{
## Setup
d <- 3 # dimension (3 stock prices)
library(copula)
tcop <- tCopula(iTau(tCopula(), tau = 0.5), dim = d, df = 4) # t_4 copula
vola <- seq(0.05, 0.20, length.out = d) # volatilities sigma
r <- 0.01 # risk-free interest rate
drift <- r - vola^2/2 # marginal drifts
init <- seq(10, 100, length.out = d) # initial stock prices
N <- 1000 # number of replications
n <- 250 # number of time steps

## Simulate N paths of a cross-sectionally dependent d-dimensional
## geometric Brownian motion over n time steps
set.seed(271)
U <- rCopula(N * n, copula = tcop)
res <- rBrownian(N, t = 1:n, d = d, U = U, drift = drift, vola = vola,
                 type = "GBM", init = init)
pairs(res[,2,], gap = 0, pch = ".") # plot after one time step

## Now based on quasi-random numbers
library(qrng)
set.seed(271)
U. <- cCopula(to_array(sobol(N, d = d * n, randomize = "digital.shift"), f = n),
              copula = tcop, inverse = TRUE)
res. <- rBrownian(N, t = 1:n, d = d, U = U., drift = drift, vola = vola,
                  type = "GBM", init = init)
pairs(res.[,2,], gap = 0, pch = ".")

## Brownian bridge
N <- 8
t <- 1:100/100
B <- rBrownian(N, t = t, type = "BB")
plot(NA, xlim = 0:1, ylim = range(B),
     xlab = "Time t", ylab = expression("Brownian bridge path"~(B[t])))
for(i in 1:N)
    lines(c(0,t), B[i,,], col = i)
}
\keyword{datagen}