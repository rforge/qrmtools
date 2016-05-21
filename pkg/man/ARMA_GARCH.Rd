\name{ARMA_GARCH}
\alias{fit_ARMA_GARCH}
\title{Fitting ARMA-GARCH Processes}
\description{
  Fail-safe componentwise fitting of univariate ARMA-GARCH processes.
}
\usage{
fit_ARMA_GARCH(x, ugarchspec.list = ugarchspec(), verbose = TRUE, ...)
}
\arguments{
  \item{x}{A \code{\link{matrix}}-like data structure, possibly an
    \code{xts} object.}
  \item{ugarchspec.list}{An object of class \code{uGARCHspec} (as returned by
    \code{ugarchspec()}) or a list of such. In case of a list, its
    length has to be equal to the number of columns of
    \code{x}. \code{ugarchspec.list} provides the ARMA-GARCH
    specifications for each of the time series (columns of \code{x}).}
  \item{verbose}{A \code{\link{logical}} indicating whether verbose
    output is given.}
  \item{\dots}{Additional arguments passed to the underlying
    \code{\link{ugarchfit}()}.}
}
\value{
  If \code{x} consists of one column only (e.g. a vector),
  \code{ARMA_GARCH()} returns the fitted object; otherwise it returns
  a list of such.
}
\author{Marius Hofert}
\examples{
library(rugarch)
library(copula)

## Read the data, build -log-returns
data(list = "SP500_const", package = "qrmdataConstSP")
stocks <- c("AAPL", "MSFT", "GOOGL")
ii <- zoo::index(SP500_const)
x <- SP500_const["2005-01-01" <= ii & ii <= "2010-12-31", stocks]
stopifnot(all(!is.na(x)))
X <- -log_returns(x)
n <- nrow(X)
d <- ncol(X)

## Fit ARMA-GARCH models to the -log-returns
## See http://unstarched.net/wp-content/uploads/2013/06/an-example-in-rugarch.pdf
uspec <- list(ugarchspec(), # ARMA(1,1)-GARCH(1,1) model for 1st component
              ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(2,2))),
              ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(2,1)),
                         mean.model = list(armaOrder = c(1,2), include.mean = TRUE)))
fitAG <- fit_ARMA_GARCH(X, ugarchspec.list = uspec)
str(fitAG, max.level = 1) # list with components fit, warning, error
stopifnot(sapply(fitAG$error, is.null)) # NULL = no error
stopifnot(sapply(fitAG$warning, is.null)) # NULL = no warning

## Pick out the residuals and plot their pseudo-observations
Z <- sapply(fitAG$fit, residuals)
U <- pobs(Z)
pairs(U, gap = 0, pch = ".")

## Fit a t copula to the residuals
fitC <- fitCopula(tCopula(dim = d, dispstr="un"), data = U, method = "mpl")

## Simulate from the fitted model overall
set.seed(271)
U. <- rCopula(n, fitC@copula) # simulate from the fitted copula
Z. <- qnorm(U.) # transform back to the residuals
## Note: coef(fitAG$fit[[.]]) contains the fitted parameters in case required
##       (for a non-standard residual distribution)
X. <- lapply(1:d, function(j)
    ugarchsim(fitAG$fit[[j]], n.sim = n, m.sim = 1, # 1 simulation of size n
              custom.dist=list(name = "sample", distfit = Z.[,j, drop=FALSE]))@simulation)
str(X., max.level = 3)
X.. <- sapply(X., `[[`, "seriesSim") # simulated time series; use 'residSim' for residuals

## Plots original vs simulated -log-returns
layout(matrix(1:(2*d), ncol = d)) # layout
ran <- range(X, X..)
for(j in 1:d) {
    plot(X[,j],   type = "l", ylim = ran, ylab = paste(stocks[j], "-log-returns"))
    plot(X..[,j], type = "l", ylim = ran, ylab = "Simulated")
}
}
\keyword{ts}