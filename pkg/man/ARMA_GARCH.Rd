\name{ARMA_GARCH}
\alias{fit_ARMA_GARCH}
\title{Fitting ARMA-GARCH Processes}
\description{
  Fail-safe componentwise fitting of univariate ARMA-GARCH processes.
}
\usage{
fit_ARMA_GARCH(x, ugarchspec.list = ugarchspec(distribution.model = "std"),
               verbose = TRUE, ...)
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
## Pick out data
data(list = "SP500_const", package = "qrmdata")
ii <- zoo::index(SP500_const)
x <- SP500_const["2005-01-01" <= ii & ii <= "2010-12-31", c("AAPL", "MSFT", "GOOGL")]
stopifnot(all(!is.na(x)))

## Fit ARMA-GARCH models
res <- fit_ARMA_GARCH(x)
str(res, max.level = 1)

## Check out the results
!sapply(res$error, is.null) # indicate errors
!sapply(res$warning, is.null) # indicate warnings
res$warning[[1]]

## Extract residuals (for components for which fitting worked)
resi <- lapply(res$fit, function(z) residuals(z@fit))
str(resi)
}
\keyword{ts}