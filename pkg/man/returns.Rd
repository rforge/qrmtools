\name{returns}
\alias{log_returns}
\title{Compute Log-Returns or the Inverse Transformation}
\description{
  Compute log-returns or data from given
  log-returns (the inverse transformation).
}
\usage{
log_returns(x, inverse = FALSE, start.value)
}
\arguments{
  \item{x}{A matrix of values to be turned into log-returns (if \code{inverse = FALSE})
    or log-returns (if \code{inverse = TRUE}).}
  \item{inverse}{A \code{\link{logical}} indicating whether the inverse transformation
    (data from given log-returns) is to be computed (if TRUE, this
    requires \code{start.value} to be specified).}
  \item{start.value}{If \code{inverse = TRUE}, the last available value
    of the time series to be constructed from the given log-returns \code{x}.}
}
\value{
  A \code{\link{matrix}} with the same number of columns as \code{x},
  just one row less (if \code{inverse = FALSE}) or one row more (if \code{inverse = TRUE}).
}
\details{
  For negative log-returns, use -\code{log_returns(x)} or
  \code{log_returns(-x, inverse = TRUE, start.value = ...)}.
}
\author{Marius Hofert}
\examples{
library(qrmdata)
data(list = "SP500", package = "qrmdata")
x <- SP500

## log-returns
X <- log_returns(x) # build log-returns
Y <- log_returns(X, inverse = TRUE, start.value = x[1]) # transform back (first date is lost then)
stopifnot(all.equal(as.numeric(x), as.numeric(Y)))

## -log-returns
X <- -log_returns(x) # build -log-returns
Y <- log_returns(-X, inverse = TRUE, start.value = x[1]) # transform back
stopifnot(all.equal(as.numeric(x), as.numeric(Y)))
}
\keyword{utilities}