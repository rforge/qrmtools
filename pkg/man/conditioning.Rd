\name{conditioning}
\alias{conditioning}
\title{}
\description{

}
\usage{

}
\arguments{
  \item{}{.}
}
\value{

}
\details{

}
\author{Marius Hofert}
\examples{
## Generate data
library(copula)
set.seed(271)
X <- qt(rCopula(1e5, copula = gumbelCopula(2, dim = 5)), df = 3.5)
Y <- conditioning(X, level = 0.9) # => internally level = c(0.9, 1)
stopifnot(nrow(Y) == 1e4) # sanity check (here exactly 1e4)

## An alternative way of calling conditioning()
Y. <- conditioning(X, level = 0.9, risk.measure = "VaR_np")
stopifnot(all.equal(Y., Y)) # sanity check

## Estimate an allocation via MC based on the sub-sample
colMeans(Y) # estimated allocation
apply(Y, 2, sd) / sqrt(nrow(Y)) # estimated standard errors
}
\keyword{utilities}