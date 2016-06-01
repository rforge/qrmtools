\name{pp_qq_plot}
\alias{pp_plot}
\alias{qq_plot}
\title{P-P and Q-Q Plots}
\description{
  Probability-probability plots and quantile-quantile plots.
}
\usage{
pp_plot(x, FUN, xlab = "Theoretical probabilities",
        ylab = "Sample probabilities", ...)
qq_plot(x, FUN, xlab = "Theoretical quantiles", ylab = "Sample quantiles", ...)
}
\arguments{
  \item{x}{The data \code{\link{vector}}.}
  \item{FUN}{A \code{\link{function}}. For
    \describe{
      \item{\code{pp_plot()}:}{The distribution function.}
      \item{\code{qq_plot()}:}{The quantile function.}
  }}
  \item{xlab}{The x-axis label.}
  \item{ylab}{The y-axis label.}
  \item{\dots}{Additional arguments passed to the underlying function
    \code{\link{plot}()}.}
}
\value{
  \code{\link{invisible}()}.
}
\details{
  Note that Q-Q plots are more widely used (as they make deviations in
  the tails more visible).
}
\author{Marius Hofert}
\examples{
## Generate data
n <- 1000
mu <- 1
sig <- 3
nu <- 3.5
set.seed(271) # set seed
x <- mu + sig * sqrt((nu-2)/nu) * rt(n, df = nu) # sample from t_nu(mu, sig^2)

## P-P plot
pF <- function(q) pt((q - mu) / (sig * sqrt((nu-2)/nu)), df = nu)
pp_plot(x, FUN = pF)

## Q-Q plot
qF <- function(p) mu + sig * sqrt((nu-2)/nu) * qt(p, df = nu)
qq_plot(x, FUN = qF)

## A comparison with R's qqplot() and qqline()
qqplot(qF(ppoints(length(x))), x) # the same (except labels)
qqline(x, distribution = qF) # slightly different (since *estimated*)
}
\keyword{hplot}