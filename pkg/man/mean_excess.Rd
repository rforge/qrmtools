\name{mean_excess}
\alias{mean_excess}
\alias{mean_excess_plot}
\title{Mean Excess}
\description{
  Sample mean excess function and sample mean excess plot.
}
\usage{
mean_excess(x, omit = 3)
mean_excess_plot(x, omit = 3,
                 xlab = "Threshold", ylab = "Mean excess over threshold", ...)
}
\arguments{
  \item{x}{\code{\link{numeric}} vector of data.}
  \item{omit}{number \eqn{\ge 1}{>= 1} of unique % differs from QRM::MEplot()!
    last observations to be omitted from the sorted data (as mean excess plot
    becomes unreliable for these observations as thresholds).}
  \item{xlab}{x-axis label.}
  \item{ylab}{y-axis label.}
  \item{\dots}{additional arguments passed to the underlying
    \code{\link{plot}()}.}
}
\value{
  \code{mean_excess()} returns a two-column matrix giving
  the sorted data without the \code{omit}-largest unique values
  (first column) and the corresponding values of the sample mean excess
  function (second column).

  \code{mean_excess_plot()} returns \code{invisible()}.
}
\details{
  Mean excess plots can be used in the peaks-over-threshold method for
  choosing a threshold. To this end, one chooses the smallest threshold
  above which the mean excess plot is roughly linear.
}
\author{Marius Hofert}
\examples{
## (Sample) mean excess function
data(fire)
ME <- mean_excess(fire)
stopifnot(dim(ME) == c(2164, 2),
          all.equal(ME[nrow(ME),], c(65.707491, 121.066231),
                    check.attributes = FALSE))
## (Sample) mean excess plot
mean_excess_plot(fire)
## => Any value in [10, 20] seems reasonable here as threshold choice
##    (one would probably go with 10 to benefit from a larger sample size).
plot(ME, xlab = "Threshold", ylab = "Mean excess over threshold") # equivalent
}
\keyword{nonparametric}
\keyword{hplot}