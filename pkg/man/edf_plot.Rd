\name{edf_plot}
\alias{edf_plot}
\title{Plot an Empirical Distribution Function}
\description{
  Plot an empirical distribution function possibly overlaid with
  a analytical distribution function.
}
\usage{
edf_plot(x, FUN = NULL, FUN.args = NULL,
         xlim = c(min(x), extendrange(x)[2]),
         q = c(max(x[1], xlim[1]), xlim[2]), length.out = 129,
         log = FALSE, xlab = "Value", ylab = "Probability", ...)
}
\arguments{
  \item{x}{data \code{\link{vector}} of which the empirical distribution
    function is to be plotted.}
  \item{FUN}{\code{\link{function}()} to be plotted additionally (unless
    \code{NULL}).}
  \item{FUN.args}{\code{\link{list}} of arguments passed to the
    underlying \code{\link{lines}()} for drawing \code{FUN}.}
  \item{xlim}{x-axis limit for determining where to evaluate the \code{ecdf(x)}.}
  \item{q}{\code{\link{numeric}} vector of quantiles where \code{FUN} is
    evaluated; if of length 2, \code{FUN} is evaluated at
    \code{length.out}-many points in \code{q[1]} to \code{q[2]} (exponential
    scale if \code{log = TRUE}).}
  \item{length.out}{number of points for constructing a sequence between
    \code{q[1]} and \code{q[2]} unless \code{q} is longer than 2.}
  \item{log}{\code{\link{logical}} indicating whether the x-axis is
    given in logarithmic scale.}
  \item{xlab}{x-axis label.}
  \item{ylab}{y-axis label.}
  \item{\dots}{additional arguments passed to the underlying
    \code{\link{plot}()}.}
}
\value{
  \code{\link{invisible}()}.
}
\details{
  Note that a Q-Q plot is typically more meaningful to check a sample
  against a theoretical distribution.
}
\author{Marius Hofert}
\examples{
data(fire)
u <- 10 # threshold
exceed <- fire[fire > u] # exceedances
excess <- exceed - u # excesses
fit <- fit_GPD_MLE(excess) # fit GPD to excesses
df <- function(q) # define fitted GPD
    pGPD(q, shape = fit$par[["shape"]], scale = fit$par[["scale"]])

## Plot empirical distribution function of excesses and exceedances
## (partially with fitted GPD overlaid)
edf_plot(excess)
edf_plot(excess, log = TRUE)
edf_plot(excess, FUN = df)
edf_plot(excess, FUN = df, log = TRUE)
edf_plot(excess, FUN = df, log = TRUE, FUN.args = list(col = "royalblue3"))
edf_plot(exceed, FUN = function(q) df(q-u), log = TRUE)

## Q-Q plots
qf <- function(p) # quantile function of df
    qGPD(p, shape = fit$par[["shape"]], scale = fit$par[["scale"]])
qq_plot(excess, FUN = qf)
qq_plot(exceed, FUN = function(p) u + qf(p))
}
\keyword{hplot}