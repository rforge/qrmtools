\name{worst_VaR}
\alias{worst_VaR_hom}
\alias{dual_bound}
\alias{crude_VaR_bounds}
\title{Worst Value-at-Risk for Given Margins}
\description{
  Compute the worst Value-at-Risk for given marginal distributions.
}
\usage{
worst_VaR_hom(alpha, d, interval=NULL, method=c("Wang", "dual", "Pareto"), ...)
dual_bound(s, d, pF, ...)
crude_VaR_bounds(alpha, d, qF, ...)
}
\arguments{
  \item{alpha}{Value-at-Risk confidence level (e.g., 0.99).}
  \item{d}{Dimension (number of risk factors).}
  \item{interval}{Initial interval (a \code{\link{numeric}(2)}). For
    \code{method="Wang"} and \code{method="Pareto"}, if not provided,
    the initivial interval is
    \eqn{[0,(1-\alpha)/d]}{[0,(1-alpha)/d]}. For \code{method="dual"},
    note that the lower endpoint of the initial interval has to be
    sufficiently large in order for the the inner root-finding algorithm
    to find a root; see Details.}
  \item{method}{Method (a \code{\link{character}}).}
  \item{pF, qF}{The loss distribution and quantile functions,
    respectively; the latter can be a \code{\link{list}} of quantile
    functions.}
  \item{s}{Dual bound evaluation point.}
  \item{\dots}{For
    \describe{
      \item{\code{worst_VaR_hom()}:}{For \code{method="Wang"} or
	\code{method="Pareto"} additional arguments passed to
	\code{\link{uniroot}()}; for \code{method="dual"} additional
	arguments passed to the underlying
	\code{\link{integrate}()} for evaluating the dual bound
	\eqn{D(s)}.}
      \item{\code{dual_bound()}:}{Additional arugments passed to the
	underlying \code{\link{integrate}()}.}
      \item{\code{crude_VaR_bounds()}:}{Additional arguments passed to
	(all provided) quantile functions.}
    }}
}
\value{
  \code{worst_VaR_hom()} returns the worst Value-at-Risk at confidence
  level \eqn{\alpha}{alpha} for \eqn{d} risks with equal distribution
  function specified by the ellipsis \code{...}.

  \code{dual_bound()} returns the value of the dual bound \eqn{D(s)} as
  given in Embrechts, Puccetti, \enc{Rüschendorf}{Rueschendorf} (2013, Eq. (12)).

  \code{crude_VaR_bounds()} returns crude lower and upper bounds for
  Value-at-Risk at confidence level \eqn{\alpha}{alpha} for any
  \eqn{d}-dimensional model with marginal quantile functions
  specified by \code{qF}.
}
\details{
  For \code{method="Wang"} and \code{method="Pareto"}, the method
  presented in Embrechts, P., Puccetti, G., \enc{Rüschendorf}{Rueschendorf}, L.,
  Wang, R., Beleraj, A. (2014, Proposition 3.1) is implemented. This
  requires one root-finding procedure.

  For \code{method="dual"}, the method presented in
  Embrechts, P., Puccetti, G., \enc{Rüschendorf}{Rueschendorf},
  L. (2013, Proposition 4) is implemented. This
  requires two (nested) root-finding procedures, and (in the general form
  implemented here) a numerical integration. For the inner root-finding
  procedure to find a root, the lower endpoint of the provided initial
  \code{interval} has to be \dQuote{sufficiently large}.
}
\author{Marius Hofert}
\references{
  Embrechts, P., Puccetti, G., \enc{Rüschendorf}{Rueschendorf}, L.,
  Wang, R., Beleraj, A. (2014).
  An Academic Response to Basel 3.5. \emph{Risks} \bold{2}(1), 25--48.

  Embrechts, P., Puccetti, G., \enc{Rüschendorf}{Rueschendorf}, L. (2013).
  Model uncertainty and VaR aggregation. \emph{Journal of Banking \&
    Finance} \bold{37}, 2750--2764.
}
\seealso{
  \code{demo(worst_VaR)}
}
\examples{
## setup
alpha <- 0.99 # VaR confidence level
d <- 8 # dimension
th <- 2 # Pareto parameter
pF <- function(q) pPar(q, theta=th) # marginal distribution function
qF <- function(p) qPar(p, theta=th) # marginal quantile function

## compare the various methods
worst.VaR.Wang <- worst_VaR_hom(alpha, d=d, method="Wang", qF=qF)
init <- crude_VaR_bounds(alpha, d=d, qF=qPar, theta=th)
worst.VaR.dual <- worst_VaR_hom(alpha, d=d, method="dual", interval=init, pF=pF)
##worst.VaR.Par <- worst_VaR_hom(alpha, d=d, method="Pareto", theta=th) # TODO Amir
##stopifnot(all.equal(worst.VaR.Wang, worst.VaR.dual, worst.VaR.Par)) # TODO: not working yet
}
\keyword{distribution}