\name{worst_VaR}
\alias{worst_VaR_hom}
\alias{dual_bound}
\alias{crude_VaR_bounds}
\title{Worst Value-at-Risk for Given Margins}
\description{
  Compute the worst Value-at-Risk for given marginal distributions.
}
\usage{
worst_VaR_hom(alpha, d, method=c("Wang", "Wang.Par", "dual"), interval=NULL, ...)
dual_bound(s, d, pF, ...)
crude_VaR_bounds(alpha, d, qF, ...)
}
\arguments{
  \item{alpha}{Value-at-Risk confidence level (e.g., 0.99).}
  \item{d}{Dimension (number of risk factors).}
  \item{method}{Method (a \code{\link{character}}). \code{method="Wang"}
    and \code{method="Wang.Par"} apply the approach of Embrechts,
    Puccetti, \enc{Rüschendorf}{Rueschendorf}, Wang and Beleraj (2014,
    Proposition 3.1); the latter assumes Pareto margins and thus does
    not require numerical integration. \code{method="dual"} applies the
    dual bound approach as in Embrechts, Puccetti,
    \enc{Rüschendorf}{Rueschendorf} (2013, Proposition 4).}
  \item{interval}{Initial interval (a \code{\link{numeric}(2)}). For
    \code{method="Wang"} and \code{method="Wang.Par"}, if not provided,
    the initivial interval is
    \eqn{[0,(1-\alpha)/d]}{[0,(1-alpha)/d]}. For \code{method="dual"},
    note that the lower endpoint of the initial interval has to be
    sufficiently large in order for the the inner root-finding algorithm
    to find a root; see Details.}
  \item{pF, qF}{The loss distribution and quantile functions,
    respectively; the latter can be a \code{\link{list}} of quantile
    functions.}
  \item{s}{Dual bound evaluation point.}
  \item{\dots}{For
    \describe{
      \item{\code{worst_VaR_hom()}:}{For \code{method="Wang"} the
	quantile function \code{qF()} and additional
	arguments passed to the underlying \code{\link{integrate}()}
	have to be provided; for \code{method="dual"} the distribution
	function \code{pF()} and additionalarguments passed to the underlying
	\code{\link{integrate}()} for computing the dual bound
	\eqn{D(s)} have to be provided; and for \code{method="Wang.Par"}
	the Pareto parameter \eqn{\theta>0}{theta>0} has to be provided.}
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
  For \code{method="Wang"} and \code{method="Wang.Par"}, the method
  presented in Embrechts, P., Puccetti, G., \enc{Rüschendorf}{Rueschendorf}, L.,
  Wang, R., Beleraj, A. (2014, Proposition 3.1) is implemented. This
  requires one \code{\link{uniroot}()} and, for the generic
  \code{method="Wang"}, \code{\link{integrate}()}. The critical part
  for these methods is the lower endpoint of the initial interval.
  If the (marginal) distribution function has finite first moment,
  this can be taken as 0. However, if it has infinite first moment,
  the lower endpoint has to be positive (and below the unknown root
  and smaller than the upper endpoint \eqn{(1-\alpha)/d}{(1-alpha)/d}).

  For \code{method="dual"}, the method presented in
  Embrechts, P., Puccetti, G., \enc{Rüschendorf}{Rueschendorf},
  L. (2013, Proposition 4) is implemented. This
  requires two (nested) \code{\link{uniroot}()}, and an
  \code{\link{integrate}()}. For the inner root-finding
  procedure to find a root, the lower endpoint of the provided initial
  \code{interval} has to be \dQuote{sufficiently large}.

  Note that both of these approaches for computing the worst
  Value-at-Risk in the homogeneous case are numerically non-trivial.
  See the source code and \code{demo(worst_VaR)} for more details.
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
I <- crude_VaR_bounds(alpha, d=d, qF=qPar, theta=th)
(worst.VaR.Wang <- worst_VaR_hom(alpha, d=d, method="Wang", qF=qF))
(worst.VaR.Par  <- worst_VaR_hom(alpha, d=d, method="Wang.Par", theta=th))
(worst.VaR.dual <- worst_VaR_hom(alpha, d=d, method="dual", interval=I, pF=pF))
stopifnot(all.equal(worst.VaR.Wang, worst.VaR.Par))
}
\keyword{distribution}