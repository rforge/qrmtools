\name{alloc}
\alias{alloc_ellip}
\title{(Capital) allocation}
\description{
  Computing allocations.
}
\usage{
alloc_ellip(total, loc, scale)
}
\arguments{
  \item{total}{total to be allocated (typically the risk measure of the
    sum of the underlying loss random variables).}
  \item{loc}{location vector of the elliptical distribution of the loss
    random vector.}
  \item{scale}{scale (covariance) matrix of the elliptical distribution of
    the loss random vector.}
}
\value{
  \eqn{d}-vector of allocated amounts (the allocation) according to the
  Euler principle under the assumption that the underlying loss random
  vector follows a \eqn{d}-dimensional elliptical distribution with
  location vector \code{loc} (\eqn{\bm{mu}}{mu} in the reference) and
  scale matrix \code{scale} (\eqn{\Sigma}{Sigma} in the reference, a
  covariance matrix) and that the risk measure is law-invariant,
  positive-homogeneous and translation invariant.
}
\details{
  The result for \code{loc = 0} can be found in McNeil et al. (2015,
  Corollary 8.43). Otherwise, McNeil et al. (2015,
  Theorem 8.28 (1)) can be used to derive the result.
}
\author{Marius Hofert}
\references{
   McNeil, A. J., Frey, R. and Embrechts, P. (2015).
   \emph{Quantitative Risk Management: Concepts, Techniques, Tools}.
   Princeton University Press.
}
\examples{
alloc_ellip # ... short enough to learn from the source
}
\keyword{models}