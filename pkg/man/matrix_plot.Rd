\name{matrix_plot}
\alias{matrix_plot}
\title{Graphical Tool for Visualizing Matrices}
\description{
  Plot of a matrix.
}
\usage{
matrix_plot(x, xlab = "Column", ylab = "Row",
            scales = list(alternating = c(1,1), tck = c(1,0), x = list(rot = 90)),
            at = NULL, colorkey = NULL, col = c("royalblue3", "white", "maroon3"),
            col.regions = NULL, ...)
}
\arguments{
  \item{x}{A \code{\link{matrix}}-like object.}
  \item{xlab}{The x-axis label.}
  \item{ylab}{The y-axis label.}
  \item{scales}{See \code{\link{levelplot}()}; if \code{\link{NULL}},
    labels and ticks are omitted.}
  \item{at}{See \code{\link{levelplot}()}. If \code{\link{NULL}}, a
    useful default is computed based on the given values in \code{x}.}
  \item{colorkey}{See \code{\link{levelplot}()}. If \code{\link{NULL}}, a
    useful default is computed based on \code{at}.}
  \item{col}{A \code{\link{vector}} of length 3 providing the color
    key's default colors.}
  \item{col.regions}{See \code{\link{levelplot}()}. If \code{\link{NULL}}, a
    useful default is computed based on \code{at}.}
  \item{\dots}{Additional arguments passed to the underlying function
    \code{\link{levelplot}()}.}
}
\value{
  The plot, a Trellis object.
}
\details{
  Plot of a matrix.
}
\author{Marius Hofert}
\examples{
## Generate a random correlation matrix
d <- 50
L <- diag(1:d)
set.seed(271)
L[lower.tri(L)] <- runif(choose(d,2))
Sigma <- L %*% t(L) # ... real, symmetric, positive definite (covariance) matrix Sigma
P <- cor(Sigma)
## Default
matrix_plot(P)
## Default if nonnegative
matrix_plot(abs(P))
## Without diagonal
P. <- abs(P)
diag(P.) <- NA
matrix_plot(P.)
## Default if nonpositive
matrix_plot(-abs(P))
## Extending the color key to [-1,1] with darker color for |rho| >> 0
## Note: When specifying 'at', one most likely also wants to provide 'col.regions'
matrix_plot(P, at = seq(-1, 1, length.out = 200),
            col.regions = grey(c(seq(0, 1, length.out = 100), seq(1, 0, length.out = 100))))
}
\keyword{hplot}