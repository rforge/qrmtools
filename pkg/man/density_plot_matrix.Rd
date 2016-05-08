\name{density_plot_matrix}
\alias{density_plot_matrix}
\title{Density Plot of the Values from a Lower Triangular Matrix}
\description{
  Density plot of all values in the lower triangular part of a matrix.
}
\usage{
density_plot_matrix(x, xlab = "Entries in the lower triangular matrix",
                    main = "", text = NULL, side = 4, line = 1, adj = 0, ...)
}
\arguments{
  \item{x}{A \code{\link{matrix}}-like object.}
  \item{xlab}{The x-axis label.}
  \item{main}{The title.}
  \item{text}{See \code{\link{mtext}()}. The \code{text = ""}, it is omitted.}
  \item{side}{See \code{\link{mtext}()}.}
  \item{line}{See \code{\link{mtext}()}.}
  \item{adj}{See \code{\link{mtext}()}.}
  \item{\dots}{Additional arguments passed to the underlying \code{\link{plot}()}.}
}
\value{
  \code{\link{invisible}()}.
}
\details{
  \code{density_plot_matrix()} is typically used for symmetric matrices
  (like correlation matrices, matrices of pairwise Kendall's tau or tail
  dependence parameters) to check the distribution of their off-diagonal
  entries.
}
\author{Marius Hofert}
\examples{
data(list = "SP500_const", package = "qrmdata")
X <- log_returns(SP500_const)
P <- cor(X, use = "na.or.complete")
density_plot_matrix(P)
}
\keyword{hplot}