\name{plot_matrix}
\alias{plot_matrix}
\title{Graphical Tool for Visualizing Square Matrices}
\description{
  Plot square matrices.
}
\usage{
plot_matrix(x, xlab="Column", ylab="Row",
            scales=list(alternating=c(1,1), tck=c(1,0), x=list(rot=90)),
            at=seq(min(x, na.rm=TRUE), max(x, na.rm=TRUE), length.out=200),
            colorkey=list(at=at),
            col.regions=colorRampPalette(c("royalblue3", "white", "maroon3"),
                                         space="Lab")(200),
            ...)
}
\arguments{
  \item{x}{A square matrix.}
  \item{xlab}{The x-axis label.}
  \item{ylab}{The y-axis label.}
  \item{scales}{See \code{\link{levelplot}()}; if \code{\link{NULL}},
    labels and ticks are omitted.}
  \item{at}{See \code{\link{levelplot}()}.}
  \item{colorkey}{See \code{\link{levelplot}()}.}
  \item{col.regions}{See \code{\link{levelplot}()}.}
  \item{\dots}{Additional arguments passed to the underlying function
    \code{\link{levelplot}()}.}
}
\value{
  The plot, a Trellis object.
}
\details{
  Plot of a square matrix.
}
\author{Marius Hofert}
\examples{
library(qrmdata)
data(SP500_const)
X <- diff(log(SP500_const))[-1,]
P <- cor(X, use="na.or.complete")
plot_matrix(P, scales=NULL, at=seq(-1, 1, length.out=200))
plot_matrix(P, scales=NULL, at=seq(-1, 1, length.out=200),
            col.regions=grey(c(seq(0, 1, length.out=100), seq(1, 0, length.out=100))))
}
\keyword{hplot}