\name{plot_NA}
\alias{plot_NA}
\title{Graphical Tool for Visualizing NAs in a Data Set}
\description{
  Plot NAs in a data set.
}
\usage{
plot_NA(x, col = c("black", "white"), xlab = "Time", ylab = "Component",
        text = "Black: NA; White: Available data",
        side = 4, line = 1, adj = 0, ...)
}
\arguments{
  \item{x}{A matrix (ideally an \code{xts} object).}
  \item{col}{A bivariate vector containing the colors for missing and
    available data, respectively.}
  \item{xlab}{The x-axis label.}
  \item{ylab}{The y-axis label.}
  \item{text}{See \code{\link{mtext}()}.}
  \item{side}{See \code{\link{mtext}()}.}
  \item{line}{See \code{\link{mtext}()}.}
  \item{adj}{See \code{\link{mtext}()}.}
  \item{\dots}{Additional arguments passed to the underlying function
    \code{\link{image}()}.}
}
\value{
  \code{\link{invisible}()}.
}
\details{
  Indicate \code{\link{NA}}s in a data set.
}
\author{Marius Hofert}
\examples{
data(list = "SP500_const", package = "qrmdata")
plot_NA(SP500_const) # some seconds...
}
\keyword{hplot}