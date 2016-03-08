\name{plot_NA}
\alias{plot_NA}
\title{Graphical Tool for Visualizing NAs in a Data Set}
\description{
  Plot NAs in a data set.
}
\usage{
plot_NA(x, xlab="Time", ylab="Component", main="",
        text="Black: NA; White: Data available", side=4, line=1, adj=0, ...)
}
\arguments{
  \item{x}{A matrix (ideally an \code{xts} object).}
  \item{xlab}{The x-axis label.}
  \item{ylab}{The y-axis label.}
  \item{main}{A title.}
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
  Plot NAs (in black) in a data set; white stands for \dQuote{data available}.
}
\author{Marius Hofert}
\examples{
library(qrmdata)
data(SP500_const)
plot_NA(SP500_const)
}
\keyword{hplot}