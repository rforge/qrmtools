\name{SnP_500_constituents}
\alias{SnP_500_constituents}
\docType{data}
\title{S&P 500 Constituents Data}
\description{
  Data set consisting of daily close prices from 2000-01-01 to
  2015-11-13 of the 505 constituents of the S\&P 500 as of 2015-11-13.
}
\usage{
data("SnP_500_constituents")
}
\format{
  A 505-column \code{xts} object (MMM, ABT, ABBV, ..., ZTS) containing
  merged historical time series data from from 2000-01-01 to
  2015-11-13 of each of the S\&P 500 constituents as of
  2015-11-13. Note that this is not the same as the
  real historical S\&P 500 constituents; the latter data is only
  available through specialized databases. % e.g., WRDS; see https://bizlib247.wordpress.com/2014/10/06/historical-index-constituents-e-g-sp-500/
}
\source{
  The data was drawn from Yahoo! Finance
  (\url{http://finance.yahoo.com/}) on 2015-11-13 based on the tickers
  found on Wikipedia (\url{https://en.wikipedia.org/wiki/List_of_S\%26P_500_companies}).
}
\examples{
data("SnP_500_constituents") # use as.data.frame(SnP_500_constituents) for converting it to a data.frame

## Plot NAs
dat <- SnP_500_constituents
z <- is.na(dat)
image(seq_len(nrow(dat)), seq_len(ncol(dat)), z, col=c("white", "black"),
      xlab="Time t", ylab="Stock j",
      main="Plot of NAs among the S&P 500 constituents")
mtext("Black: NA; White: Data available", side=4, line=1, adj=0)
}
\keyword{datasets}
