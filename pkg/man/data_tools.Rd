\name{data_tools}
\alias{get_data}
\title{Tools for Getting and Working with Data}
\description{
  Download (and possibly) merge data from freely available databases.
}
\usage{
get_data(x, from=NULL, to=NULL, src="yahoo", FUN=NULL, verbose=TRUE, ...)
}
\arguments{
  \item{x}{A vector of ticker symbols (e.g. \code{"^GSPC"} if
    \code{src="yahoo"} or \code{"EUR/USD"} if \code{src="oanda"}).}
  \item{from}{start date as a \code{Date} object or character string
    (in international date format \code{"yyyy-mm-dd"}); if \code{NULL},
    the earliest date with available data is picked.}
  \item{to}{end date as a \code{Date} object or character string
    (in international date format \code{"yyyy-mm-dd"}); if \code{NULL},
    the last date with available data is picked.}
  \item{src}{A character string specifying the data source
    (e.g. \code{"yahoo"} for stocks or \code{"oanda"} for FX data);
    see \code{\link{getSymbols}()}.}
  \item{FUN}{A \code{\link{function}} to be applied to the data before
    being returned. This can be
    \describe{
      \item{the identity}{if the data could not be retrieved (and is
	thus replaced by \code{\link{NA}});}
      \item{the given \code{FUN}}{if \code{FUN} has been provided;}
      \item{a useful default}{if \code{FUN=NULL}; the default uses
	the adjusted close price \code{\link{Ad}()} if
	\code{src="yahoo"},
	the close price \code{\link{Cl}()} if
	\code{src="google"} and the identity otherwise.}
    }}
  \item{verbose}{A \code{\link{logical}} indicating whether progress
    monitoring should be done.}
  \item{\dots}{Additional arguments passed to the underlying function
    \code{\link{getSymbols}()} from \pkg{quantmod}.}
}
\value{
  An \code{\link{xts}} object containing the data with column name(s)
  adjusted to be (close to) the ticker symbol.
}
\details{
  \code{FUN} is typically one of \pkg{quantmod}'s \code{\link{Op}},
  \code{\link{Hi}}, \code{\link{Lo}}, \code{\link{Cl}}, \code{\link{Vo}},
  \code{\link{Ad}} or one of the combined functions \code{\link{OpCl}},
  \code{\link{ClCl}}, \code{\link{HiCl}}, \code{\link{LoCl}}, \code{\link{LoHi}},
  \code{\link{OpHi}}, \code{\link{OpLo}}, \code{\link{OpOp}}.
}
\author{Marius Hofert}
\examples{
## Note: The following needs a working internet connection!
## Get stock and volatility data (all available trading days)
dat <- get_data(c("^GSPC", "^VIX"))
require(zoo)
plot.zoo(dat, screens=1, main="", xlab="Trading day", ylab="Value")
## Alternative: plot.xts() from xtsExtra

## Get FX data
## Note: get_data() internally recursively calls itself to download more than
##       the normal limit of 5y of daily FX data from OANDA
EUR_USD <- get_data("EUR/USD", from="2000-01-01", to="2015-12-11", src="oanda")
plot(EUR_USD, xlab="Trading day", ylab="EUR/USD")
## Getting an empty xts object (FX data in a period not available);
\dontrun{ # Note: more time-consuming
    EUR_USD <- get_data("EUR/USD", from="1990-01-01", to="1993-12-31", src="oanda") # < 5y
    EUR_USD <- get_data("EUR/USD", from="1990-01-01", to="1995-12-31", src="oanda") # >= 5y
}
}
\keyword{manip}