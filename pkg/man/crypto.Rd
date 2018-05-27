\name{crypto}
\alias{crypto}
\docType{data}
\title{Cryptocurrencies}
\description{
  Bitcoin, Ethereum, Litecoin and Ripple prices in USD from their
  first available date to 2017-12-31.
}
\usage{
data("crypto")
}
\format{
  \code{\link{xts}} object containing cryptocurrency prices in USD
  Bitcoin, Ethereum, Litecoin and Ripple until 2015-12-31.
}
\source{
  The data was obtained from Yahoo Finance on 2018-05-27
  via the function \code{get_data()} from \CRANpkg{qrmtools}.
%% Created via:
%% library(xts)
%% library(qrmtools)
%% ## Get Cryptocurrencies (1 unit in USD)
%% ## See https://www.visualcapitalist.com/wp-content/uploads/2017/09/bitcoin-ethereum-other-cryptocurrencies.html
%% end <- "2017-12-31"
%% BTC_USD <- get_data("BTC-USD", to = end) # Bitcoin
%% ETH_USD <- get_data("ETH-USD", to = end) # Ethereum
%% LTC_USD <- get_data("LTC-USD", to = end) # Litecoin
%% XRP_USD <- get_data("XRP-USD", to = end) # Ripple
%% crypto <- do.call(merge, list(BTC_USD, ETH_USD, LTC_USD, XRP_USD))
%% colnames(crypto) <- c("BTC", "ETH", "LTC", "XRP")
%% save(crypto, file = "crypto.rda", compress = "xz")
}
\author{Marius Hofert}
\examples{
data("crypto")
str(crypto)
library(xts)
plot.zoo(crypto, main = "Cryptocurrencies in USD", xlab = "Time")
}
\keyword{datasets}
