### Tools for working with data sets ###########################################

##' @title Download data via quantmod's getSymbols()
##' @param x A vector of, for example, ticker symbols (if src="yahoo") or
##'        "EUR/USD" if (src="oanda")
##' @param from start date as character string (e.g. 2015-01-01); if NULL,
##'        the earliest available date is picked
##' @param to end date; today unless otherwise specified
##' @param src The source of the data ("yahoo", "oanda", "google", etc.)
##' @param FUN A function to apply to the downloaded data:
##'        - if data is NA (could not be retrieved): none
##'        - if provided: the given function
##'        - if not provided: Ad() if src="yahoo"; Cl() if src="google";
##'          none otherwise (e.g. if src="oanda")
##' @param verbose A logical indicating whether progress monitoring is done
##' @param ... Additional arguments passed to getSymbols()
##' @return (n, d)-matrix of data (an xts object)
##' @author Marius Hofert
##' @note - One could do...
##'         nenv <- new.env()
##'         getSymbols(c("^GSPC", "AAPL"), env=nenv, from="2015-11-28", to="2015-12-07")
##'         ... but we would not get progress monitoring then
##'       - In case of failure to get the data, note that we don't need to
##'         throw a warning, as getSymbols() produces them (unfortunately even
##'         if 'warnings=FALSE')
get_data <- function(x, from=NULL, to=NULL, src="yahoo", FUN=NULL,
                     verbose=TRUE, ...)
{
    ## Checking
    if(is.data.frame(x)) x <- as.matrix(x)
    if(is.matrix(x)) x <- as.vector(x)
    stopifnot((d <- length(x)) >= 1)
    if(is.null(from)) from <- "1900-01-01" # to get all data available
    if(is.null(to)) to <- as.character(Sys.Date()) # to get all data available
    stopifnot(is.character(from) || inherits(from, "Date"),
              is.character(to) || inherits(to, "Date"),
              is.character(src), is.logical(verbose))

    ## Distinguish univariate/multivariate data
    res <- if(d == 1) {
        ## Get data
        start <- as.Date(from)
        end <- as.Date(to)
        time.diff <- difftime(end, start, units="days")[[1]] # in days
        num.blocks.5y <- time.diff/(5*365) # number of blocks of 5y to fit in [start, end]
        if(src=="oanda" && num.blocks.5y > 1) { # get the data in blocks (5y max)
            periods <- seq(start, end, length.out=ceiling(num.blocks.5y)+1) # start,..., end (>= 2 dates of class Date)
            num.periods <- length(periods)-1 # >= 1 periods to get data from
            dat <- NULL
            for(k in seq_len(num.periods)) {
                from. <- as.character(periods[k])
                to. <- as.character(if(k==num.periods) periods[k+1] else periods[k+1]-1)
                if(verbose) cat("Getting", x, "from", from., "to", to., "\n")
                dat. <- get_data(x, from=from., to=to., src=src, FUN=FUN, verbose=verbose, ...)
                if(!(length(dat.)==1 && is.na(dat.)))
                    dat <- rbind(dat, dat.) # exclude NAs (if no data available for that time period)
            }
        } else { # get the data all at once
            if(!is.character(from)) from <- as.character(from)
            if(!is.character(to)) to <- as.character(to)
            dat <- tryCatch(getSymbols(x, from=from, to=to, src=src, auto.assign=FALSE, ...),
                            error=function(e) e)
        }
        ## Return NA in case of an error
        if(is(dat, "simpleError")) {
            NA
        } else { # getting the data worked fine
            ## Apply FUN
            if(is.null(FUN)) { # if not given, apply a useful default
                if(src=="yahoo") {
                    Ad(dat) # does a grep on colnames(x)
                } else if(src=="google") {
                    Cl(dat) # no Ad() available, use Cl()
                } else dat # nothing (important for src="oanda")
            } else { # FUN given
                stopifnot(is.function(FUN))
                dat. <- FUN(dat)
                if(ncol(dat.)!=1)
                    stop("'FUN' has to return a single time series")
                dat.
            }
        }
    } else { # d > 1
        res <- vector("list", length=d)
        for(j in 1:d) {
            if(verbose) cat("Getting", x[j], "\n")
            res[[j]] <- get_data(x[j], from=from, to=to, src=src, FUN=FUN, verbose=verbose, ...)
        }
        do.call(merge, res)
    }

    ## Use ticker symbols as column names
    names(res) <- x
    res
}
