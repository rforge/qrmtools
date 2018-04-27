### Mean excess function and plot ##############################################

##' @title (Sample) Mean Excess Function
##' @param x numeric vector of data.
##' @param omit number of unique last observations to be omitted
##'        from the sorted data (as mean excess estimator becomes unreliable for
##'        these observations as thresholds)
##' @return Two-column matrix giving the sorted data without the omit-largest
##'         unique values (first column) and the corresponding values of the
##'         sample mean excess function (second column).
##' @author Marius Hofert
##' @note - Main idea:
##'         1)              e_n(u) = (sum_{i=1}^n (X_{(i)}-u) I(X_{(i)} > u)) / sum_{i=1}^n I(X_{(i)} > u)
##'            (*) => e_n(X_{(k)}) = sum_{i=k+1}^n (X_{(i)}-X_{(k)}) / (n-k)
##'            ... but (*) holds only if there are *no* ties.
##'         2) Make the data unique, then apply (*), then assign each duplicated value
##'            the corresponding value of the mean excess function computed via (*).
mean_excess <- function(x, omit = 3)
{
    stopifnot(omit >= 1) # last one has to be omitted (otherwise e_n(X_{(i)}) not defined)
    ## Sort data and make it unique
    x <- sort(as.numeric(x)) # sorted data
    runs <- rle(x) # runs (make sense here because of x is sorted)
    xu <- runs$values # unique x values (sorted)
    nu <- length(xu) # how many unique x values
    ## Omit the last so-many unique values
    nuo <- nu - omit # how many we consider
    if(nuo < 1)
        stop("Not enough unique x values for computing the sample mean excess function.")
    ## Compute e_n(X_{k}) for k = 1,..,n based on unique X values and
    ## n = nu as described above
    enu <- vapply(1:nuo, function(k) mean(xu[(k+1):nu]-xu[k]), NA_real_) # note: nuo = nu - omit <= nu - 1 => nuo + 1 <= nu
    ## Expand the 'enu' values (to plot all points according to their correct number of appearances)
    xu.times <- runs$lengths[1:nuo] # numbers of how often the first nuo-many unique values appear
    en <- rep(enu, times = xu.times) # expand (corresponding mean excesses e_n(X_{(i)}))
    ## Return
    cbind("data" = x[1:sum(xu.times)], mean.excess = en) # X_{(i)}'s up to the last unique 'omit'-many data points
}

##' @title (Sample) Mean Excess Plot
##' @param x numeric vector of data.
##' @param omit number of unique last observations to be omitted
##'        from the sorted data (as mean excess plot becomes unreliable for
##'        these observations as thresholds)
##' @param xlab x-axis label; see plot()
##' @param ylab y-axis label; see plot()
##' @param ... additional arguments passed to the underlying plot()
##' @return invisible()
##' @author Marius Hofert
##' @note - Note that QRM::MEplot() only considers *unique* values. In particular,
##'         each value in the plot then appears equally often (which is actually
##'         wrong when using alpha transparency). This is also visible for the
##'         Danish data: str(MEplot(danish, omit = 3)) => 1645 points and not 2164
##'         as we obtain.
##'       - evd::mrlplot(QRM::danish) essentially uses type = "l" with pointwise
##'         asymptotic CIs, but does not seem too useful. qqtest's idea for CIs
##'         would be good but then probably too slow.
mean_excess_plot <- function(x, omit = 3, xlab = "Threshold",
                             ylab = "Mean excess over threshold", ...)
{
    plot(mean_excess(x, omit = omit), xlab = xlab, ylab = ylab, ...)
    invisible()
}


### Plot of fitted GPD shape as a function of the threshold ####################

##' @title Fitted GPD Shape as a Function of the Threshold
##' @param x numeric vector of data
##' @param thresholds numeric vector of thresholds for which to fit a GPD to
##'        the excesses
##' @param conf.level confidence level of the confidence intervals
##' @param lines.args list of arguments passed to underlying lines() for
##'        drawing the confidence intervals
##' @param xlab2 label of the secondary x-axis
##' @param xlab x-axis label
##' @param ylab y-axis label
##' @param ... additional arguments passed to the underlying plot()
##' @return invisible()
##' @author Marius Hofert
GPD_shape_plot <- function(x, thresholds = seq(quantile(x, 0.5), quantile(x, 0.99), length.out = 33),
                           conf.level = 0.95, lines.args = list(lty = 2),
                           xlab2 = "Excesses", xlab = "Threshold",
                           ylab = "Estimated GPD shape parameter with confidence intervals", ...)
{
    ## Checks
    stopifnot(length(thresholds) >= 2, 0 <= conf.level, conf.level <= 1)
    ## Fit GPD models to the given thresholds
    x <- as.numeric(x)
    fits <- lapply(thresholds, function(u) fit_GPD_MLE(x[x > u] - u))
    ## Extract the fitted shape parameters and compute CIs
    xi <- sapply(fits, function(f) f$par[["shape"]])
    xi.SE <- sapply(fits, function(f) f$SE[["shape"]])
    q <- qnorm(1 - (1 - conf.level)/2)
    xi.CI.low <- xi - xi.SE * q
    xi.CI.up  <- xi + xi.SE * q
    ## Plot
    ylim <- range(xi, xi.CI.low, xi.CI.up)
    plot(thresholds, xi, type = "l", ylim = ylim,
         xlab = xlab, ylab = ylab, ...)
    do.call(lines, args = c(list(x = thresholds, y = xi.CI.low), lines.args))
    do.call(lines, args = c(list(x = thresholds, y = xi.CI.up),  lines.args))
    pu <- pretty(thresholds)
    axis(3, at = pu, labels = sapply(pu, function(u) sum(x > u)))
    mtext(xlab2, side = 3, line = 3)
    invisible()
}
