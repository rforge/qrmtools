### Computing risk measures ####################################################

### 1 Value-at-Risk ############################################################

##' @title Nonparametric VaR estimator
##' @param x vector of losses
##' @param alpha confidence level
##' @param names see ?quantile
##' @param type 'type' used (1 = inverse of empirical df); see ?quantile
##' @param ... see ?quantile
##' @return nonparametric VaR_alpha estimate
##' @author Marius Hofert
##' @note We use the conservative type = 1 here as for sufficiently large alpha,
##'       type = 7 (quantile()'s default) would interpolate between the two
##'       largest losses (to be continuous) and thus return a(n even) smaller
##'       VaR_alpha estimate.
VaR_np <- function(x, alpha, names = FALSE, type = 1, ...)
    quantile(x, probs = alpha, names = names, type = type, ...) # vectorized in x and alpha

##' @title Value-at-Risk for normal and t distributions
##' @param alpha confidence level
##' @param mu location
##' @param sigma scale
##' @param df degrees of freedom; Inf for the normal distribution
##' @return Value-at-Risk
##' @author Marius Hofert
VaR_t <- function(alpha, mu = 0, sigma = 1, df = Inf)
{
    stopifnot(0 <= alpha, alpha <= 1, sigma > 0, df > 0)
    mu + sigma * if(identical(df, Inf)) qnorm(alpha) else qt(alpha, df = df)
}

##' @title Value-at-Risk for the Pareto distribution
##' @param alpha confidence level
##' @param theta Pareto parameter
##' @param kappa Pareto parameter
##' @return Value-at-Risk
##' @author Marius Hofert
VaR_Par <- function(alpha, theta, kappa = 1) qPar(alpha, theta = theta, kappa = kappa)


### 2 Expected shortfall #######################################################

##' @title Nonparametric expected shortfall estimator
##' @param x vector of losses
##' @param alpha confidence level
##' @param method method
##' @param verbose logical indicating whether verbose output is provided in
##'        case the mean is taken over (too) few losses
##' @param ... additional arguments passed VaR_np()
##' @return nonparametric ES_alpha estimate (derived under the assumption of continuity)
##' @author Marius Hofert
##' @note - Vectorized in x and alpha
##'       - ">" : Mathematically correct for discrete dfs, but
##'               produces NaN for alpha > (n-1)/n (=> F^-(alpha) = x_{(n)} but
##'               there is no loss strictly beyond x_{(n)})
##'         ">=": mean() will always include the largest loss (so no NaN appears),
##'               but might be computed just based on this one loss.
ES_np <- function(x, alpha, method = c(">", ">="), verbose = FALSE, ...)
{
    stopifnot(0 < alpha, alpha < 1)
    VaR <- VaR_np(x, alpha = alpha, ...) # length(alpha)-vector
    method <- match.arg(method)
    vapply(VaR, function(v) { # v = VaR value for one alpha
        ind <- if(method == ">") x > v else x >= v
        if(verbose) {
            num <- sum(ind)
            if(num == 0) {
                warning("No loss ",method," VaR; NaN returned instead")
            } else if(num == 1){
                warning("Only ",num," loss ",method," VaR")
            } else if(num <= 5) {
                warning("Only ",num," losses ",method," VaR")
            }
        }
        mean(x[ind]) # mean over all losses >(=) VaR
    }, NA_real_)
}

##' @title Expected shortfall for normal and t distributions
##' @param alpha confidence level
##' @param mu location
##' @param sigma scale
##' @param df degrees of freedom; Inf for the normal distribution
##' @return Expected shortfall
##' @author Marius Hofert
ES_t <- function(alpha, mu = 0, sigma = 1, df = Inf)
{
    stopifnot(0 <= alpha, alpha <= 1, sigma > 0, df > 0)
    mu + (sigma/(1-alpha)) * if(identical(df, Inf)) dnorm(qnorm(alpha)) else
    dt(qt(alpha, df = df), df = df) * (df + qt(alpha, df = df)^2) / (df-1)
}

##' @title Expected shortfall for the Pareto distribution
##' @param alpha confidence level
##' @param theta Pareto parameter
##' @param kappa Pareto parameter
##' @return Expected shortfall
##' @author Marius Hofert
ES_Par <- function(alpha, theta, kappa = 1)
{
    stopifnot(0 <= alpha, alpha <= 1, theta > 1, kappa > 0)
    kappa * ((theta / (theta-1)) * (1-alpha)^(-1/theta) - 1)
}
