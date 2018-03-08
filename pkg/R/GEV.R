### GEV(xi, mu, sigma) distribution ############################################

##' @title Density of the GEV(xi, mu, sigma) distribution
##' @param x evaluation points
##' @param xi parameter xi (real)
##' @param mu parameter mu (real)
##' @param sigma parameter sigma (also real here; density is 0 if sigma <= 0)
##' @param log logical indicating whether the log density is computed
##' @return density of the GEV(xi, mu, sigma) distribution
##' @author Marius Hofert
dGEV <- function(x, xi, mu = 0, sigma = 1, log = FALSE)
{
    if(sigma <= 0)
        return(if(log) -Inf else 0)
    y <- (x - mu) / sigma
    if(xi == 0) { # xi == 0
        res <- -log(sigma) - (y + exp(-y))
    } else { # xi != 0
        res <- rep(-Inf, length(y)) # correctly extend log-density
        xiy <- xi * y
        ii <- 1 + xiy > 0
        res[ii] <- -log(sigma) + (-1/xi - 1) * log1p(xiy[ii]) - (1 + xiy[ii])^(-1/xi)
    }
    if(log) res else exp(res)
}

##' @title Distribution function of the GEV(xi, mu, sigma) distribution (vectorized in q)
##' @param q quantile
##' @param xi parameter xi
##' @param mu parameter mu
##' @param sigma parameter sigma
##' @param lower.tail logical indicating whether lower/upper tail is used
##' @param log.p logical indicating whether probabilities are given as log()
##' @return distribution function of the GEV(xi, mu, sigma) distribution
##' @author Marius Hofert
pGEV <- function(q, xi, mu = 0, sigma = 1, lower.tail = TRUE, log.p = FALSE)
{
    stopifnot(sigma > 0)
    y <- (q-mu)/sigma
    if(xi == 0) { # xi == 0
        if(lower.tail)
            if(log.p) -exp(-y) else exp(-exp(-y))
        else if(log.p) log1p(-exp(-exp(-y))) else 1-exp(-exp(-y))
    } else { # xi != 0
        xiy <- pmax(xi*y, -1) # see dGEV()
        if(lower.tail) {
            if(log.p) {
                -(1+xiy)^(-1/xi) # log H
            } else {
                exp(-(1+xiy)^(-1/xi)) # H
            }
        } else {
            if(log.p) {
                log1p(-exp(-(1+xiy)^(-1/xi))) # log(bar{H})
            } else {
                1-exp(-(1+xiy)^(-1/xi)) # bar{H}
            }
        }
        ## Formerly (also works):
        ## xiy <- xi*y
        ## ii <- 1+xiy > 0
        ## if(lower.tail) {
        ##     if(log.p) {
        ##         res <- if(xi < 0) rep(0, length(q)) else rep(-Inf, length(q)) # log('see case below')
        ##         res[ii] <- -(1+xiy[ii])^(-1/xi) # log H
        ##     } else {
        ##         res <- if(xi < 0) rep(1, length(q)) else rep(0, length(q)) # correctly extend; think of case '1+xi*x = 0' in both cases
        ##         res[ii] <- exp(-(1+xiy[ii])^(-1/xi)) # H
        ##     }
        ## } else {
        ##     if(log.p) {
        ##         res <- if(xi < 0) rep(-Inf, length(q)) else rep(0, length(q)) # log('see case below')
        ##         res[ii] <- log1p(-exp(-(1+xiy[ii])^(-1/xi))) # log(bar{H})
        ##     } else {
        ##         res <- if(xi < 0) rep(0, length(q)) else rep(1, length(q)) # 1 - 'case above'
        ##         res[ii] <- 1-exp(-(1+xiy[ii])^(-1/xi)) # bar{H}
        ##     }
        ## }
        ## res
    }
}

##' @title Quantile function of the GEV(xi, mu, sigma) distribution (vectorized in p)
##' @param p probability
##' @param xi parameter xi
##' @param mu parameter mu
##' @param sigma parameter sigma
##' @param lower.tail logical indicating whether lower/upper tail is used
##' @param log.p logical indicating whether probabilities are given as log()
##' @return quantile function of the GEV(xi, mu, sigma) distribution
##' @author Marius Hofert
qGEV <- function(p, xi, mu = 0, sigma = 1, lower.tail = TRUE, log.p = FALSE)
{
    stopifnot(sigma > 0)
    p <- if(log.p) pmin(p, 0) else pmin(pmax(p, 0), 1) # correctly extend
    if(xi == 0) { # xi == 0
        res <- if(lower.tail)
                   if(log.p) -log(-p) else -log(-log(p))
               else if(log.p) -log(-log1p(-exp(p))) else -log(-log1p(-p))
    } else { # xi != 0
          res <- if(lower.tail)
                     if(log.p) ((-p)^(-xi)-1)/xi else ((-log(p))^(-xi)-1)/xi
                 else if(log.p) ((-log1p(-exp(p)))^(-xi)-1)/xi else ((-log1p(-p))^(-xi)-1)/xi
          res <- if(xi < 0) pmin(res, -1/xi) else pmax(res, -1/xi)
    }
    mu + sigma * res
}

##' @title Generating random numbers from a GEV(xi, mu, sigma) distribution
##' @param n sample size n
##' @param xi parameter xi
##' @param mu parameter mu
##' @param sigma parameter sigma
##' @return n-vector containing GEV(xi, mu, sigma) random variates
##' @author Marius Hofert
rGEV <- function(n, xi, mu = 0, sigma = 1)
    qGEV(runif(n), xi = xi, mu = mu, sigma = sigma)


### Fitting GEV(xi, mu, sigma) #################################################

## ## Other packages
## - QRM:
##   + fit.GEV() uses hard-coded value for xi and mu and sigma from case xi = 0
##     (see below)
## - Renext:
##   + also based on optim()
##   + fGEV.MAX() -> calls parIni.MAX(); detailed in "Renext Computing Details"
##     (but not available)
##   + some sort of weird procedure based on a regression idea of some sort
##   + good 'caution': log-lik = Inf for xi <= -1 but optimization
##     done in unconstrained way. Typically xi > -1, but if xi <= -1, then
##     meaningless.
##     For xi < -0.5, treat with care.
## - evd:
##   + also based on optim()
##   + fextreme(): requires 'start' to be provided
##   + fgev(): requires 'start' to be provided
## - evir:
##   + ./R -> bmax.R -> gev(): hardcoded xi = 0.1 as 'QRM'
## - extRemes:
##   + fevd(): huge, uses optim() with 'init.pars', can provide 'initial' or
##     determine them
##   + 'initial' is determine (on 'find.init') via L moments or moments and
##     hard-coded (xi = 0.01, mu = 0, sig = 1) on failure of both methods
## - fExtremes:
##   + gevFit() -> .gevFit() -> .gevmleFit(): uses optim(), hardcoded initial
##     values as 'QRM' (refers to evir for that)
## - ismev:
##   + gev.fit(): uses optim(), hardcoded initial values as in 'QRM'
## - lmom:
##   + pelp(): based on optim(), 'start' needs to be provided
## - texmex:
##   + evm() -> evm.default() -> evmFit(): based on optim(); uses
##     'family$start(data)' if start not provided (see also .pdf)
##   + 'start()' is a function and for the GEV found in ./R/gev.R
##   + start() returns a longer vector (unclear why) but seems to use mean()
##     as initial value for mu and log(IQR(data)/2) as initial value for xi;
##     quite unclear
##   + points out that MLE will often fail with small sample size

##' @title Computing Initial Values for MLE of the GEV
##' @param x numeric vector of data. In the block maxima method, these are the
##'        block maxima (based on block size n).
##' @return numeric(3) specifying the initial values for xi, mu, sigma.
##' @author Marius Hofert
##' @note Used to have a nice quantile-matching idea here, but this (as any
##'       other method) could lead to 1 + xi (x - mu) / sigma <= 0 for some
##'       data x and thus log-likelihood = -Inf, so optim() stops with error:
##'       "Error in optim(init, fn = function(param) logLik_GEV(param, x = x),
##'        hessian = estimate.cov: function cannot be evaluated at initial
##'        parameters"
fit_GEV_init <- function(x)
{
    ## Idea: - Use xi = 0 here => mu and sigma explicit in terms of mean and
    ##         variance (method-of-moment estimator for mu, sigma).
    ##         This guarantees that 1 + xi (x - mu) / sigma > 0, so a finite
    ##         log-likelihood.
    ##       - As in 'QRM', 'evir', 'fExtremes', 'ismev'
    sig.init <- sqrt(6 * var(x)) / pi # var for xi = 0 is (\sigma\pi)^2 / 6 => \sigma
    c(0, mean(x) - 0.5772157 * sig.init, sig.init) # mean for xi is mu + sig * gamma => mu; gamma = -digamma(1) = Euler--Mascheroni constant
}

##' @title Log-likelihood of the GEV
##' @param param numeric(3) giving xi, mu, sigma (all real here; if sigma <= 0
##'        dGEV(, log = TRUE) returns -Inf)
##' @param x numeric vector of data. In the block maxima method, these are the
##'        block maxima (based on block size n)
##' @return log-likelihood at xi, mu, sigma
##' @author Marius Hofert
logLik_GEV <- function(param, x)
    sum(dGEV(x, xi = param[1], mu = param[2], sigma = param[3], log = TRUE))

##' @title MLE for GEV Parameters
##' @param x numeric vector of data. In the block maxima method, these are the
##'        block maxima (based on block size n)
##' @param init numeric(3) giving the initial values for xi, mu, sigma; if NULL,
##'        determined by a procedure to find initial values
##' @param estimate.cov logical indicating whether the asymptotic covariance
##'        matrix of the parameter estimators is to be estimated
##'        (inverse of observed Fisher information (negative Hessian
##'        of log-likelihood evaluated at MLE))
##' @param control see ?optim
##' @param ... additional arguments passed to the underlying optim()
##' @return list with the return value of optim() with the estimated asymptotic
##'         covariance matrix of the parameter estimators
##'         (Cramer--Rao bound = inverse of Fisher information = negative Hessian)
##'         of the parameter estimators appended if estimate.cov.
##' @author Marius Hofert
##' @note - similar to copula:::fitCopula.ml()
##'       - careful for xi <= -0.5 (very short, bounded upper tail):
##'         MLE doesn't have standard asymptotic properties.
fit_GEV <- function(x, init = NULL, estimate.cov = TRUE, control = list(), ...)
{
    ## Checks
    stopifnot(is.numeric(x), is.null(init) || (length(init) == 3),
              is.logical(estimate.cov), is.list(control))
    ## Compute initial values
    if(is.null(init)) init <- fit_GEV_init(x)
    ## Fit
    control <- c(as.list(control), fnscale = -1) # maximization (overwrites possible additionally passed 'fnscale')
    fit <- optim(init, fn = function(param) logLik_GEV(param, x = x),
                 hessian = estimate.cov, control = control, ...)
    ## Estimate of the asymptotic covariance matrix of the parameter estimators
    ## Note: manually; see QRM::fit.GEV and also copula:::fitCopula.ml():
    ##       fisher <- -hessian(logLik_GEV(fit$par, x = x)) # observed Fisher information (estimate of Fisher information)
    ##       Cov <- solve(fisher)
    ##       std.err <- sqrt(diag(Cov)) # see also ?fit_GEV
    Cov <- if(estimate.cov) {
               negHessianInv <- catch(solve(-fit$hessian))
               if(is(negHessianInv, "error")) {
                   warning("Hessian matrix not invertible: ", negHessianInv$error)
                   matrix(NA_real_, 0, 0)
               } else negHessianInv$value # result on warning or on 'worked'
           } else matrix(NA_real_, 0, 0)
    ## Return (could create an object here)
    if(estimate.cov) c(fit, list(Cov = Cov)) else fit
}
