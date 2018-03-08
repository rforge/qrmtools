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
##' @param method character indicating the method to be used:
##'        - zero.xi: uses xi = 0 and method-of-moment estimators based on that case
##'        - quantile.matching: chooses xi, mu, sigma such that the empirical
##'          1/4-, 1/2-, 3/4-quantile are matched
##'        - prob.weighted.moments: Matching probability-weighted moments;
##'          see Hosking, Wallis, Wood (1985)
##' @return numeric(3) specifying the initial values for xi, mu, sigma.
##' @author Marius Hofert
##' @note Used to have a nice quantile-matching idea here, but this (as any
##'       other method) could lead to 1 + xi (x - mu) / sigma <= 0 for some
##'       data x and thus log-likelihood = -Inf, so optim() stops with error:
##'       "Error in optim(init, fn = function(param) logLik_GEV(param, x = x),
##'        hessian = estimate.cov: function cannot be evaluated at initial
##'        parameters"
fit_GEV_init <- function(x, method = c("zero.xi", "quantile.matching", "prob.weighted.moments"))
{
    method <- match.arg(method)
    switch(method,
    "zero.xi" = {
        ## Idea: - Use xi = 0 here => mu and sigma explicit in terms of mean and
        ##         variance (method-of-moment estimator for mu, sigma).
        ##         This guarantees that 1 + xi (x - mu) / sigma > 0, so a finite
        ##         log-likelihood.
        ##       - As in 'QRM', 'evir', 'fExtremes', 'ismev'
        sig.init <- sqrt(6 * var(x)) / pi # var for xi = 0 is (\sigma\pi)^2 / 6 => \sigma
        c(0, mean(x) - 0.5772157 * sig.init, sig.init) # mean for xi is mu + sig * gamma => mu; gamma = -digamma(1) = Euler--Mascheroni constant
    },
    "prob.weighted.moments" = {
        ## b_r = estimator / sample version of M_{1,r,0}
        x. <- sort(x)
        n <- length(x.)
        i <- 1:n
        b0 <- mean(x)
        b1 <- mean((i-1)/(n-1) * x.)
        b2 <- mean((i-1)*(i-2)/((n-1)*(n-2)) * x.)
        ## Match probability-weighted moments
        y <- (3*b2-b0) / (2*b1-b0) # evaluation point of h^{-1}
        y <- max(y, 1) # sanity
        h.prime.m1 <- (log(3) * 3^(-1) * (2^(-1) - 1) - log(2) * 2^(-1) * (3^(-1) - 1)) / (2^(-1) - 1)^2 # h'(-1)
        xi.init <- if(y <= 4/3 + h.prime.m1 * (1.455495 + 1)) { # y <= h approximation at 1.455495
                       (y - 4/3) / h.prime.m1 - 1 # invert tangent in -1
                   } else { # invert (3/2)^xi
                       log(y) / log(3/2)
                   }
        if(xi.init >= 1) # should not happen (according to Hosking, Wallis, Wood (1985))
            xi.init <- 0.95 # set here due to gamma() function call below
        sig.init <- if(xi.init == 0) {
                        (2 * b1 - b0) / log(2)
                    } else {
                        (2 * b1 - b0) * xi.init / (gamma(1 - xi.init) * (2^xi.init - 1))
                    }
        mu.init <- if(xi.init == 0) {
                       b0 - sig.init * 0.5772157 # ... Euler--Mascheroni constant
                   } else {
                       b0 - sig.init / xi.init * (gamma(1-xi.init) - 1)
                   }
        ## Force initial values to produce non-zero density
        sig.bound <- max(-xi.init * (x - mu.init))
        if(sig.init <= sig.bound)
            sig.init <- 1.05 * sig.bound # blow up by 5%
        ## Return
        c(xi.init, mu.init, sig.init)
    },
    "quantile.matching" = {
        p <- c(0.25, 0.5, 0.75) # probabilities whose quantiles are matched (could be made an argument)
        cutoff <- 3 # value after which exp(-x) is truncated to 0 (could be made an argument)
        q <- quantile(x, probs = p, names = FALSE) # empirical p-quantiles
        if(length(unique(q)) < 3) {
            fit_GEV_init(x, method = "zero.xi")
        }
        q.diff <- diff(q)
        y <- q.diff[1]/q.diff[2] # (q[2] - q[1]) / (q[3] - q[2]) = (H^{-1}(p2) - H^{-1}(p1)) / (H^{-1}(p3) - H^{-1}(p2))
        l <- log(-log(p)) # decreasing in p
        a <- rev(diff(rev(l))) # l[1] - l[2]; l[2] - l[3]
        a. <- a[1]/a[2]
        ## Initial value for xi
        xi.init <- if(y < 1/expm1(cutoff/a.)) {
                       log1p(1/y)/a[2]
                   } else if(y <= a.) {
                       m1 <- a[2]/cutoff * log(a./(expm1(a.*cutoff)))
                       log(y/a.) / m1
                   } else if(y <= expm1(a.*cutoff)) {
                       m2 <- -a[1]/cutoff * log(a.*expm1(cutoff/a.))
                       log(y/a.) / m2
                   } else -log1p(y)/a[1]
        ## Initial value for sigma
        sig.init <- if(xi.init == 0) {
                        q.diff[1] / (-l[2] + l[1])
                    } else {
                        xi.init * q.diff[1] / ((-log(p[2]))^(-xi.init) - (-log(p[1]))^(-xi.init))
                    }
        ## Initial value for mu
        mu.init <- if(xi.init == 0) {
                       q[2] + sig.init * l[2]
                   } else {
                       q[2] - sig.init/xi.init * ((-log(p[2]))^(-xi.init)-1)
                   }
        ## Force initial values to produce non-zero density
        sig.bound <- max(-xi.init * (x - mu.init))
        if(sig.init <= sig.bound)
            sig.init <- 1.05 * sig.bound # blow up by 5%

        ## Return
        c(xi.init, mu.init, sig.init)
    },
    stop("Wrong 'method'"))
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
