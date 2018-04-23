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
##'        - xi0: uses xi = 0 and method-of-moment estimators based on that case
##'        - PWM: Matching probability-weighted moments;
##'          see Hosking, Wallis, Wood (1985)
##'        - quant.match: chooses xi, mu, sigma such that the empirical
##'          1/4-, 1/2-, 3/4-quantile are matched
##' @return numeric(3) specifying the initial values for xi, mu, sigma.
##' @author Marius Hofert
##' @note 1) No method except 'xi0' guarantees that 1 + xi (x - mu) / sigma > 0
##'          and thus that the log-likelihood > -Inf for all data points; if
##'          log-likelihood == 0, optim() stops with error:
##'          "Error in optim(init, fn = function(param) logLik_GEV(param, x = x),
##'          hessian = estimate.cov: function cannot be evaluated at initial
##'          parameters"
##'          => below we guarantee to get a finite log-likelihood
##'             (by doubling sigma) no matter what the input
##'       2) Other options as default method are 'xi0' with xi = 0 (exact) but
##'          method = "BFGS", or 'PWM"
fit_GEV_init <- function(x, method = c("xi0", "PWM", "quant.match"))
{
    method <- match.arg(method)
    init <- switch(method,
    "xi0" = {
        ## Idea: - Use xi = 0 here => mu and sigma explicit in terms of mean and
        ##         variance (method-of-moment estimator for mu, sigma).
        ##         This guarantees that 1 + xi (x - mu) / sigma > 0, so a finite
        ##         log-likelihood.
        ##       - As in 'QRM', 'evir', 'fExtremes', 'ismev'
        ##       - Note that xi = 0 can fail for method = "Nelder-Mead"
        ##         (seen for the Black Monday example)
        ##         => .Machine$double.eps works
        sig.hat <- sqrt(6 * var(x)) / pi # var for xi = 0 is (\sigma\pi)^2 / 6 => \sigma
        c(.Machine$double.eps, mean(x) - 0.5772157 * sig.hat, sig.hat) # mean for xi is mu + sig * gamma => mu; gamma = -digamma(1) = Euler--Mascheroni constant
    },
    "PWM" = {
        ## b_r = estimator (unbiased according to Landwehr and Wallis (1979)) / a sample version of M_{1,r,0} = E(X F(X)^r)
        x. <- sort(as.numeric(x)) # Caution! If is.xts(x), then sort won't do anything!
        n <- length(x.)
        k <- 1:n
        b0 <- mean(x)
        b1 <- mean(x. * (k-1)/(n-1))
        b2 <- mean(x. * (k-1)*(k-2)/((n-1)*(n-2)))
        ## Match probability-weighted moments
        y <- (3*b2-b0) / (2*b1-b0) # evaluation point of h^{-1}
        y <- max(y, 1) # sanity
        xi.hat <- if(y <= (2/3) * (2-log(3/4)) * (1.455495 + 1)) { # y <= h approximation at 1.455495
                       (2-(3/2)*y) / log(3/4) - 1 # invert tangent in -1
                   } else { # invert (3/2)^xi
                       log(y) / log(3/2)
                   }
        if(xi.hat >= 1) # should not happen (according to Hosking, Wallis, Wood (1985))
            xi.hat <- 0.95 # set here due to gamma() function call below
        sig.hat <- if(xi.hat == 0) {
                        (2 * b1 - b0) / log(2)
                    } else {
                        (2 * b1 - b0) * xi.hat / (gamma(1 - xi.hat) * (2^xi.hat - 1))
                    }
        mu.hat <- if(xi.hat == 0) {
                       b0 - sig.hat * 0.5772157 # ... Euler--Mascheroni constant
                   } else {
                       b0 - sig.hat / xi.hat * (gamma(1-xi.hat) - 1)
                   }
        ## Return
        c(xi.hat, mu.hat, sig.hat)
    },
    "quant.match" = {
        p <- c(0.25, 0.5, 0.75) # probabilities whose quantiles are matched (could be made an argument)
        cutoff <- 3 # value after which exp(-x) is truncated to 0 (could be made an argument)
        q <- quantile(x, probs = p, names = FALSE) # empirical p-quantiles
        if(length(unique(q)) < 3) {
            fit_GEV_init(x, method = "xi0")
        }
        q.diff <- diff(q)
        y <- q.diff[1]/q.diff[2] # (q[2] - q[1]) / (q[3] - q[2]) = (H^{-1}(p2) - H^{-1}(p1)) / (H^{-1}(p3) - H^{-1}(p2))
        l <- log(-log(p)) # decreasing in p
        a <- rev(diff(rev(l))) # l[1] - l[2]; l[2] - l[3]
        a. <- a[1]/a[2]
        ## Initial value for xi
        xi.hat <- if(y < 1/expm1(cutoff/a.)) {
                       log1p(1/y)/a[2]
                   } else if(y <= a.) {
                       m1 <- a[2]/cutoff * log(a./(expm1(a.*cutoff)))
                       log(y/a.) / m1
                   } else if(y <= expm1(a.*cutoff)) {
                       m2 <- -a[1]/cutoff * log(a.*expm1(cutoff/a.))
                       log(y/a.) / m2
                   } else -log1p(y)/a[1]
        ## Initial value for sigma
        sig.hat <- if(xi.hat == 0) {
                        q.diff[1] / (-l[2] + l[1])
                    } else {
                        xi.hat * q.diff[1] / ((-log(p[2]))^(-xi.hat) - (-log(p[1]))^(-xi.hat))
                    }
        ## Initial value for mu
        mu.hat <- if(xi.hat == 0) {
                       q[2] + sig.hat * l[2]
                   } else {
                       q[2] - sig.hat/xi.hat * ((-log(p[2]))^(-xi.hat)-1)
                   }
        ## Return
        c(xi.hat, mu.hat, sig.hat)
    },
    stop("Wrong 'method'"))
    ## Force initial values to produce finite log-likelihood
    ## Formerly:
    ## if(method != "xi0") { # only for 'xi0', a non-zero density is guaranteed
    ##     sig.bound <- max(-init[1] * (x - init[2]))
    ##     if(init[3] <= sig.bound)
    ##         init[3] <- 1.05 * sig.bound # blow up by 5%
    ## }
    while(!is.finite(logLik_GEV(init, x = x)) && is.finite(init[3])) init[3] <- init[3] * 2
    ## Note: if !is.finite(init[3]), there's nothing we can do...
    ## Return
    init
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
fit_GEV_MLE <- function(x, init = NULL, estimate.cov = TRUE, control = list(), ...)
{
    ## Checks
    stopifnot(is.numeric(x), is.null(init) || is.character(init) || (length(init) == 3),
              is.logical(estimate.cov), is.list(control))
    ## Compute initial values
    if(is.null(init)) {
        init <- fit_GEV_init(x)
    } else if(is.character(init)) { # easter egg
        init <- fit_GEV_init(x, method = match.arg(init, choices = eval(formals(fit_GEV_init)$method)))
    }
    ## Fit
    control <- c(as.list(control), fnscale = -1) # maximization (overwrites possible additionally passed 'fnscale')
    fit <- optim(init, fn = function(param) logLik_GEV(param, x = x),
                 hessian = estimate.cov, control = control, ...)
    names(fit$par) <- c("xi", "mu", "sigma")
    ## Estimate of the asymptotic covariance matrix and standard errors
    ## of the parameter estimators
    ## Note: manually; see QRM::fit.GEV and also copula:::fitCopula.ml():
    ##       fisher <- -hessian(logLik_GEV(fit$par, x = x)) # observed Fisher information (estimate of Fisher information)
    ##       Cov <- solve(fisher)
    ##       std.err <- sqrt(diag(Cov)) # see also ?fit_GEV_MLE
    if(estimate.cov) {
        negHessianInv <- catch(solve(-fit$hessian))
        if(is(negHessianInv, "error")) {
            warning("Hessian matrix not invertible: ", negHessianInv$error)
            Cov <- matrix(NA_real_, 0, 0)
            SE <- numeric(0)
        } else {
            Cov <- negHessianInv$value # result on warning or on 'worked'
            rownames(Cov) <- c("xi", "mu", "sigma")
            colnames(Cov) <- c("xi", "mu", "sigma")
            SE <- sqrt(diag(Cov))
            names(SE) <- c("xi", "mu", "sigma")
        }
    } else {
        Cov <- matrix(NA_real_, 0, 0)
        SE <- numeric(0)
    }
    ## Return (could create an object here)
    c(fit, list(Cov = Cov, SE = SE))
}
