### Fitting GPD(xi, beta) ######################################################

## ## Other packages
## - QRM:
##   + fit.GPD() has various options but, per default, uses probability-weighted
##     moments to get initial values and then MLE based on optim() *with*
##     provided gradient and observed information
## - Renext:
##   + fGPD(); uses a bit of a weird procedure based on the special cases
##     lomax and maxlo
##   + detailed in "Renext Computing Details" (but not available)
## - evd:
##   + fpot() -> fpot.norm()
##   + based on optim() with start values shape = 0, scale = 0
## - evir:
##   + ./R -> pot.R -> gpd(): based on optim() with initial values
##     s2 <- var(excess)
##     xi0 <- -0.5 * (((xbar * xbar)/s2) - 1)
##     beta0 <- 0.5 * xbar * (((xbar * xbar)/s2) + 1)
##     theta <- c(xi0, beta0)
##   + with hessian and observed information
## - fExtremes:
##   + gpdFit() -> .gpdmleFit(): as evir
## - ismev:
##   + gpd.fit(): uses optim(), hardcoded initial values
##     (shape = 0.1, scale = sqrt(6 * var(xdat))/pi)
## - texmex:
##   + evm() -> evm.default() -> evmFit(): based on optim(); uses
##     'family$start(data)' if start not provided (see also .pdf)
##   + 'start()' is a function and for the GPD found in ./R/gpd.R
##   + start() returns a longer vector (unclear why) but seems to use log(mean())
##     as initial value for xi and 1e-05 for beta; quite unclear

##' @title Computing Initial Values for MLE of the GPD
##' @param x numeric vector of data. In the POT method, these are the excesses
##'        over a sufficiently high threshold.
##' @param method character indicating the method to be used:
##'        - PWM: Matching probability-weighted moments; see Hosking and Wallis (1987).
##'               Findings of the latter:
##'               + MLE requires n ~>= 500 to be efficient
##'               + MoM reliable for xi > 0.2
##'               + PWM recommended for xi > 0.
##'               Note that their '-k' and 'alpha' is our 'xi' and 'beta'.
##'        - MoM: method-of-moments estimator
##'        - xi0: uses xi = 0 and method-of-moment estimators based on that case
##' @return numeric(2) specifying the initial values for xi, beta.
##' @author Marius Hofert
##' @note - No method guarantees a support of IR for the density and thus, for
##'         some data x, can lead to log-likelihood = -Inf, so optim() stops with error:
##'         "Error in optim(init, fn = function(param) logLik_GPD(param, x = x),
##'          hessian = estimate.cov: function cannot be evaluated at initial
##'          parameters"
##'       - Could be taken apart and exported (fit_GPD_PWM(), for example)
fit_GPD_init <- function(x, method = c("PWM", "MoM", "xi0"))
{
    method <- match.arg(method)
    init <- switch(method,
    "PWM" = {
        ## a_s = estimator (unbiased according to Landwehr, Matalas, Wallis (1979)) / a sample version of M_{1,0,s} = E(X (1-F(X))^s)
        x. <- sort(x)
        n <- length(x.)
        k <- 1:n
        a0 <- mean(x)
        a1 <- mean(x. * (n-k)/(n-1))
        a2 <- mean(x. * (n-k)*(n-k-1)/((n-1)*(n-2)))
        ## Estimators of xi and beta based on a0, a1, a2
        c(2 - a0/(a0 - 2 * a1) , (2 * a0 * a1) / (a0 - 2 * a1)) # xi, beta
    },
    "MoM" = {
        mu.hat <- mean(x)
        sig2.hat <- var(x)
        xi.hat <- (1-mu.hat^2/sig2.hat)/2
        c(xi.hat, mu.hat*(1-xi.hat)) # xi, beta
    },
    "xi0" = {
        ## Idea: use xi = 0
        c(0, mean(x)) # mean for xi = 0 is beta => beta
    },
    stop("Wrong 'method'"))
    ## Force initial values to produce finite log-likelihood if xi < 0
    ## (all x need to be in [0, -beta/xi) => choose beta = -xi * max(x) * 1.01)
    if(init[1] < 0) {
        mx <- max(x)
        if(mx >= -init[2]/init[1]) init[2] <- -init[1] * mx * 1.01 # beta = -xi * max(x) * 1.01
    }
    ## Return
    init
}

##' @title Log-likelihood of the GPD
##' @param param numeric(2) giving xi, beta (all real here; if beta <= 0
##'        dGPD(, log = TRUE) returns -Inf)
##' @param x numeric vector of data. In the POT method, these are the excesses
##'        over a sufficiently high threshold.
##' @return log-likelihood at xi, beta
##' @author Marius Hofert
logLik_GPD <- function(param, x)
    sum(dGPD(x, xi = param[1], beta = param[2], log = TRUE))

##' @title MLE for GPD Parameters
##' @param x numeric vector of data. In the POT method, these are the excesses
##'        over a sufficiently high threshold.
##' @param init numeric(2) giving the initial values for xi, beta; if NULL,
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
##' @note - Note that for no initial xi is the support of the density IR and thus
##'         the log-likelihood could be -Inf is (some) x are out of the support.
##'         This would lead to a failure of optim() at init (log-likelihood needs
##'         to be finite).
##'       - Most of the results generalize to the case where the GPD has another
##'         parameter for the location (like mu for the GEV distribution). If this
##'         mu would be chosen <= min(x), then one could probably guarantee a
##'         finite log-likelihood.
##'       - Similar to copula:::fitCopula.ml()
##'       - *Expected* information is available, too; see QRM::fit.GPD()
fit_GPD_MLE <- function(x, init = NULL, estimate.cov = TRUE, control = list(), ...)
{
    ## Checks
    stopifnot(is.numeric(x), is.null(init) || is.character(init) || (length(init) == 2),
              is.logical(estimate.cov), is.list(control))
    if(any(x < 0))
        stop("The support for the two-parameter GPD distribution is >= 0.")
    ## Compute initial values
    if(is.null(init)) {
        init <- fit_GPD_init(x)
    } else if(is.character(init)) { # easter egg
        init <- fit_GPD_init(x, method = match.arg(init, choices = eval(formals(fit_GPD_init)$method)))
    }
    ## Fit
    control <- c(as.list(control), fnscale = -1) # maximization (overwrites possible additionally passed 'fnscale')
    fit <- optim(init, fn = function(param) logLik_GPD(param, x = x),
                 hessian = estimate.cov, control = control, ...)
    ## Note: Could incorporate the gradient of the log-likelihood (w.r.t. xi, beta)
    ##       for the methods "BFGS", "CG" and "L-BFGS-B".
    ##       It is: (0, (-n+sum(x)/beta)/beta) for xi != 0 and
    ##       ( (1/xi^2)*sum(log1p(xi*x/beta))-(xi+1)*sum(xi*x/(beta+xi*x)),
    ##         -n/beta + (1+1/xi) * sum((xi*x/beta)/(beta+xi*x)) )
    ##
    ## Estimate of the asymptotic covariance matrix of the parameter estimators
    ## Note: manually:
    ##       fisher <- -hessian(logLik_GPD(fit$par, x = x)) # observed Fisher information (estimate of Fisher information)
    ##       Cov <- solve(fisher)
    ##       std.err <- sqrt(diag(Cov)) # see also ?fit_GPD_MLE
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
