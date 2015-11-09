##' @title Conditional expectation (\bar{I}(a, b)) for computing best/worst VaR as
##'        in Embrechts, Puccetti, Rueschendorf, Wang, Beleraj (2014, Prop. 3.1)
##' @param a lower evaluation point
##' @param b upper evaluation point
##' @param alpha confidence level alpha
##' @param d dimension d
##' @param method character string giving the method
##'        generic = numerical integration; Wang.Par = Pareto distibution (explicit)
##' @param ... ellipsis argument passed to integrate()
##'        => must contain qF for "generic" and "theta" for "Wang.Par"
##' @return \bar{I}(a, b) = 1/(b-a)\int_a^b qF(y) dy =(subs) IE[L|L\in [qF(a), aF(b)]]
##' @author Marius Hofert
Wang_Ibar <- function(a, b, alpha, d, method=c("generic", "Wang.Par"), ...)
{
    stopifnot(length(a) == length(b), 0 <= a, a < b, b <= 1)
    ddd <- list(...)
    method <- match.arg(method)
    switch(method,
           "generic" = {
               stopifnot(length(a) == 1, length(b) == 1) # not vectorized (due to integrate())
               qF <- ddd$qF # grab out provided 'qF()'
               ddd$qF <- NULL # rm qF from '...'
               h <- function(...)
                   integrate(qF, lower=a, upper=b, ...)$value / (b-a)
               do.call(h, ddd) # call integrate() on the remaining arguments in '...'
           },
           "Wang.Par" = { # vectorized
               th <- ddd$theta # use provided 'theta'
               if(th == 1) log((1-a)/(1-b))/(b-a) - 1
               else (th/(1-th))*((1-b)^(1-1/th)-(1-a)^(1-1/th))/(b-a) - 1
           },
           stop("Wrong method"))
}

##' @title Right-hand side term in the objective function for computing the worst VaR
##'        as in Embrechts, Puccetti, Rueschendorf, Wang, Beleraj (2014, Prop. 3.1)
##' @param c evaluation point
##' @param alpha confidence level alpha
##' @param d dimension d
##' @param method character string giving the method
##'        generic = numerical integration; Wang.Par = Pareto distibution
##' @param ... ellipsis argument containing theta (for method="Wang.Par")
##'        or qF (for method="generic")
##' @return Right-hand side term in Prop. 3.1
##' @author Marius Hofert
##' @note for the correct 'c', this is the conditional expectation
Wang_h_aux <- function(c, alpha, d, method=c("generic", "Wang.Par"), ...)
{
    ddd <- list(...)
    method <- match.arg(method)
    qF <- if(method=="Wang.Par") function(y) qPar(y, theta=ddd$theta) else ddd$qF
    a <- alpha + (d-1)*c
    b <- 1-c
    qF(a)*(d-1)/d + qF(b)/d
}

##' @title Objective function for computing the worst VaR as in
##'        Embrechts, Puccetti, Rueschendorf, Wang, Beleraj (2014, Prop. 3.1)
##' @param c evaluation point
##' @param alpha confidence level alpha
##' @param d dimension d
##' @param method character string giving the method
##' @param ... ellipsis argument passed to Wang_h_aux() and Wang_Ibar()
##' @return objective function for computing the worst VaR
##' @author Marius Hofert
Wang_h <- function(c, alpha, d, method=c("generic", "Wang.Par"), ...)
{
    stopifnot(0 <= c, c <= (1-alpha)/d) # sanity check (otherwise b > a)
    method <- match.arg(method)
    ddd <- list(...)
    ## Properly deal with limit c=(1-alpha)/d
    Ib <- if(c == (1-alpha)/d) {
        if(method=="generic") { # qF() needs to be provided
            ddd$qF((d-1+alpha)/d)
        } else { # theta needs to be provided
            qPar((d-1+alpha)/d, theta=ddd$theta)
        }
    } else {
        Wang_Ibar(a=alpha+(d-1)*c, b=1-c, alpha=alpha, d=d, method=method, ...)
    }
    ## Return
    Ib - Wang_h_aux(c, alpha=alpha, d=d, method=method, ...)
}

## This has various flaws:
## - It needs a smaller default tolerance for uniroot() to work properly
## - For large d and both small and large thetas, the *theoretically* correct
##   uniroot() initial interval lower endpoint does not lead to opposite signs of h
##   at the interval endpoints (h is >0 as well at the lower endpoint)
##   => possible to fix for the theta > 1, but difficult to fix for theta < 1 (see below)
## => work on transformed scale!
"Wang.Par" = { # this requires a smaller default tolerance for uniroot()

    ## Check 'theta'
    theta <- NULL # make CRAN check happy
    if(!hasArg(theta))
        stop("Method 'Wang.Par' requires the parameter theta")
    th <- list(...)$theta
    stopifnot(length(th) == 1, th > 0) # check theta here

    ## Compute lower uniroot() boundaries and check
    if(is.null(interval)) {
        low <- if(th > 1) {
            ## We use an even smaller lower interval endpoint than the *theoretical*
            ## (1-alpha)/((1+d/(th-1))^th + d-1) to avoid that h is *numerically* > 0
            ## for large d (e.g., for d=542, theta=4, alpha=0.99).
            ## Note that (1-alpha)/((1+d/(th-1))^th + d-1) >= (1-alpha)/d^theta
            (1-alpha)/d^th
        } else if(th == 1) {
            e <- exp(1)
            (1-alpha)/((d+1)^(e/(e-1))+d-1)
        } else {
            ## Again we would like to use an even smaller lower interval endpoint than the *theoretical*
            ## (1-th)*(1-alpha)/d here to avoid that h is *numerically > 0 for large d and small theta
            ## (e.g., d=302, theta=0.1, alpha=0.99)
            ## => difficult to fix! as numbers are >>0, not just >0 => heavy numerical problem
            (1-th)*(1-alpha)/d
        }
        up <- if(th == 1) (1-alpha)/(3*d/2-1) else
                                                  (1-alpha)*(d-1+th)/((d-1)*(2*th+d))
        interval <- c(low, up)
    } else {
        if(interval[1] < 0) stop("interval[1] needs to be >= 0")
        if(interval[1] > (1-alpha)/d) stop("interval[2] needs to be <= (1-alpha)/d")
        if(interval[1] >= interval[2]) stop("interval[1] needs to be smaller than interval[2]")
    }
    if(th <= 1 && interval[1] == 0)
        stop("If theta <=1, interval[1] has to be > 0 as otherwise the internal Wang_h() is NaN")

    ## Root-finding on 'interval'
    c. <- uniroot(function(c)
        Wang_h(c, alpha=alpha, d=d, method="Wang.Par", ...),
        interval=interval, tol=tol)$root
    d * Wang_h_aux(c., alpha=alpha, d=d, method="Wang.Par", theta=th)

},
