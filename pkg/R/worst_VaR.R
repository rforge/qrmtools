### Tools for computing the worst VaR_alpha for given margins ##################

### Crude VaR bounds (for both best and worst VaR) #############################

## Crude VaR bounds in the homogeneous case with marginal quantile function qF
## The lower bound is obtained by the Poincare--Sylvester sieve formula
## The upper bound is obtained by the lower Frechet--Hoeffding bound
VaR_hom_crude <- function(alpha, d, qF, ...)
    d * c(qF(alpha, ...), qF(1-(1-alpha)/d, ...))


### Dual bound #################################################################

## D(s,t) = d \int_{t}^{s-(d-1)t} \bar{F}(x) dx / (s-dt)
## s is any real number, d > 2, t < s/d. If t -> s/d-, l'Hospital's Rule shows
## that D(s, s/d) = d\bar{F}(s/d)
dual_bound_2 <- function(s, t, d, pF, ...)
{
    stopifnot(length(t) == 1)
    if(t > s/d) stop("t must be <= s/d")
    ## use D(s,t) = d( 1-\int_{t}^{s-(d-1)t} F(x) dx/(s-d*t) ) in this case
    if(t == s/d) d*(1-pF(s/d)) else
    d * (1 - (1/(s-d*t)) * integrate(pF, lower=t, upper=s-(d-1)*t, ...)$value)
}

## Return function \bar{F}(t) + (d-1) * \bar{F}(s-(d-1)*t)
dual_bound_2_deriv_term <- function(s, t, d, pF)
    1-pF(t) + (d-1)*(1-pF(s-(d-1)*t))

## Dual bound D(s)
## note: - f will be optimized in t
##       - f depends on s (not suitable for vectorization)
## The "first-order condition" (second equality in (14) in 2)) comes from the
## fact that
## (d/dt) D(s,t) = [ (-d)[\bar{F}(s-(d-1)t)(d-1)+\bar{F}(t)](s-dt) +
##                 d^2 \int_{t}^{s-(d-1)t} \bar{F}(x) dx ] / (s-dt)^2 = 0
## if and only if
## d (\int_{t}^{s-(d-1)t} \bar{F}(x) dx) / (s-dt) = \bar{F}(s-(d-1)t)(d-1)-\bar{F}(t)
## => solving d (\int_{t}^{s-(d-1)t} \bar{F}(x) dx) / (s-dt) -
##          (\bar{F}(s-(d-1)t)(d-1)-\bar{F}(t)) = 0
##    as a function in t for sufficiently large s leads to D(s)
dual_bound <- function(s, d, pF, ...)
{
    stopifnot(length(s) == 1, s >= 0)
    if(s > 0) {
        f <- function(t) dual_bound_2(s, t=t, d=d, pF=pF, ...) -
            dual_bound_2_deriv_term(s, t=t, d=d, pF=pF)
        ## note: f(t) -> 0 for t -> s/d- this is bad for uniroot() which will
        ##       simply stop with the root t=s/d => we thus set f.upper > 0
        ##       alternatively, one could use (the more invasive)
        ##       interval = c(0, s/d - .Machine$double.eps^0.25)
        f1 <- sign(-f(0)) * .Machine$double.xmin # guarantee that uniroot() doesn't fail due to root s/d
        t. <- uniroot(f, interval=c(0, s/d), f.upper=f1)$root # optimal t in Equation (12) [= arginf]
        dual_bound_2_deriv_term(s, t=t., d=d, pF=pF) # dual bound D(s) in Equation (12) [= inf]
    } else {
        ## if s = 0, then t in [0, s/d] requires t to be 0 *and* f(0) = 0, so
        ## 0 is a root (as s/d). Furthermore, at t=0 (and with s=0),
        ## dual_bound_2_deriv_term(...) = d
        d
    }
}


### Wang's methods #############################################################

## Integral qF() for computing the worst VaR according to
## Embrechts, Puccetti, Rueschendorf, Wang, Beleraj (2014, Prop. 3.1)
Wang_int <- function(c, alpha, d, ...)
{
    ## extract qF from '...'
    ddd <- list(...)
    qF <- ddd$qF
    ddd$qF <- NULL
    f <- function(...) # sapply() here to vectorize
        sapply(c, function(c.) integrate(qF, lower=alpha+(d-1)*c.,
                                         upper=1-c., ...)$value)
    do.call(f, ddd) # call integrate() on the remaining arguments in ...
}

## Right-hand side term for computing the worst VaR according to
## Embrechts, Puccetti, Rueschendorf, Wang, Beleraj (2014, Prop. 3.1)
Wang_obj_aux <- function(c, alpha, d, qF)
    ((1-alpha-d*c)/d) * ((d-1)*qF(alpha+(d-1)*c) + qF(1-c))

## Objective function for computing the worst VaR according to
## Embrechts, Puccetti, Rueschendorf, Wang, Beleraj (2014, Prop. 3.1)
Wang_obj <- function(c, alpha, d, ...)
    Wang_int(c, alpha=alpha, d=d, ...) -
    Wang_obj_aux(c, alpha=alpha, d=d, qF=list(...)$qF)


### Main function for computing the worst VaR in the homogeneous case ##########

## Compute the worst VaR_\alpha in the homogeneous case with
## 1) Embrechts, Puccetti, Rueschendorf, Wang, Beleraj (2014, Prop. 3.1) [best]
## 2) Embrechts, Puccetti, Rueschendorf (2013, Proposition 4) [via dual bound]
## The Pareto case uses Wang's method but specifically for the Pareto distribution
##
## Non-optimal things for 2):
## we need a "sufficiently large s" (otherwise the inner root finding algorithm
## fails as there is no root below (1-alpha)/d) => good initial interval needed
##
## Assumptions:
## - for 1): F needs to be continuous with ultimately decreasing density
## - for 2): F needs to be continuous with unbounded support and and ultimately
##   decreasing density, F(0) = 0 (otherwise, 0 as a lower bound for uniroot()
##   in dual_bound() is not valid)
worst_VaR_hom <- function(alpha, d, interval=NULL,
                          method=c("Wang", "dual", "Pareto"), ...)
{
    stopifnot(0 < alpha, alpha < 1, d >= 3)
    method <- match.arg(method)
    switch(method,
           "Wang" = {
               if(!hasArg(qF)) # just check that 'qF' has been provided
                   stop("Method 'Wang' requires the quantile function qF of F")
               if(is.null(interval)) interval <- c(0, (1-alpha)/d)
               f1 <- .Machine$double.xmin # guarantees that uniroot() doesn't fail due to root at (1-alpha)/d
               c. <- uniroot(function(c) Wang_obj(c, alpha=alpha, d=d, ...),
                             interval=interval, f.upper=f1)$root
               (d/(1-alpha-d*c.)) * Wang_obj_aux(c., alpha=alpha, d=d,
                                                 qF=list(...)$qF)
           },
           "dual" = {
               pF <- NULL # make CRAN check happy
               if(!hasArg(pF))
                   stop("Method 'dual' requires the distribution function pF")
               if(!hasArg(interval))
                   stop("Method 'dual' requires an initial interval")
               uniroot(function(s) dual_bound(s, d=d, ...)-(1-alpha),
                       interval=interval)$root # s interval
               ## note: we can't pass arguments to the inner root-finding (yet)
           },
           "Pareto" = {
               theta <- NULL # make CRAN check happy
               if(!hasArg(theta))
                   stop("Method 'Pareto' requires the parameter theta")
               qF <- function(p) qPar(p, theta=list(...)$theta)
               stop("TODO AM")
           },
           stop("Wrong method"))
}

