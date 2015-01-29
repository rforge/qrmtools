### Tools for computing the worst VaR_alpha for given margins ##################

### 1) Crude VaR bounds (for both best and worst VaR) ##########################

## Crude bounds for any VaR_alpha
crude_VaR_bounds <- function(alpha, d, qF, ...)
{
    if(is.function(qF))
        d * c(qF(alpha/d, ...), qF((d-1+alpha)/d, ...))
    else { # ... are passed to *all* qF()
        if(!is.list(qF))
            stop("qF has to be a (quantile) function or list of such")
        qF.low <- sapply(qF, function(qF.) qF.(alpha/d, ...))
        qF.up  <- sapply(qF, function(qF.) qF.((d-1+alpha)/d, ...))
        d * c(min(qF.low), max(qF.up))
    }
}


### 2) Explicit worst VaR in the homogeneous case ##############################

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
        ## h(s, t)
        h <- function(t) dual_bound_2(s, t=t, d=d, pF=pF, ...) -
            dual_bound_2_deriv_term(s, t=t, d=d, pF=pF)
        ## note: f(t) -> 0 for t -> s/d- this is bad for uniroot() which will
        ##       simply stop with the root t=s/d => we thus set f.upper > 0
        ##       alternatively, one could use (the more invasive)
        ##       interval = c(0, s/d - .Machine$double.eps^0.25)
        h.up <- sign(-h(0)) * .Machine$double.xmin # guarantee that uniroot() doesn't fail due to root s/d
        t. <- uniroot(h, interval=c(0, s/d), f.upper=h.up)$root # optimal t in Equation (12) [= arginf]
        dual_bound_2_deriv_term(s, t=t., d=d, pF=pF) # dual bound D(s) in Equation (12) [= inf]
    } else {
        ## if s = 0, then t in [0, s/d] requires t to be 0 *and* f(0) = 0, so
        ## 0 is a root (as s/d). Furthermore, at t=0 (and with s=0),
        ## dual_bound_2_deriv_term(...) = d
        d
    }
}


### Wang's methods #############################################################

## Conditional expectation (\bar{I}(c)) for computing the worst VaR according
## to Embrechts, Puccetti, Rueschendorf, Wang, Beleraj (2014, Prop. 3.1)
Wang_Ibar <- function(c, alpha, d, method=c("generic", "Wang.Par"), ...)
{
    a <- alpha + (d-1)*c
    b <- 1-c
    ddd <- list(...)
    method <- match.arg(method)
    switch(method,
           "generic" = {
               stopifnot(length(c) == 1) # not vectorized (due to integrate())
               qF <- ddd$qF # use provided 'qF()'
               ## properly deal with the limit c == (1-alpha)/d (empty integration interval)
               if(a == b) { # or c == (1-alpha)/d
                   qF((d-1+alpha)/d)
               } else { # now the non-empty case
                   ddd$qF <- NULL
                   h <- function(...)
                       integrate(qF, lower=a, upper=b, ...)$value / (b-a)
                   do.call(h, ddd) # call integrate() on the remaining arguments in '...'
               }
           },
           "Wang.Par" = {
               th <- ddd$theta # use provided 'theta'
               ## properly deal with the limit c == (1-alpha)/d (empty integration interval)
               ## (vectorized)
               res <- rep(NA, length(c))
               i <- a == b # or c == (1-alpha)/d
               if(any(i))   res[i] <- rep(qPar((d-1+alpha)/d, theta=th), sum(i))
               if(any(!i)) res[!i] <- if(th == 1) log((1-a[!i])/(1-b[!i]))/(b[!i]-a[!i]) - 1
                   else (th/(1-th))*((1-b[!i])^(1-1/th)-(1-a[!i])^(1-1/th))/(b[!i]-a[!i]) - 1
               res
           },
           stop("Wrong method"))
}

## Right-hand side term in the objective function for computing the worst VaR
## according to Embrechts, Puccetti, Rueschendorf, Wang, Beleraj (2014, Prop. 3.1)
## for the correct 'c', this is the conditional expectation
Wang_h_aux <- function(c, alpha, d, method=c("generic", "Wang.Par"), ...)
{
    ddd <- list(...)
    method <- match.arg(method)
    qF <- if(method=="Wang.Par") function(y) qPar(y, theta=ddd$theta) else ddd$qF
    a <- alpha + (d-1)*c
    b <- 1-c
    qF(a)*(d-1)/d + qF(b)/d
}

## Objective function for computing the worst VaR according to
## Embrechts, Puccetti, Rueschendorf, Wang, Beleraj (2014, Prop. 3.1)
Wang_h <- function(c, alpha, d, method=c("generic", "Wang.Par"), ...)
{
    stopifnot(0 <= c, c <= (1-alpha)/d) # sanity check (otherwise b > a)
    Wang_Ibar(c, alpha=alpha, d=d, method=method, ...) -
    Wang_h_aux(c, alpha=alpha, d=d, method=method, ...)
}


### Main wrapper function for computing the worst VaR in the homogeneous case ##

## Compute the worst VaR_\alpha in the homogeneous case with:
## 1) d=2: Embrechts, Puccetti, Rueschendorf (2013, Proposition 2)
## 2) d>=3:
##    "Wang": Embrechts, Puccetti, Rueschendorf, Wang, Beleraj (2014, Prop. 3.1)
##    "dual": Embrechts, Puccetti, Rueschendorf (2013, Proposition 4)
##
## Assumptions:
## - d=2: ultimately decreasing density (for x >= x0), alpha >= F(x0)
## - "Wang": F needs to live on [0, Inf), admitting a positive density which is
##           ultimately decreasing (for x >= x0), alpha >= F(x0)
## - "dual": F needs to be continuous with unbounded support and and ultimately
##           decreasing density, F(0) = 0 (otherwise, 0 as a lower bound for
##           uniroot() in dual_bound() is not valid)
worst_VaR_hom <- function(alpha, d, method=c("Wang", "Wang.Par", "dual"),
                          interval=NULL, tol=NULL, ...)
{
    stopifnot(0<alpha, alpha<1, d>=2)
    method <- match.arg(method)
    if(is.null(tol)) # use higher tolerance (MATLAB default) for Wang's approach (matters!)
        tol <- if(method!="dual") 2.2204e-16 else .Machine$double.eps^0.25

    ## d == 2
    if(d == 2) { # see Embrechts, Puccetti, Rueschendorf (2013, Proposition 2)
        qF <- NULL # make CRAN check happy
        if(!hasArg(qF))
            stop("The bivariate case requires the quantile function qF of F")
        return(2*list(...)$qF((1+alpha)/2))
    }

    ## d >= 3
    switch(method,
           "Wang" = {
               ## check qF()
               qF <- NULL # make CRAN check happy
               if(!hasArg(qF))
                   stop("Method 'Wang' requires the quantile function qF of F")
               ## check 'interval'
               if(is.null(interval)) interval <- c(0, (1-alpha)/d)
               else if(!(0 <= interval[1] && interval[2] <= (1-alpha)/d &&
                         interval[1]<interval[2]))
                   stop("interval[1] needs to be >= 0, interval[2] needs to be <= (1-alpha)/d")
               ## c_u: guarantee that uniroot() doesn't fail due to root at (1-alpha)/d
               h.up <- .Machine$double.xmin
               ## c_l: compute and check whether we have opposite sign
               ##      idea: give good error message
               h.low <- Wang_h(interval[1], alpha=alpha, d=d, ...)
               if(is.na(h.low))
                   stop("Objective function at interval[1] is not a number. Provide a larger interval[1].")
               if(h.up * h.low >= 0)
                   stop("Objective function at end points of 'interval' not of opposite sign. Provide a smaller interval[1].")
               ## root-finding on 'interval'
               c. <- uniroot(function(c) Wang_h(c, alpha=alpha, d=d, ...),
                             interval=interval, f.lower=h.low, f.upper=h.up, tol=tol)$root
               d * Wang_h_aux(c., alpha=alpha, d=d, qF=list(...)$qF)
           },
           "Wang.Par" = { # only kept for ruling out problems due to numerical integration
               ## check 'theta'
               theta <- NULL # make CRAN check happy
               if(!hasArg(theta))
                   stop("Method 'Wang.Par' requires the parameter theta")
               th <- list(...)$theta
               stopifnot(length(th) == 1, th > 0) # check theta here
               ## check 'interval'
               if(is.null(interval)) interval <- c(0, (1-alpha)/d)
               else if(!(0 <= interval[1] && interval[2] <= (1-alpha)/d &&
                         interval[1]<interval[2]))
                   stop("interval[1] needs to be >= 0, interval[2] needs to be <= (1-alpha)/d")
               ## c_u: guarantee that uniroot() doesn't fail due to root at (1-alpha)/d
               h.up <- .Machine$double.xmin
               ## c_l: compute and check whether we have opposite sign
               ##      idea: give good error message
               h.low <- Wang_h(interval[1], alpha=alpha, d=d, method="Wang.Par", ...)
               if(is.na(h.low))
                   stop("Objective function at interval[1] is not a number. Provide a larger interval[1].")
               if(h.up * h.low >= 0)
                   stop("Objective function at end points of 'interval' not of opposite sign. Provide a smaller interval[1].")
               ## root-finding on 'interval'
               c. <- uniroot(function(c)
                             Wang_h(c, alpha=alpha, d=d, method="Wang.Par", ...),
                             interval=interval, f.lower=h.low, f.upper=h.up, tol=tol)$root
               d * Wang_h_aux(c., alpha=alpha, d=d, method="Wang.Par", theta=th)
           },
           "dual" = {
               pF <- NULL # make CRAN check happy
               if(!hasArg(pF))
                   stop("Method 'dual' requires the distribution function pF")
               if(!hasArg(interval))
                   stop("Method 'dual' requires an initial interval")
               uniroot(function(s) dual_bound(s, d=d, ...)-(1-alpha),
                       interval=interval, tol=tol)$root # s interval
               ## note: we can't pass arguments to the inner root-finding (yet)
           },
           stop("Wrong method"))
}


### 3) Worst VaR in the inhomogeneous case #####################################

### Rearrangement algorithm ####################################################

##' @title Iteratively oppositely order a column of a matrix with respect to
##'        the sum of all other colums
##' @param x matrix
##' @return a matrix as x after one iteration over all columns
##' @author Marius Hofert
##' @note 'Iteratively' here means that the jth column already involves the
##'       changed values from the (j-1)th column
oppositely_order <- function(x)
{
    ## stopifnot(is.matrix(x), (d <- ncol(x)) >= 2); no checking here due to speed
    for(j in seq_len(d)) x[,j] <- sort(x[,j], decreasing=TRUE)[rank(rowSums(x[,-j]))]
    x
}

## TODO: - comment here + export + document all + use below as demo
quantile_matrix <- function(alpha, qmargins, N, method=c("worst", "best")) {
    stopifnot(0 <= alpha, alpha <= 1, is.list(qmargins),
              (d <- length(qmargins)) >= 2, N >= 2)
    method <- match.arg(method)
    ## compute quantiles
    p <- if(method=="worst") alpha + (1-alpha)*0:N/N else alpha*0:N/N
    res <- sapply(qmargins, function(qF) qF(p))
    ## adjust those that are +/- Inf
    ## use alpha+(1-alpha)*(N-1+N)/(2*N) = alpha+(1-alpha)*(1-1/(2*N)) instead of 1 quantile
    if(method == "worst")
        res[N+1,] <- sapply(1:d, function(j) if(is.infinite(res[N+1,j]))
                            qmargins[[j]](alpha+(1-alpha)*(1-1/(2*N))) else res[N+1,j])
    else # method="best"; use alpha*(0+1)/2 = alpha/2 instead of 0 quantile
        res[1,] <- sapply(1:d, function(j) if(is.infinite(res[1,j]))
                            qmargins[[j]](alpha/2) else res[1,j])
    res
}
quantile_matrix(0.9, qmargins=list(qF1=qnorm, qF2=function(p) qPar(p, theta=2)),
                N=10)


##' @title Computing lower/upper Bounds for the Worst VaR
##' @param x (N, d)-matrix of quantiles containing
##'        F_j^{-1}(alpha + (1-alpha)*i/N), i=0,..,N (if method="worst");
##'        F_j^{-1}(alpha*i/N), i=0,..,N (if method="best")
##' @param alpha confidence level
##' @param eps epsilon error to determine convergence
##' @param method character indicating which VaR is approximated (worst/best)
##' @param verbose logical indicating whether additional output is written
##' @return lower and upper bounds for the worst (or best) VaR
##' @author Marius Hofert
##' @note - Notation is from p. 2757 in Embrechts, Puccetti, Rueschendorf (2013);
##'         variables are named according to the 'worst' VaR case.
##'       - We use "<= eps" to determine convergence instead of "< eps" as
##'         this then also nicely works with "= 0" (if eps=0) which stops in
##'         case the matrices are identical (no change at all).
RA <- function(x, alpha, eps=0, method=c("worst", "best"), verbose=FALSE)
{
    ## Checks and Step 1 (get N, eps)
    stopifnot(is.matrix(x), (nr <- nrow(x)) >= 3, (d <- ncol(x)) >= 2,
              0 < alpha, alpha < 1, eps >= 0)
    if(!all(is.finite(x))) stop("x must contain only finite numerical values")
    method <- match.arg(method)
    N <- nr-1 # number of discretization points N; N >= 2

    ## Steps 2, 3 (build/get \underline{X}^\alpha, \overline{X}^\alpha and permute their columns)
    X.low <- apply(x[1:(N-1),], 2, sample) # = \underline{X}^\alpha; for lower bound
    X.up  <- apply(x[2:N,],     2, sample) # = \overline{X}^\alpha; for upper bound

    ## Steps 4, 5 (determine \underline{X}^*)
    ## repeat oppositely ordering \underline{X}^\alpha until there is only an
    ## eps change in the min (method=="worst") or max (method=="best") row sum
    fun <- if(method=="worst") min else max
    while (TRUE) {
        Y <- oppositely_order(X.low)
        mrs <- c(fun(rowSums(Y)), fun(rowSums(X.low))) # minimal/maximal row sum
        if(verbose)
            if(method=="worst")
                cat("Min of row sums of Y: ", mrs[1],"\n")
            else
                cat("Max of row sums of Y: ", mrs[1],"\n")
        if(abs(mrs[1]-mrs[2]) <= eps) break # we use equality here as it entails eps=0
        else X.low <- Y
    }
    X.low.star <- Y # take the last rearranged matrix Y as \underline{X}^*
    s.low <- mrs[1] # with corresponding minimal/maximal row sum

    ## Step 6 (determine \overline{X}^*)
    ## repeat oppositely ordering \overline{X}^\alpha until there is only an
    ## eps change in the min (method=="worst") or max (method=="best") row sum
    while (TRUE) {
        Y <- oppositely_order(X.up)
        mrs <- c(fun(rowSums(Y)), fun(rowSums(X.up))) # minimal/maximal row sum
        if(verbose)
            if(method=="worst")
                cat("Min of row sums of Y: ", mrs[1],"\n")
            else
                cat("Max of row sums of Y: ", mrs[1],"\n")
        if(abs(mrs[1]-mrs[2]) <= eps) break # we use equality here as it entails eps=0
        else X.up <- Y
    }
    X.up.star <- Y # take the last rearranged matrix Y as \overline{X}^*
    s.up <- mrs[1] # with corresponding minimal/maximal row sum

    ## Step 7 (return \underline{s}_N, \overline{s}_N)
    c(lower=s.low, upper=s.up) # return the range; in practice often very narrow and ~= worst/best VaR
}


