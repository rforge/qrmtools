### Tools for computing the worst VaR_alpha for given margins ##################

### 1) Crude VaR bounds (for both best and worst VaR) ##########################

##' @title Crude bounds for any VaR_alpha
##' @param alpha confidence level
##' @param d dimension
##' @param qF (list of) marginal quantile functions
##' @param ... ellipsis argument passed to qF()
##' @return 2-vector containing crude VaR_alpha bounds
##' @author Marius Hofert
crude_VaR_bounds <- function(alpha, d, qF, ...)
{
    if(is.function(qF))
        d * c(qF(alpha/d, ...), qF((d-1+alpha)/d, ...))
    else { # ... are passed to *all* qF()
        if(!is.list(qF))
            stop("qF has to be a (quantile) function or list of such")
        stopifnot(length(qF) == d)
        qF.low <- sapply(qF, function(qF.) qF.(alpha/d, ...))
        qF.up  <- sapply(qF, function(qF.) qF.((d-1+alpha)/d, ...))
        d * c(min(qF.low), max(qF.up))
    }
}


### 2) Explicit worst VaR in the homogeneous case ##############################

### Dual bound #################################################################

##' @title D(s,t) = d \int_{t}^{s-(d-1)t} \bar{F}(x) dx / (s-dt)
##' @param s real number
##' @param t real number < s/d
##' @param d dimension (integer > 2)
##' @param pF marginal distribution function (same for all d)
##' @param ... ellipsis argument passed to integrate()
##' @return D(s,t)
##' @author Marius Hofert
##' @note If t -> s/d-, l'Hospital's Rule shows that D(s, s/d) = d\bar{F}(s/d)
dual_bound_2 <- function(s, t, d, pF, ...)
{
    stopifnot(length(t) == 1)
    if(t > s/d) stop("t must be <= s/d")
    ## use D(s,t) = d( 1-\int_{t}^{s-(d-1)t} F(x) dx/(s-d*t) ) in this case
    if(t == s/d) d*(1-pF(s/d)) else
    d * (1 - (1/(s-d*t)) * integrate(pF, lower=t, upper=s-(d-1)*t, ...)$value)
}

##' @title Auxiliary function \bar{F}(t) + (d-1) * \bar{F}(s-(d-1)*t)
##' @param s real number
##' @param t real number < s/d
##' @param d dimension (integer > 2)
##' @param pF marginal distribution function (same for all d)
##' @return \bar{F}(t) + (d-1) * \bar{F}(s-(d-1)*t)
##' @author Marius Hofert
dual_bound_2_deriv_term <- function(s, t, d, pF)
    1-pF(t) + (d-1)*(1-pF(s-(d-1)*t))

##' @title Dual bound D(s)
##' @param s real number
##' @param d dimension (integer > 2)
##' @param pF marginal distribution function (same for all d)
##' @param ... ellipsis argument passed to dual_bound_2()'s integrate()
##' @return D(s)
##' @author Marius Hofert
##' @note The "first-order condition" (second equality in (14) in 2)) comes from the
##'       fact that
##'       (d/dt) D(s,t) = [ (-d)[\bar{F}(s-(d-1)t)(d-1)+\bar{F}(t)](s-dt) +
##'                         d^2 \int_{t}^{s-(d-1)t} \bar{F}(x) dx ] / (s-dt)^2 = 0
##'       if and only if
##'       d (\int_{t}^{s-(d-1)t} \bar{F}(x) dx) / (s-dt) = \bar{F}(s-(d-1)t)(d-1)-\bar{F}(t)
##'       => solving d (\int_{t}^{s-(d-1)t} \bar{F}(x) dx) / (s-dt) -
##'                  (\bar{F}(s-(d-1)t)(d-1)-\bar{F}(t)) = 0
##'          as a function in t for sufficiently large s leads to D(s)
dual_bound <- function(s, d, pF, tol=.Machine$double.eps^0.25, ...)
{
    stopifnot(length(s) == 1, s >= 0)
    if(s > 0) {
        ## h(s, t)
        h <- function(t) dual_bound_2(s, t=t, d=d, pF=pF, ...) -
            dual_bound_2_deriv_term(s, t=t, d=d, pF=pF)
        ## Note: h(t) -> 0 for h -> s/d- which is bad for uniroot() as the
        ##       latter will simply stop with the root t=s/d => we thus set f.upper > 0
        h.up <- -h(0) # guarantee that uniroot() doesn't fail due to root s/d
        t. <- uniroot(h, interval=c(0, s/d), f.upper=h.up, tol=tol)$root # optimal t in Equation (12) [= arginf]
        dual_bound_2_deriv_term(s, t=t., d=d, pF=pF) # dual bound D(s) in Equation (12) [= inf]
    } else {
        ## If s = 0, then t in [0, s/d] requires t to be 0 *and* f(0) = 0, so
        ## 0 is a root (as s/d). Furthermore, at t=0 (and with s=0),
        ## dual_bound_2_deriv_term(...) = d
        d
    }
}


### Wang's methods #############################################################

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


### Main wrapper function for computing the best/worst VaR in the homogeneous case

## Assumptions:
## - d=2: ultimately decreasing density (for x >= x0), alpha >= F(x0)
## - "Wang": F needs to live on [0, Inf), admitting a positive density which is
##           ultimately decreasing (for x >= x0), alpha >= F(x0)
## - "dual": F needs to be continuous with unbounded support and and ultimately
##           decreasing density, F(0) = 0 (otherwise, 0 as a lower bound for
##           uniroot() in dual_bound() is not valid)

##' @title Compute the best/worst VaR_\alpha in the homogeneous case with:
##'        1) d=2: Embrechts, Puccetti, Rueschendorf (2013, Proposition 2)
##'        2) d>=3:
##'           "Wang": Embrechts, Puccetti, Rueschendorf, Wang, Beleraj (2014, Prop. 3.1)
##'                   Integral evaluated numerically; needs smaller default
##'                   tolerance for uniroot()!
##'           "Wang.Par": The same, just with explicit formula for the integral
##'                       in the Pareto case; needs smaller default tolerance
##'                       for uniroot()!
##'           "Wang.Par.trafo": The same, just transforming the problem to
##'                             a different scale; for best VaR this is not needed
##'           "dual": Embrechts, Puccetti, Rueschendorf (2013, Proposition 4)
##'                   Numerically less stable; no formula for best VaR known (=> NA)
##' @param alpha confidence level
##' @param d dimension
##' @param method character string giving the method
##' @param interval initial interval
##' @param tol uniroot() x-tolerance
##' @param ... ellipsis arguments passed to Wang_h()
##' @return (best VaR, worst VaR) in the homogeneous case
##' @author Marius Hofert
VaR_bounds_hom <- function(alpha, d, method=c("Wang", "Wang.Par",
                          "Wang.Par.trafo", "dual"), interval=NULL, tol=NULL, ...)
{
    stopifnot(0<alpha, alpha<1, d>=2)
    method <- match.arg(method)

    ## Deal with d==2 first
    if(d==2) { # See Embrechts, Puccetti, Rueschendorf (2013, Prop. 2)
        qF <- NULL # make CRAN check happy
        if(!hasArg(qF))
            stop("The bivariate case requires the quantile function qF of F")
        qF <- list(...)$qF
        return(c(qF(alpha), 2*qF((1+alpha)/2)))
    }

    ## Best VaR for d >= 3 #####################################################

    best <-
        switch(method,
               "Wang" = {
                   qF <- NULL # make CRAN check happy
                   if(!hasArg(qF))
                       stop("Method 'Wang' requires the quantile function qF of F")
                   max((d-1)*qF(0)+qF(alpha), # Note: Typo in Wang, Peng, Yang (2013)
                       d*Wang_Ibar(a=0, b=alpha, alpha=alpha, d=d, ...))
               },
               "Wang.Par" =, "Wang.Par.trafo" = {
                   theta <- NULL # make CRAN check happy
                   if(!hasArg(theta))
                       stop("Method 'Wang.Par' and 'Wang.Par.trafo' require the parameter theta")
                   max((d-1)*qF(0)+qF(alpha), # Note: Typo in Wang, Peng, Yang (2013)
                       d*Wang_Ibar(a=0, b=alpha, alpha=alpha, d=d, method="Wang.Par", ...))
               },
               "dual" = { # "dual" only provides worst VaR
                   NA
               },
               stop("Wrong method"))

    ## Worst VaR for d >= 3  ###################################################

    if(is.null(tol)) # use smaller tolerance if required (matters; see also demo)
        tol <- if(method=="Wang" || method=="Wang.Par")
                   2.2204e-16 # MATLAB default
               else .Machine$double.eps^0.25 # uniroot() default
    worst <- switch(method,
           "Wang" = {

               ## Check qF()
               qF <- NULL # make CRAN check happy
               if(!hasArg(qF))
                   stop("Method 'Wang' requires the quantile function qF of F")
               ## Check 'interval'
               if(is.null(interval)) interval <- c(0, (1-alpha)/d)
               else {
                   if(interval[1] < 0) stop("interval[1] needs to be >= 0")
                   if(interval[1] > (1-alpha)/d) stop("interval[2] needs to be <= (1-alpha)/d")
                   if(interval[1] >= interval[2]) stop("interval[1] needs to be smaller than interval[2]")
               }

               ## Compute (and adjust) function values at endpoints
               h.low <- Wang_h(interval[1], alpha=alpha, d=d, ...)
               if(is.na(h.low))
                   stop("Objective function at interval[1] is NA or NaN. Provide a larger interval[1].")
               h.up <- -h.low # avoid that uniroot() fails due to 0 at upper interval endpoint

               ## Root-finding on 'interval'
               c. <- uniroot(function(c) Wang_h(c, alpha=alpha, d=d, ...),
                             interval=interval, f.lower=h.low, f.upper=h.up, tol=tol)$root
               d * Wang_h_aux(c., alpha=alpha, d=d, qF=list(...)$qF)

           },
           "Wang.Par" = {

               ## Check 'theta'
               theta <- NULL # make CRAN check happy
               if(!hasArg(theta))
                   stop("Method 'Wang.Par' requires the parameter theta")
               th <- list(...)$theta
               stopifnot(length(th) == 1, th > 0) # check theta here

               ## Compute lower uniroot boundaries and check
               if(is.null(interval)) {
                   low <- if(th > 1) {
                       (1-alpha)/((1+d/(th-1))^th+d-1)
                   } else if(th == 1) {
                       e <- exp(1)
                       (1-alpha)/((d+1)^(e/(e-1))+d-1)
                   } else { (1-th)*(1-alpha)/d }
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
           "Wang.Par.trafo" = { # here we compute the root on a different scale

               ## Check 'theta'
               theta <- NULL # make CRAN check happy
               if(!hasArg(theta))
                   stop("Method 'Wang.Par' requires the parameter theta")
               th <- list(...)$theta
               stopifnot(length(th) == 1, th > 0) # check theta here

               ## Compute uniroot boundaries (has to be in [1,Inf) here) and check
               if(is.null(interval)) {
                   low <- if(th == 1) d/2 else (d-1)*(1+th)/(d-1+th)
                   up <- if(th > 1) {
                       (1+d/(th-1))^th
                   } else if(th == 1) {
                       e <- exp(1)
                       (d+1)^(e/(e-1))
                   } else { d*th/(1-th)+1 }
                   interval <- c(low, up)
               } else {
                   if(interval[1] < 1) stop("interval[1] needs to be >= 1")
                   if(interval[1] >= interval[2]) stop("interval[1] needs to be smaller than interval[2]")
               }

               ## Define objective function (\tilde{\tilde{h}})
               h.tt <- if(th == 1) {
                   function(x) x^2 + x*(-d*log(x)+d-2)-(d-1)
               } else {
                   function(x)
                       (d/(1-th)-1)*x^(-1/th + 1) - (d-1)*x^(-1/th) + x - (d*th/(1-th) + 1)
               }

               ## Root-finding on 'interval'
               x. <- uniroot(h.tt, interval=interval, tol=tol)$root
               c. <- (1-alpha)/(x.+d-1) # convert back to c
               d * Wang_h_aux(c., alpha=alpha, d=d, method="Wang.Par", theta=th)

           },
           "dual" = {

               pF <- NULL # make CRAN check happy
               if(!hasArg(pF))
                   stop("Method 'dual' requires the distribution function pF")
               if(!hasArg(interval))
                   stop("Method 'dual' requires an initial interval c(s_l, s_u) to be given")
               uniroot(function(s) dual_bound(s, d=d, tol=tol, ...)-(1-alpha),
                       interval=interval, tol=tol)$root # s interval
               ## Note: We can't pass arguments to the inner root-finding

           },
           stop("Wrong method"))

           ## Return
           c(best, worst)
}


### 3) Worst VaR in the inhomogeneous case #####################################

##' @title Compute the number of columns oppositely ordered to the sum of all others
##' @param X (N, d)-matrix
##' @return Number of columns oppositely ordered to the sum of all others
##' @author Marius Hofert
##' @note Same time-saving tricks as in rearrange()
num_of_opp_ordered_cols <- function(X) {
    X.rs <- .rowSums(X, nrow(X), ncol(X)) # faster than rowSums()
    X.lst <- split(X, col(X)) # to avoid indexing the jth column, we work with a list!
    X.lst.sorted <- lapply(X.lst, sort.int, decreasing=TRUE) # sorting is only necessary once!
    sum(vapply(seq_len(ncol(X)),
               function(j) {
                   xj <- X.lst[[j]]
                   all(X.lst.sorted[[j]][rank(X.rs - xj)]
                       == xj)
               }, NA))
}

##' @title Basic rearrangement function for (A)RA
##' @param X (N, d)-matrix (either \underline{X}^\alpha or \overline{X}^\alpha)
##' @param tol Tolerance to determine (the individual) convergence;
##'        if NULL, the iteration is done until the matrix doesn't change
##' @param tol.type Character string indicating the tolerance function used
##'        ("relative" or "absolute")
##' @param maxiter Maximal number of iterations
##' @param method Character indicating which VaR is approximated (worst/best)
##'        determines optimizing function (min for worst VaR; max
##'        for best VaR)
##' @return List containing the
##'         1) Computed (lower or upper [depending on X]) bound for (worst or
##'            best [depending on method]) VaR
##'         2) (Individual) tolerance reached
##'         3) The (optimally) rearranged (N, d)-matrix X
##'         4) (N, .)-matrix of row sums (one column for each iteration)
##' @author Marius Hofert and Kurt Hornik
##' @note - We use "<= tol" to determine convergence instead of "< tol" as
##'         this then also nicely works with "= 0" (if tol=0) which stops in
##'         case the matrices are identical (no change at all).
##'       - No checking here due to speed
rearrange <- function(X, tol, tol.type=c("relative", "absolute"), maxiter=Inf,
                      method=c("worst", "best"))
{
    N <- nrow(X)
    d <- ncol(X)

    ## Define helper functions
    optim.fun <- if(method=="worst") min else max
    tol.fun <- if(tol.type=="absolute") {
        function(x, y) abs(x-y)
    } else {
        function(x, y) abs((x-y)/y)
    }

    ## Setup before major loop
    row.sums <- matrix(, nrow=nrow(X), ncol=0) # (N, 0)-matrix of row sums
    mrs.old <- optim.fun(.rowSums(X, N, d)) # old minimal row sum (currently: the one of X)

    ## Define lists
    X.rs <- .rowSums(X, N, d) # faster than rowSums()
    X.lst <- split(X, col(X)) # to avoid indexing the jth column all the time, we work with a list!
    X.lst.sorted <- lapply(X.lst, sort.int, decreasing=TRUE) # sorting is only necessary once!

    ## Loop through the columns
    while (TRUE) {

        ## Oppositely order X (=> Y)
        Y.lst <- X.lst
        Y.rs <- X.rs # row sum of Y
        unchanged <- TRUE # defined even if tol is not NULL
        for(j in 1:d) { # one iteration over all columns of the matrix
            yj <- Y.lst[[j]]
            rs <- Y.rs - yj # sum over all other columns (but the jth)
            yj <- X.lst.sorted[[j]][rank(rs)]
            if(unchanged && any(yj != Y.lst[[j]])) unchanged <- FALSE # check whether yj has changed
            Y.rs <- rs + yj # update row sum of Y
            Y.lst[[j]] <- yj # update list with rearranged jth column
        }

        ## Compute row sums and minimal/maximal row sums
        row.sums <- cbind(row.sums, Y.rs) # append the new row sums
        mrs.new <- optim.fun(Y.rs) # compute new minimal/maximal row sum

        ## Check convergence (we use "<= tol" as it entails tol=0)
        if((ncol(row.sums) == maxiter) ||
           (if(is.null(tol)) unchanged else tol.fun(mrs.new, mrs.old) <= tol))
        {
            tol <- tol.fun(mrs.new, mrs.old) # compute the (individual) tolerance
            break
        } else {
            mrs.old <- mrs.new # update mrs.old
            X.lst <- Y.lst
            X.rs <- Y.rs
        }

    }

    ## Return
    colnames(row.sums) <- NULL # remove column names so that they don't appear in output
    list(bound=mrs.new, # computed bound (\underline{s}_N or \overline{s}_N); required
         tol=tol, # tolerance for the computed bound; not costly
         X.rearranged=matrix(unlist(Y.lst), ncol=length(Y.lst)), # the rearranged X; computed anyways
         row.sums=row.sums) # the computed row sums after each iteration through all cols; computed anyways

}

##' @title Computing lower/upper bounds for the worst VaR with the RA
##' @param alpha confidence level
##' @param d dimension
##' @param qF a marginal quantile function (homogeneous case) or a d-list of such;
##'        note that each marginal quantile function has to be vectorized in p
##' @param N number of discretization points
##' @param abstol absolute convergence tolerance (to determine convergence);
##'        if NULL, then the iteration is done until the matrix doesn't change
##'        anymore, i.e., each column is oppositely ordered with respect to the
##'        sum of all others
##' @param maxiter maximal number of iterations
##' @param method character indicating which VaR is approximated (worst/best)
##' @param sample logical indicating whether each column of the two working
##'        matrices are sampled before the iteration begins
##' @return List containing the
##'         1) Computed lower and upper bound for (worst or best) VaR
##'         2) The relative rearrangement gap
##'            "|(upper bound - lower bound) / upper bound|"
##'         3) (Individual) tolerances reached (for each bound)
##'         4) List of (N, d) input matrices X (for each bound)
##'         5) List of rearranged Xs (for each bound)
##'         6) Number of iterations through the matrix columns used
##'         7) List of (N, .)-matrices of row sums (one column for each iteration;
##'            for each bound)
##'         8) Vectors of minimal [for worst VaR] or maximal [for best VaR] row sums
##'            (for each bound)
##' @author Marius Hofert
##' @note Notation is from p. 2757 in Embrechts, Puccetti, Rueschendorf (2013);
##'       variables are named according to the 'worst' VaR case.
RA <- function(alpha, d, qF, N, abstol=NULL, maxiter=Inf,
               method=c("worst", "best"), sample=TRUE)
{
    ## Checks and Step 1 (get N, abstol)
    stopifnot(0 < alpha, alpha < 1, is.null(abstol) || abstol >= 0,
              length(N) >= 1, N >= 2, maxiter >= 1, is.logical(sample))
    method <- match.arg(method)
    ## Checking of d, qF
    if(missing(d))
        stopifnot(is.list(qF), sapply(qF, is.function), (d <- length(qF)) >= 2)
    else { # extend the given quantile function to a list of functions
        stopifnot(d >= 2, is.function(qF))
        qF <- rep(list(qF), d)
    }

    ## Compute lower bound

    ## Step 2 (build \underline{X}^\alpha)
    p <- if(method=="worst") alpha + (1-alpha)*0:(N-1)/N else alpha*0:(N-1)/N # N-vector of prob.
    X <- sapply(qF, function(qF) qF(p))
    ## adjust those that are -Inf (for method="best")
    ## use alpha*((0+1)/2 / N) = alpha/(2N) instead of 0 quantile
    if(method == "best")
        X[1,] <- sapply(1:d, function(j) if(is.infinite(X[1,j])) qF[[j]](alpha/(2*N)) else X[1,j])

    ## Step 3 (randomly permute each column of \underline{X}^\alpha)
    if(sample) X <- apply(X, 2, sample)

    ## Steps 4, 5 (determine \underline{X}^*)
    ## repeat oppositely ordering \underline{X}^\alpha until there is only an
    ## abstol change in the min (method="worst") or max (method="best") row sum
    ## or until we reached maxiter number of iterations
    res.low <- rearrange(X, tol=abstol, tol.type="absolute", maxiter=maxiter,
                         method=method)

    ## Compute upper bound

    ## Step 2 (build \overline{X}^\alpha)
    p <- if(method=="worst") alpha + (1-alpha)*1:N/N else alpha*1:N/N # N-vector of prob.
    X. <- sapply(qF, function(qF) qF(p))
    ## adjust those that are Inf (for method="worst")
    ## use alpha+(1-alpha)*(N-1+N)/(2*N) = alpha+(1-alpha)*(1-1/(2*N)) instead of 1 quantile
    if(method == "worst") X.[N,] <- sapply(1:d, function(j)
        if(is.infinite(X.[N,j])) qF[[j]](alpha+(1-alpha)*(1-1/(2*N))) else X.[N,j])

    ## Step 3 (randomly permute each column of \overline{X}^\alpha)
    if(sample) X. <- apply(X., 2, sample)

    ## Step 6 (determine \overline{X}^*)
    ## repeat oppositely ordering \overline{X}^\alpha until there is only an
    ## abstol change in the min (method="worst") or max (method="best") row sum
    ## or until we reached maxiter number of iterations
    res.up <- rearrange(X., tol=abstol, tol.type="absolute", maxiter=maxiter,
                        method=method)

    ## Step 7 (return \underline{s}_N, \overline{s}_N and other non-time critical info)
    optim.fun <- if(method=="worst") min else max
    list(bounds=c(res.low$bound, res.up$bound), # (\underline{s}_N, \overline{s}_N)
         rel.ra.gap=abs((res.up$bound-res.low$bound)/res.up$bound), # relative RA gap
         individual.tol=c(res.low$tol, res.up$tol), # individual tolerances
         X=list(low=X, up=X.), # input matrices X (low, up)
         X.rearranged=list(low=res.low$X.rearranged, up=res.up$X.rearranged), # rearranged Xs (low, up)
         num.iter=c(ncol(res.low$row.sums), ncol(res.up$row.sums)), # number of iterations (low, up)
         row.sums=list(low=res.low$row.sums, up=res.up$row.sums), # row sums (low, up)
         m.row.sums=list(low=apply(res.low$row.sums, 2, optim.fun),
                         up=apply(res.up$row.sums, 2, optim.fun))) # optimal row sums (low, up)
}

##' @title Computing lower/upper bounds for the worst VaR with the ARA
##' @param alpha confidence level
##' @param d dimension
##' @param qF a marginal quantile function (homogeneous case) or a d-list of such;
##'        note that each marginal quantile function has to be vectorized in p
##' @param N.exp vector of exponents of 2 used as discretization points
##' @param reltol 2-vector of relative convergence tolerances
##'        for determining the individual relative tolerance (i.e., the relative
##'        tolerance in the minimal/maximal row sum for each of the bounds) and
##'        the joint relative tolerance (i.e., the relative
##'        tolerance between the computed lower and upper bounds). reltol[1] can be
##'        NULL (=> iteration is done until matrices don't change anymore, i.e.,
##'        until each column is oppositely ordered to the sum of all others,
##'        but only if maxiter hasn't been reached)
##' @param maxiter maximal number of iterations per N
##' @param method character indicating which VaR is approximated (worst/best)
##' @param sample logical indicating whether each column of the two working
##'        matrices are sampled before the iteration begins
##' @return List containing the
##'          1) Computed lower and upper bound for (worst or best) VaR
##'          2) The relative rearrangement gap
##'             "|(upper bound - lower bound) / upper bound|"
##'          3) (Individual) tolerances reached (for each bound)
##'          4) The joint tolerance reached between the two bounds
##'          5) List of (N, d) input matrices X (for each bound)
##'          6) List of rearranged Xs (for each bound)
##'          7) The number of discretization points used
##'          8) Number of iterations through the matrix columns used
##'          9) List of (N, .)-matrices of row sums (one column for each iteration;
##'             for each bound)
##'         10) Vectors of minimal [for worst VaR] or maximal [for best VaR] row sums
##'             (for each bound)
##' @author Marius Hofert
ARA <- function(alpha, d, qF, N.exp=seq(8, 20, by=1), reltol=c(0.001, 0.01),
                maxiter=12, method=c("worst", "best"), sample=TRUE)
{
    ## Checks and Step 1 (get N, reltol)
    stopifnot(0 < alpha, alpha < 1, length(reltol) == 2,
              is.null(reltol[1]) || reltol[1] >=0, reltol[2] >= 0,
              length(N.exp) >= 1, N.exp >= 1, maxiter >= 1, is.logical(sample))
    method <- match.arg(method)
    ## Checking of d, qF
    if(missing(d))
        stopifnot(is.list(qF), sapply(qF, is.function), (d <- length(qF)) >= 2)
    else { # extend the given quantile function to a list of functions
        stopifnot(d >= 2, is.function(qF))
        qF <- rep(list(qF), d)
    }

    ## Loop over N
    for(N. in 2^N.exp) {

        ## Compute lower bound

        ## Step 2 (build \underline{X}^\alpha)
        p <- if(method=="worst") alpha + (1-alpha)*0:(N.-1)/N. else alpha*0:(N.-1)/N. # N.-vector of prob.
        X <- sapply(qF, function(qF) qF(p))
        ## adjust those that are -Inf (for method="best")
        ## use alpha*((0+1)/2 / N.) = alpha/(2N.) instead of 0 quantile
        if(method == "best")
            X[1,] <- sapply(1:d, function(j) if(is.infinite(X[1,j])) qF[[j]](alpha/(2*N.)) else X[1,j])

        ## Step 3 (randomly permute each column of \underline{X}^\alpha)
        if(sample) X <- apply(X, 2, sample)

        ## Steps 4, 5 (determine \underline{X}^*)
        ## repeat oppositely ordering \underline{X}^\alpha until there is only an
        ## reltol[1] change in the min (method="worst") or max (method="best") row sum
        ## or until we reached maxiter number of iterations
        res.low <- rearrange(X, tol=reltol[1], tol.type="relative",
                             maxiter=maxiter, method=method)

        ## Compute upper bound

        ## Step 2 (build \overline{X}^\alpha)
        p <- if(method=="worst") alpha + (1-alpha)*1:N./N. else alpha*1:N./N. # N.-vector of prob.
        X. <- sapply(qF, function(qF) qF(p))
        ## adjust those that are Inf (for method="worst")
        ## use alpha+(1-alpha)*(N.-1+N.)/(2*N.) = alpha+(1-alpha)*(1-1/(2*N.)) instead of 1 quantile
        if(method == "worst") X.[N.,] <- sapply(1:d, function(j)
            if(is.infinite(X.[N.,j])) qF[[j]](alpha+(1-alpha)*(1-1/(2*N.))) else X.[N.,j])

        ## Step 3 (randomly permute each column of \overline{X}^\alpha)
        if(sample) X. <- apply(X., 2, sample)

        ## Step 6 (determine \overline{X}^*)
        ## repeat oppositely ordering \overline{X}^\alpha until there is only an
        ## reltol[1] change in the min (method="worst") or max (method="best") row sum
        ## or until we reached maxiter number of iterations
        res.up <- rearrange(X., tol=reltol[1], tol.type="relative",
                            maxiter=maxiter, method=method)

        ## Determine convergence (individual + joint)
        joint.tol <- abs((res.low$bound-res.up$bound)/res.up$bound)
        if( (max(res.low$tol, res.up$tol) <= reltol[1]) &&
            joint.tol <= reltol[2] ) break

    }

    ## Step 7 (return \underline{s}_N, \overline{s}_N and other non-time critical info)
    optim.fun <- if(method=="worst") min else max
    list(bounds=c(res.low$bound, res.up$bound), # (\underline{s}_N, \overline{s}_N)
         rel.ra.gap=abs((res.up$bound-res.low$bound)/res.up$bound), # relative RA gap
         individual.tol=c(res.low$tol, res.up$tol), # individual tolerances
         joint.tol=joint.tol, # joint tolerance
         X=list(low=X, up=X.), # input matrices X (low, up)
         X.rearranged=list(low=res.low$X.rearranged, up=res.up$X.rearranged), # rearranged Xs (low, up)
         N.used=N., # number of discretization points used
         num.iter=c(ncol(res.low$row.sums), ncol(res.up$row.sums)), # # of iterations (low, up) over all cols
         row.sums=list(low=res.low$row.sums, up=res.up$row.sums), # row sums (low, up) after each iteration over all cols
         m.row.sums=list(low=apply(res.low$row.sums, 2, optim.fun),
                         up=apply(res.up$row.sums, 2, optim.fun))) # optimal row sums (low, up)
}
