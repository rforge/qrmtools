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
