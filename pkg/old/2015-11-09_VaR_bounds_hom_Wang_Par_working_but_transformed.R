           "Wang.Par" = {

               ## Here we compute the root on a different scale as numerically more
               ## challenging otherwise

               ## Check 'theta'
               theta <- NULL # make CRAN check happy
               if(!hasArg(theta))
                   stop("Method 'Wang.Par' requires the parameter theta")
               th <- list(...)$theta
               stopifnot(length(th) == 1, th > 0) # check theta here

               ## Compute uniroot() initial interval and check it
               ## Note: c-bound = (1-alpha)/(x-bound + d-1) or
               ##       x-bound = (1-alpha)/c-bound - (d-1)
               ##       (and lower bound changes to upper bound and vice versa)
               if(is.null(interval)) {
                   ## Lower bound in x-scale (= upper bound in c-scale)
                   low <- if(th == 1) {
                       d/2
                   } else {
                       (d-1)*(1+th)/(d-1+th)
                   }
                   ## Upper bound in x-scale (= lower bound in c-scale)
                   up <- if(th > 1) {
                       ## Theoretical bound: (1+d/(th-1))^th
                       ## => We use an even larger numerical bound here to avoid that
                       ##    uniroot() fails due to h not being of opposite sign at
                       ##    interval endpoints.
                       2 * (1+d/(th-1))^th # instead of (1+d/(th-1))^th
                   } else if(th == 1) {
                       e <- exp(1)
                       (d+1)^(e/(e-1))
                   } else { # th < 0
                       d*th/(1-th)+1
                   }
                   interval <- c(low, up)
               } else {
                   ## The user provides 'interval' in 'c-scale' (c in [0, (1-alpha)/d])
                   if(interval[1] < 0) stop("interval[1] needs to be >= 0")
                   if(interval[1] > (1-alpha)/d) stop("interval[2] needs to be <= (1-alpha)/d")
                   if(interval[1] >= interval[2]) stop("interval[1] needs to be smaller than interval[2]")
                   ## We convert 'interval' to the 'x-scale' (x in [1,Inf)) to determine the root
                   interval <- (1-alpha)/rev(interval) - (d-1)
               }

               ## Define objective function (transformed 'h')
               h. <- if(th == 1) {
                   function(x) x^2 + x*(-d*log(x)+d-2)-(d-1)
               } else {
                   function(x)
                       (d/(1-th)-1)*x^(-1/th + 1) - (d-1)*x^(-1/th) + x - (d*th/(1-th) + 1)
               }

               ## Root-finding on 'interval'
               x <- uniroot(h., interval=interval, tol=tol)$root
               c <- (1-alpha)/(x+d-1) # convert back to c-scale
               d * Wang_h_aux(c, alpha=alpha, d=d, method="Wang.Par", theta=th)

           },
