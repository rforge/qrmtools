require(qrmtools)

## Setup
alpha <- 0.99
d <- seq(2, 1002, by=20)
theta <- c(0.1, 0.5, 1, 5, 10, 20)
palette <- colorRampPalette(c("darkorange2", "maroon3", "royalblue3", "black"), space="Lab")
cols <- palette(length(theta))
grid <- expand.grid(d=d, theta=theta)

## Compute best/worst VaR
bounds <- apply(grid, 1, function(x)
    VaR_bounds_hom(alpha, d=x[1], method="Wang.Par", theta=x[2]))
r <- cbind(grid, best=bounds[1,], worst=bounds[2,])
stopifnot(all(is.finite(r[,"best"])), all(is.finite(r[,"worst"])))
bVaR <- split(r[,"best"], r[,"theta"]) # best VaR
wVaR <- split(r[,"worst"], r[,"theta"]) # worst VaR

## Plot of best/worst VaR
plot(NA, xlim=range(d), ylim=range(bVaR, wVaR, na.rm=TRUE), log="xy",
     xlab=expression(italic(d)),
     ylab=as.expression(substitute("Best (dashed) and worst (solid) "~VaR[a]~
                                   "for Par("*theta*") margins", list(a=alpha))))
for(k in 1:length(theta)) {
    lines(d, wVaR[[k]], col=cols[k])
    lines(d, bVaR[[k]], col=cols[k], lty=2)
}
legend("topleft", bty="n", lty=rep(1,length(theta)), col=cols,
       legend=as.expression(lapply(1:length(theta),
       function(k) substitute(theta==k, list(k=theta[k])))))

## For debugging and bug-fixing purposes
if(FALSE) {
    ## Systematically search for numerical issues
    ii <- c()
    for(i in 1:nrow(grid)) {
        res  <- tryCatch(VaR_bounds_hom(alpha, d=grid[i,1], method="Wang.Par", theta=grid[i,2]), error=function(e) e)
        if(is(res, "simpleError")) {
            warning("Error for d=",grid[i,1],", theta=",grid[i,2],": ", conditionMessage(res))
            ii <- c(ii, i)
        }
    }

    ## Problem 1): d=1002, theta=10 (h2 not of opposite sign at endpoints)

    d <- 1002
    th <- 10

    ## Lower bound in x-scale (=> upper bound in c-scale)
    low <- if(th == 1) {
        d/2
    } else {
        (d-1)*(1+th)/(d-1+th)
    }
    ## Upper bound in x-scale (=> lower bound in c-scale)
    up <- if(th > 1) {
        (1+d/(th-1))^th
    } else if(th == 1) {
        e <- exp(1)
        (d+1)^(e/(e-1))
    } else { # th < 0
        d*th/(1-th)+1
    }
    (interval <- c(low, up))

    ## Internal auxiliary function to find root of
    h2 <- if(th == 1) {
        function(x) x^2 + x*(-d*log(x)+d-2)-(d-1)
    } else {
        function(x)
        (d/(1-th)-1)*x^(-1/th + 1) - (d-1)*x^(-1/th) + x - (d*th/(1-th) + 1)
    }

    ## Evaluate the auxiliary function h2 at the endpoints to see what fails
    h2(interval[1]) # should be negative
    h2(interval[2]) # should be positive => fails => adjust in this case


    ## Problem 2): d=500, theta=20 (worst VaR = Inf)

    d <- 500
    th <- 20

    ## Lower bound in x-scale (=> upper bound in c-scale)
    low <- if(th == 1) {
        d/2
    } else {
        (d-1)*(1+th)/(d-1+th)
    }
    ## Upper bound in x-scale (=> lower bound in c-scale)
    up <- if(th > 1) {
        (1+d/(th-1))^th
    } else if(th == 1) {
        e <- exp(1)
        (d+1)^(e/(e-1))
    } else { # th < 0
        d*th/(1-th)+1
    }
    (interval <- c(low, up))

    ## Internal auxiliary function to find root of
    h <- if(th == 1) {
        function(x) x^2 + x*(-d*log(x)+d-2)-(d-1)
    } else {
        function(x)
        (d/(1-th)-1)*x^(-1/th + 1) - (d-1)*x^(-1/th) + x - (d*th/(1-th) + 1)
        ## function(x)
        ## exp(x) - (d/(th-1)+1)*exp((1-1/th)*x) - (d-1)*exp(-x/th) +  (d*th/(th-1) - 1)
    }

    ## Evaluate the auxiliary function h2 at the endpoints to see what fails
    h(interval[1]) # should be negative
    h(interval[2]) # should be positive => fine

    ## Plot h2 there
    x <- seq(interval[1], interval[2], length.out=256)
    y <- h(x)
    plot(x, y, type="l")

    ## Root-finding on 'interval'
    tol <- .Machine$double.eps^0.25
    (x <- uniroot(h, interval=interval, tol=tol)$root) # => right end point of interval!
    ## => h(there) >> 0!
    x <- exp(x)
    c <- (1-alpha)/(x+d-1) # convert back to c-scale

    a <- alpha + (d-1)*c
    b <- 1-c # => problem: is 1
    qPar(a, theta=th)*(d-1) + qPar(b, theta=th)
    (d-1)*(1-a)^(-1/th)-d+c^(-1/th) # => works
}
