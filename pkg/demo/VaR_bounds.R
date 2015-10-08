### Worst Value-at-Risk under known margins ####################################

require(qrmtools)
doPDF <- !dev.interactive(orNone=TRUE)


### 1) Homogeneous case (all margins equal) ####################################

## Setup
alpha <- 0.99 # confidence level
d <- 8 # dimension
qF <- function(p, th=2) qPar(p, theta=th) # Pareto quantile function
pF <- function(q, th=2) pPar(q, theta=th) # Pareto distribution function


### 1.1) Checks for method="dual" ##############################################

## Investigating h(s,t) = D(s,t) - ..., the function for the inner
## root-finding to compute D(s); see dual_bound()
s <- c(1, 5, 10, 100, 500, 1000)
t <- sapply(seq_along(s), function(i) {
    res <- exp(seq(log(1e-3), log(s[i]/d), length.out=257))
    res[length(res)] <- s[i]/d # to avoid numerical issues (t > s/d)
    res
})
f <- sapply(seq_along(s), function(i)
            sapply(t[,i], function(t.)
                   qrmtools:::dual_bound_2(s[i], t=t., d=d, pF=pF) -
                   qrmtools:::dual_bound_2_deriv_term(s[i], t=t., d=d, pF=pF)))
palette <- colorRampPalette(c("maroon3", "darkorange2", "royalblue3"), space="Lab")
cols <- palette(6)
if(doPDF)
    pdf(file=(file <- paste0("fig_worst_VaR_",alpha,"_hom_dual_h_Par=2_d=",d,".pdf")),
        width=6, height=6)
par(pty="s")
plot(t[,1], f[,1], type="l", log="x", xlim=range(t), ylim=range(f), col=cols[1],
     xlab="t", ylab=expression("h(s,t) for d = 8 and F being Par(2)"))
lines(t[,2], f[,2], col=cols[2])
lines(t[,3], f[,3], col=cols[3])
lines(t[,4], f[,4], col=cols[4])
lines(t[,5], f[,5], col=cols[5])
lines(t[,6], f[,6], col=cols[6])
abline(h=0, lty=2)
legend("topright", lty=rep(1,6), col=cols,
       bty="n", legend=as.expression(lapply(1:6,
           function(i) substitute(s==s., list(s.=s[i])))))
if(doPDF) dev.off.pdf(file=file)
## Conclusion: As we know, h(s, s/d) = 0. We also see that s has to be
##             sufficiently large in order to find a root h(s, t) = 0 for t < s/d

## Plot dual bound D(s) for various theta (concerns the outer root-finding)
theta <- c(0.5, 1, 2, 4) # theta values
s <- seq(48, 2000, length.out=257) # s values
D <- sapply(theta, function(th)
            sapply(s, function(s.)
                   dual_bound(s., d=8, pF=function(q) pPar(q, theta=th)))) # (s, theta) matrix
if(doPDF)
    pdf(file=(file <- paste0("fig_worst_VaR_",alpha,"_hom_dual_D_s_Par=",
                             paste0(theta, collapse="_"),"_d=",d,".pdf")),
        width=6, height=6)
par(pty="s")
plot(s, D[,1], type="l", ylim=range(D), col="maroon3",
     ylab=expression("Dual bound D(s) for d = 8 and F being Par("*theta*")"))
lines(s, D[,2], col="darkorange2")
lines(s, D[,3], col="royalblue3")
lines(s, D[,4], col="black")
legend("topright", lty=rep(1,4),
       col=c("maroon3", "darkorange2", "royalblue3", "black"),
       bty="n", legend=as.expression(lapply(1:4,
           function(j) substitute(theta==j, list(j=theta[j])))))
if(doPDF) dev.off.pdf(file=file)


### 1.2) Checks for method="Wang"/"Wang.Par" ###################################

### 1.2.1) Check of auxiliary functions *with* numerical integration (for theta=2)

## Check Wang_Ibar()
c <- seq(0, (1-alpha)/d-1e-8, length.out=129) # initial interval for root finding
Ib <- sapply(c, function(c.)
    qrmtools:::Wang_Ibar(a=alpha+(d-1)*c., b=1-c., alpha=alpha, d=d, qF=qF))
par(mar=c(5.1, 6.1, 4.1, 2.1)) # more space for the y-axis label
plot(c, Ib, type="l", xlab="c (in initial interval)",
     ylab=expression(bar(I)(c)))

## Check Wang_h_aux()
h.aux <- qrmtools:::Wang_h_aux(c=c, alpha=alpha, d=d, qF=qF)
plot(c, h.aux, type="l", xlab="c (in initial interval)",
     ylab=expression(frac(d-1,d)~{F^{-1}}(a[c])+frac(1,d)~{F^{-1}}(b[c])))

## Check objective function h(c) (Wang_h() with numerical integration) for
## the default theta
h <- sapply(c, function(c.) qrmtools:::Wang_h(c., alpha=alpha, d=d, qF=qF))
if(doPDF)
    pdf(file=(file <- paste0("fig_worst_VaR_",alpha,"_hom_Wang_h_Par=2_d=",d,"_num.pdf")),
        width=6, height=6)
par(pty="s")
plot(c, h, type="l", xlab="c (in initial interval)",
     ylab=expression("h(c) for"~~alpha~"= 0.99, d = 8 and F being Par(2)"))
abline(h=0, lty=2)
if(doPDF) dev.off.pdf(file=file)

## Check endpoints of objective function for root-finding
sapply(c(0, (1-alpha)/d), function(c.)
       qrmtools:::Wang_h(c., alpha=alpha, d=d, qF=qF)) # -Inf, 0
## -Inf is not a problem for root finding (for theta > 1; for theta <= 1 it is NaN),
## but the 0 at the right endpoint is => A proper initial interval [c_l,u_l]
## with 0<c_l<=c_u<(1-alpha)/d is required (containing the root).


### 1.2.2) Check of h(c) *without* numerical integration (for a range of thetas)

## Check objective function h(c) (Wang_h() without numerical integration)
d <- 8 # or d <- 100
c <- seq(0, (1-alpha)/d, length.out=2^13+1)
## => They all go further down to 0 if length.out is increased.
##    Smaller theta thus corresponds to a larger derivative in the root
##    Root-finding thus requires higher precision for smaller theta
h <- matrix(, nrow=length(c), ncol=length(theta))
for(j in 1:length(theta))
    h[,j] <- sapply(c, function(c.)
        qrmtools:::Wang_h(c., alpha=alpha, d=d, method="Wang.Par", theta=theta[j]))
z <- h
z[z <= 0] <- NA # > 0 => makes log-scale possible

## Plot
if(doPDF)
    pdf(file=(file <- paste0("fig_worst_VaR_",alpha,"_hom_Wang_h_Par_d=",d,".pdf")),
        width=6, height=6)
par(pty="s")
plot(c, z[,1], type="l", log="y", ylim=c(1e-4, 2e7), # works for d = 8 or 100; range(z, na.rm=TRUE)
     xlab="c", col="maroon3",
     ylab=substitute("h(c) for"~~alpha~"= 0.99, d ="~d.~"and F being Par("*theta*")",
                     list(d.=d)))
abline(h=0, lty=2)
lines(c, z[,2], col="darkorange2")
lines(c, z[,3], col="royalblue3")
lines(c, z[,4], col="black")
legend("topleft", lty=rep(1,4),
       col=c("maroon3", "darkorange2", "royalblue3", "black"),
       bty="n", legend=as.expression(lapply(1:4,
           function(j) substitute(theta==j, list(j=theta[j])))))
if(doPDF) dev.off.pdf(file=file)


### 1.2.3) Compute corresponding (via Wang.Par) worst VaR_alpha ################

## Note: - for theta in (0, 1], c_l has to be > 0
##       - we focus on worst VaR only
alpha. <- 1-2^seq(-0.001, -10, length.out=128) # alpha values (in (0,1); concentrated near 1)
d <- 8 # or d <- 100
w.VaR.Wang  <- sapply(theta, function(th)
    sapply(alpha., function(a) {
               VaR_bounds_hom(a, d=d,
                              method="Wang.Par", # or Wang.Par.trafo
                              theta=th)[2] # worst VaR
           })) # (alpha, theta) matrix

## Plot
if(doPDF)
    pdf(file=(file <- paste0("fig_worst_VaR_hom_Wang_Par_d=",d,".pdf")),
        width=6, height=6)
par(pty="s")
plot(1-alpha., w.VaR.Wang[,1], type="l", log="xy", col="maroon3",
     ylim=c(4, 6e10), # works for d = 8 or 100; range(w.VaR.Wang)
     xlab=expression(1-alpha),
     ylab=substitute(bar(VaR)[alpha]*group("(",L^{"+"},")")~
     "for d ="~d.~"and F being Par("*theta*")", list(d.=d)))
lines(1-alpha., w.VaR.Wang[,2], col="darkorange2")
lines(1-alpha., w.VaR.Wang[,3], col="royalblue3")
lines(1-alpha., w.VaR.Wang[,4], col="black")
legend("topright", lty=rep(1,4),
       col=c("maroon3", "darkorange2", "royalblue3", "black"),
       bty="n", legend=as.expression(lapply(1:4,
           function(j) substitute(theta==j, list(j=theta[j])))))
if(doPDF) dev.off.pdf(file=file)


### 2) Compare 'crude', "Wang" (numerical), "Wang.Par", "dual", 'RA' ###########

## Setup (as before)
alpha <- 0.99 # confidence level
d <- 8 # dimension


### 2.1) Why we need (at least) a lower endpoint for the root-finding procedure
###      if theta in (0,1] in the Pareto case (infinite mean)
interval <- c(0, (1-alpha)/d) # endpoints
method <- "Wang.Par" # note: this also holds for (the numerical) method="Wang"
th <- 0.99
qrmtools:::Wang_h(interval[1], alpha=alpha, d=d, method=method, theta=th) # NaN => uniroot() fails
## Note: Wang_h() is actually already NaN for c <= 1e-17
qrmtools:::Wang_Ibar(a=alpha+(d-1)*interval[1], b=1-interval[1], alpha=alpha,
                     d=d, method=method, theta=th) # Inf
qrmtools:::Wang_h_aux(interval[1], alpha=alpha, d=d, method=method, theta=th) # Inf


### 2.2) Graphical comparison ##################################################

## Setup
alpha <- 0.99 # confidence level
d <- 8 # or d <- 100
n.th <- 64 # number of thetas
th.low <- 0.2
th.up <- 5
th <- seq(th.low, th.up, length.out=n.th) # thetas
qFs <- lapply(th, function(th.) {th.; function(p) qPar(p, theta=th.)}) # n.th-vector of Pareto quantile functions
pFs <- lapply(th, function(th.) {th.; function(q) pPar(q, theta=th.)}) # n.th-vector of Pareto dfs
N <- 1e4 # number of discretization points for RA(); N=1e5 does not improve the situation

## Compute values (~30s (d=8; d=100: ~15min); all with default tol, see implementation)
res <- matrix(, nrow=n.th, ncol=7)
colnames(res) <- c("Wang", "Wang.Par", "Wang.Par.uniroot.tol",
                   "Wang.Par.trafo", "dual", "RA.low", "RA.up")
pb <- txtProgressBar(max=n.th, style=if(isatty(stdout())) 3 else 1) # setup progress bar
on.exit(close(pb)) # on exit, close progress bar
for(i in seq_len(n.th)) {
    ## "Wang" (numerical integration with smaller uniroot() tolerance; numerically critical)
    Wang.num.res <- tryCatch(VaR_bounds_hom(alpha, d=d, qF=qFs[[i]])[2], error=function(e) e)
    if(is(Wang.num.res, "simpleError")) { # warnings appeared
        warning("there was an error: ", conditionMessage(Wang.num.res), " (will use NA as result)")
        Wang.num.res <- NA
    }
    res[i,"Wang"] <- Wang.num.res
    ## "Wang.Par" (with smaller uniroot() tolerance)
    res[i,"Wang.Par"] <- VaR_bounds_hom(alpha, d=d, method="Wang.Par", theta=th[i])[2]
    ## "Wang.Par" (with uniroot()'s default tolerance)
    res[i,"Wang.Par.uniroot.tol"] <- VaR_bounds_hom(alpha, d=d, method="Wang.Par", theta=th[i],
                                                    tol=.Machine$double.eps^0.25)[2]
    ## "Wang.Par.trafo" (with uniroot()'s default tolerance)
    res[i,"Wang.Par.trafo"] <- VaR_bounds_hom(alpha, d=d, method="Wang.Par.trafo", theta=th[i])[2]
    ## "dual" (with uniroot()'s default tolerance)
    res[i,"dual"] <- VaR_bounds_hom(alpha, d=d, method="dual",
                                    interval=crude_VaR_bounds(alpha, d=d, qF=qFs[[i]]),
                                    pF=pFs[[i]])[2]
    ## Rearrangement Algorithm
    set.seed(271) # use the same random permutation for each theta
    RA. <- RA(alpha, d=d, qF=qFs[[i]], N=N, abstol=0.0001)
    res[i,"RA.low"] <- RA.$bounds[1]
    res[i,"RA.up"]  <- RA.$bounds[2]
    ## Progress
    setTxtProgressBar(pb, i) # update progress bar
}

## Plot
res. <- res/res[,"Wang.Par.trafo"] # standardize
if(doPDF) {
    file <- paste0("fig_worst_VaR_",alpha,"_hom_comparison_d=",d,"_N=",N,".pdf")
    pdf(file=(file), width=6, height=6)
}
par(pty="s")
plot(th, res.[,"Wang"], type="l", ylim=range(res., na.rm=TRUE),
     xlab=expression(theta), ylab=substitute(bar(VaR)[0.99]*group("(",L^{"+"},")")~
     "comparison (standardized) for d ="~d.~"and F being Par("*theta*")", list(d.=d)),
     col="gray80", lty=2, lwd=5) # Wang (with numerical integration)
lines(th, res.[,"Wang.Par"], col="gray50", lty=2, lwd=3) # Wang Pareto (wo num. integration)
lines(th, res.[,"Wang.Par.uniroot.tol"], col="maroon3", lty=1, lwd=1) # not as bad as without initial interval [c_l,c_u], still...
lines(th, res.[,"dual"], col="royalblue3", lty=2, lwd=2)
lines(th, res.[,"RA.low"], col="black", lty=3, lwd=1)
lines(th, res.[,"RA.up"], col="gray20", lty=2, lwd=1)
legend("topright", bty="n",
       col=c("gray80", "gray50", "maroon3", "royalblue3", "black", "gray20"),
       lty=c(2,2,1,2,3,2), lwd=c(5,3,1,2,1,1),
       legend=c("Wang (num. int.)", "Wang Pareto (wo num. int.)",
                "Wang Pareto (uniroot() tol.)", "Dual bound",
                "lower RA bound", "upper RA bound"))
if(doPDF) dev.off.pdf(file=file)
