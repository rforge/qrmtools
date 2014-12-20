### Worst Value-at-Risk under known margins ####################################

require(qrmtools)
doPDF <- !dev.interactive(orNone=TRUE)


### 1) Homogeneous case (all margins equal) ####################################

## setup
alpha <- 0.99
d <- 8
th <- 2
qF <- function(p) qPar(p, theta=th)
pF <- function(q) pPar(q, theta=th)


### 1.1) Checks for method="Wang" ##############################################

## check Wang_I()
c <- seq(0, (1-alpha)/d, length.out=129) # initial interval for root finding
yc <- sapply(c, function(c.) qrmtools:::Wang_I(c., alpha=alpha, d=d, qF=qF))
par(mar=c(5.1, 6.1, 4.1, 2.1)) # more space for the y-axis label
plot(c, yc, type="l", xlab="c  (initial interval)",
     ylab=expression(I(c)-(b[c]-a[c])~
     bgroup("(",frac(d-1,d)~{F^{-1}}(a[c])+frac(1,d)~{F^{-1}}(b[c]),")")~~
     "for"~F~"being Par(2)"))

## check Wang_h_aux()
yc. <- qrmtools:::Wang_h_aux(c=c, alpha=alpha, d=d, qF=qF)
plot(c, yc., type="l", xlab="c  (initial interval)",
     ylab=expression((b[c]-a[c])~
     bgroup("(",frac(d-1,d)~{F^{-1}}(a[c])+frac(1,d)~{F^{-1}}(b[c]),")")~~
     "for"~F~"being Par(2)"))

## check objective function Wang_h()
yc.. <- qrmtools:::Wang_h(c, alpha=alpha, d=d, qF=qF) # TODO: from here; why not 0 anymore?
if(doPDF) pdf(file=(file <- "fig_worst_VaR_0.99_hom_Wang_Par=2_d=8.pdf"),
              width=6.5, height=6.5)
par(pty="s")
plot(c, yc.., type="l", xlab="c  (initial interval)", ylab="h(c)")
abline(h=0, lty=2)
if(doPDF) dev.off.pdf(file=file)

## check endpoints of objective function for root-finding
qrmtools:::Wang_h(c(0, (1-alpha)/d), alpha=alpha, d=d, qF=qF) # -Inf, 0
## -Inf is not a problem for root finding, but the 0 at the right endpoint is.
## We take care of this by adjusting f.upper in the root finding

## plot worst VaR_alpha for various Par(theta)
theta <- c(0.5, 1, 2, 4) # theta values
alpha <- seq(0.1, 0.9, length.out=128) # alpha values (in (0,1))
worst.VaR.Wang  <- sapply(theta, function(th)
                          sapply(alpha, function(a)
                                 worst_VaR_hom(a, d=8,
                                               qF=function(p) qPar(p, theta=th)))) # (alpha, theta) matrix
TODO: fix
if(doPDF) pdf(file=(file <- "fig_worst_VaR_hom_Wang_Par_d=8.pdf"),
              width=6.5, height=6.5)
par(pty="s")
plot(alpha, worst.VaR.Wang[,1], type="l", ylim=range(worst.VaR.Wang[,1]),
     ylab=expression("Worst"~~VaR[alpha]))
lines(alpha, worst.VaR.Wang[,2], col="blue")
lines(alpha, worst.VaR.Wang[,3], col="orange")
lines(alpha, worst.VaR.Wang[,4], col="red")
legend("topright", inset=0.02, lty=rep(1,4), y.intersp=1.2,
       col=c("black", "blue", "orange", "red"),
       bty="n", legend=as.expression(lapply(1:4,
           function(i) substitute(theta==i, list(i=theta[i])))))
if(doPDF) dev.off.pdf(file=file)


### 1.2) Checks for method="dual" ##############################################

## plot dual bound D(s) for various theta
## this concerns the outer root-finding
theta <- c(0.5, 1, 2, 4) # theta values
s <- seq(48, 2000, length.out=257) # s values
D <- sapply(theta, function(th)
            sapply(s, function(s.)
                   dual_bound(s., d=8, pF=function(q) pPar(q, theta=th)))) # (s, theta) matrix
if(doPDF)
    pdf(file=(file <- "fig_worst_VaR_0.99_hom_dual_D_s_Par=0.5_1_2_4_d=8.pdf"),
        width=6, height=6)
par(pty="s")
plot(s, D[,1], type="l", ylim=range(D), ylab="Dual bound D(s)")
lines(s, D[,2], col="blue")
lines(s, D[,3], col="orange")
lines(s, D[,4], col="red")
legend("topright", inset=0.02, lty=rep(1,4), y.intersp=1.2,
       col=c("black", "blue", "orange", "red"),
       bty="n", legend=as.expression(lapply(1:4,
           function(i) substitute(theta==i, list(i=theta[i])))))
if(doPDF) dev.off.pdf(file=file)

## investigating h(s,t) = D(s,t) - ..., the function for the inner
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
palette <- colorRampPalette(c("red", "orange", "blue"), space="Lab")
cols <- palette(6)
if(doPDF)
    pdf(file=(file <- "fig_worst_VaR_0.99_hom_dual_D_st_Par=2_d=8.pdf"),
        width=6, height=6)
par(pty="s")
plot(t[,1], f[,1], type="l", log="x", xlim=range(t), ylim=range(f), col=cols[1],
     xlab="t", ylab=expression(h(s,t)))
lines(t[,2], f[,2], col=cols[2])
lines(t[,3], f[,3], col=cols[3])
lines(t[,4], f[,4], col=cols[4])
lines(t[,5], f[,5], col=cols[5])
lines(t[,6], f[,6], col=cols[6])
abline(h=0, lty=2)
legend("topright", inset=0.02, lty=rep(1,6), col=cols,
       bty="n", legend=as.expression(lapply(1:6,
           function(i) substitute(s==s., list(s.=s[i])))), y.intersp=1.2)
if(doPDF) dev.off.pdf(file=file)
## Conclusion: As we know, h(s, s/d) = 0. We also see that s has to be
##             sufficiently large in order to find a root h(s, t) = 0 for t < s/d


### 1.3) Compare various methods of VaR_bound_hom() ############################

## TODO AM: quite different values; also check with Wang's method but Pareto case
init <- crude_VaR_bounds(alpha, d=d, qF=qPar, theta=th)
(worst.VaR.Wang <- worst_VaR_hom(alpha, d=d, qF=function(p) qPar(p, theta=th)))
(worst.VaR.dual <- worst_VaR_hom(alpha, d=d, method="dual", interval=init, pF=pF))
(worst.VaR.Par  <- worst_VaR_hom(alpha, d=d, method="Pareto", theta=th))

