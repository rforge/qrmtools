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

## check Wang_int()
c <- seq(0, (1-alpha)/d, length.out=128) # initial interval for root finding
yc <- qrmtools:::Wang_int(c, alpha=alpha, d=d, qF=qF)
par(mar=c(5.1, 6.1, 4.1, 2.1)) # more space for the y-axis label
plot(c, yc, type="l", xlab="c  (initial interval)",
     ylab=expression(integral({F^{-1}}(y)~dy, alpha+(d-1)*c, 1-c)~~
     "for"~F~"being Par(2)"))

## check Wang_obj_aux()
yc. <- qrmtools:::Wang_obj_aux(c=c, alpha=alpha, d=d, qF=qF)
plot(c, yc., type="l", xlab="c  (initial interval)",
     ylab=expression(frac(1-alpha-d*c,d)~((d-1)*{F^{-1}}(alpha+(d-1)*c)+{F^{-1}}(1-c))~~
     "for"~F~"being Par(2)"))

## check Wang_obj()
yc.. <- qrmtools:::Wang_obj(c, alpha=alpha, d=d, qF=qF)
if(doPDF) pdf(file=(file <- "fig_worst_VaR_0.99_hom_Wang_Par=2_d=8.pdf"),
              width=6.5, height=6.5)
par(pty="s")
plot(c, yc.., type="l", xlab="c  (initial interval)",
     ylab=expression(integral({F^{-1}}(y)~dy, alpha+(d-1)*c, 1-c)-
     frac(1-alpha-d*c,d)~((d-1)*{F^{-1}}(alpha+(d-1)*c)+{F^{-1}}(1-c))))
abline(h=0, lty=2)
if(doPDF) dev.off.pdf(file=file)

## check endpoints of objective function for root-finding
qrmtools:::Wang_obj(c(0, (1-alpha)/d), alpha=alpha, d=d, qF=qF) # -Inf, 0
## -Inf is not a problem for root finding, but the 0 at the right endpoint is.
## We take care of this by adjusting f.upper in the root finding


### 1.2) Checks for method="dual" ##############################################

## plot dual bound D(s) for various theta
## this concerns the outer root-finding
theta <- c(0.5, 1, 2, 4) # theta values
s <- seq(48, 2000, length.out=256) # s values
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

## investigating f(t) = D(s,t) - ..., the function for the inner
## root-finding to compute D(s); see dual_bound()
s <- c(1, 5, 10, 100, 500, 1000)
t <- sapply(seq_along(s), function(i) {
    res <- exp(seq(log(1e-3), log(s[i]/d), length.out=256))
    res[length(res)] <- s[i]/d # to avoid numerical issues (t > s/d)
    res
})
f <- sapply(seq_along(s), function(i)
            sapply(t[,i], function(t.)
                   dual_bound_2(s[i], t=t., d=d, pF=pF) -
                   dual_bound_2_deriv_term(s[i], t=t., d=d, pF=pF)))
palette <- colorRampPalette(c("red", "orange", "blue"), space="Lab")
cols <- palette(6)
if(doPDF)
    pdf(file=(file <- "fig_worst_VaR_0.99_hom_dual_D_st_Par=2_d=8.pdf"),
        width=6, height=6)
par(pty="s")
plot(t[,1], f[,1], type="l", log="x", xlim=range(t), ylim=range(f), col=cols[1],
     xlab="t", ylab=expression(f(s,t)))
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
## Conclusion: As we know, f(s, s/d) = 0. We also see that s has to be
##             sufficiently large in order to find a root f(s, t) = 0 for t < s/d


### 1.3) Checks for method="Pareto" ############################################

## TODO

### 1.4) Compare various methods of VaR_bound_hom() ############################

## TODO AM: quite different values; also check with Wang's method but Pareto case
wc.VaR.Wang <- VaR_bound_hom(alpha, d=d, method="Wang", qF=function(p) qPar(p, theta=th))
init <- VaR_bound_hom_crude(alpha, d=d, qF=qPar, theta=th)
wc.VaR.dual <- VaR_bound_hom(alpha, d=d, method="dual", interval=init, pF=pF)
##wc.VaR.Par <- VaR_bound_hom(alpha, d=d, method="Pareto", theta=th) # TODO AM => implement
##stopifnot(all.equal(wc.VaR.Wang, wc.VaR.dual, wc.VaR.Par))

