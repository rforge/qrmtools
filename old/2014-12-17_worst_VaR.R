## removed from demo worst_VaR:

## check D(s,t)
qrmtools:::dual_bound_2(s=32, t=1, d=d, pF=pF) # works
qrmtools:::dual_bound_2(s=0.5, t=0, d=d, pF=pF) # works
db.up  <- qrmtools:::dual_bound_2(s=100, t=100/d, d=d, pF=pF) # works at upper end point for t
db.up. <- d*pPar(100/d, theta=th, lower.tail=FALSE) # theoretical value at upper end point
stopifnot(all.equal(db.up, db.up.)) # => same

## plot dual bound D(s)
s <- seq(0, 1000, length.out=256)
D <- sapply(s, function(s.) dual_bound(s., d=d, pF=pF))
if(doPDF)
    pdf(file=(file <- "fig_wc_VaR_0.99_hom_dual_D_s_Par=2_d=8.pdf"),
        width=6, height=6)
par(pty="s")
plot(s, D, type="l", log="y", ylab="Dual bound D(s)")
abline(h=1-alpha, lty=2)
legend("topright", inset=0.02, lty=1:2, col=c("black", "black"),
       bty="n", legend=c("D(s)", expression(1-alpha)), y.intersp=1.2)
if(doPDF) dev.off.pdf(file=file)
