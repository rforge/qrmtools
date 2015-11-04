## svn co svn+ssh://khornik@r-forge.r-project.org/svnroot/qrmtools


### Setup ######################################################################

require(qrmtools)

## Create input matrix X
alpha <- 0.99
d <- 1000
N <- 2^14 # = 16384
qF <- function(p, th=2) qPar(p, theta=th)
qF <- rep(list(qF), d)
p <- alpha + (1-alpha)*(0:(N-1))/N # for 'worst' (= largest) VaR
X <- sapply(qF, function(qF) qF(p))


### Profiling and measuring run time ###########################################

if(FALSE)
    Rprof(profiling <- tempfile(), line.profiling=TRUE)
system.time(res <- rearrange(X, tol=0, sample=FALSE, is.sorted=TRUE)) # ~ 23s (~ 14s with sample=TRUE) => sampling makes a big difference as ordering the values is the most time consuming part
if(FALSE) {
    Rprof(NULL) # disable profiling
    summaryRprof(profiling, lines="both")
}

## Result correct?
VaR_bounds_hom(alpha, d=d, method="Wang.Par", theta=2)[2]
res$bound


### Apply ARA() ################################################################

set.seed(271)
system.time(res. <- ARA(alpha, qF=qF, reltol=c(0.001, 0.01)))
res.$bounds
res.$N.used
## ~ 10s for d=1000, reltol=c(0.001, 0.01) (of course the N is different from above)
## ~ 20s for d=1000, reltol=c(0, 0.01)


### Is this sorting mechanism stable (on ties)? ################################

A <- matrix(c(1:4, 2*(1:4)-1, 2^(0:3)), ncol=3)
r <- rearrange(A, tol=NULL, sample=FALSE, is.sorted=TRUE, trace=TRUE)
## It first seems to be stable (see 2nd column on first rearrangement),
## but then it isn't (see 1st column on second rearrangement). Nevertheless,
## the procedure converges on all input matrices (see extra tester script)

B <- matrix(rep(1:3, 3), ncol=3)
r <- rearrange(B, tol=NULL, sample=FALSE, is.sorted=TRUE, trace=TRUE)

