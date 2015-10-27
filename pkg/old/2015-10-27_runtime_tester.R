## svn co svn+ssh://khornik@r-forge.r-project.org/svnroot/qrmtools
require(qrmtools)

alpha <- 0.99
d <- 1280 # 128 or 1280
N <- 2^14 # = 16384
qF <- function(p, th=2) qPar(p, theta=th)
method <- "worst"
maxiter <- Inf

qF <- rep(list(qF), d)
p <- if(method=="worst") alpha + (1-alpha)*(N-1):0/N else alpha*(N-1):0/N
X <- sapply(qF, function(qF) qF(p))

system.time(res <- rearrange(X))
## with C: 1.172s (11.19s for d=1280)
## with rank(): 1.18s (14.7s for d=1280)
## with order(order()): 1.48s (19.5s for d=1280)
## with R_orderVector(): 3.47s (34.13s for d=1280)
## with code from revision 55 (summer '15): 30s (4061s for d=1280) => C improvement factor: 95% (for d=1280: 99.7%)

## ARA()
set.seed(271)
system.time(res <- ARA(alpha, qF=qF))
## with C: 1.316s (12.8s for d=1280)
## with rank(): 1.70s (16.2s for d=1280)
## with order(order()): 2s (19.9s for d=1280)
## with R_orderVector(): 3.14s (31.27s for d=1280)

## Result:
res$bounds
## 2411.425 2431.257 for d=128
## 24114.26 24312.69 for d=1280
