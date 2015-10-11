## svn co svn+ssh://khornik@r-forge.r-project.org/svnroot/qrmtools
require(qrmtools)

alpha <- 0.99
d <- 128
N <- 2^14 # = 16384
qF <- function(p, th=2) qPar(p, theta=th)
method <- "worst"
maxiter <- Inf

qF <- rep(list(qF), d)
p <- if(method=="worst") alpha + (1-alpha)*(N-1):0/N else alpha*(N-1):0/N
X <- sapply(qF, function(qF) qF(p))

## rearrange()
Rprof(profiling <- tempfile())
res <- rearrange(X)
Rprof(NULL)
summaryRprof(profiling)

## ARA()
Rprof(profiling <- tempfile())
res <- ARA(alpha, qF=qF)
Rprof(NULL)
summaryRprof(profiling)



