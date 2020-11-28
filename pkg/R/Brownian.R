### Simulate dependent Brownian and related motions ############################

##' @title Simulating Paths from Brownian Motions, Geometric Brownian Motions
##'        and Brownian Bridges
##' @param N number of paths
##' @param t n-vector (t_1,...,t_n) with 0 < t_1 < ... < t_n containing the time
##'        points at which to simulate
##' @param d dimension of the stochastic process
##' @param U (N * n, d)-matrix of copula realizations used as increments of the
##'        stochastic process.
##' @param drift d-vector of drifts mu (risk-neutral drifts: r - volas^2/2 for
##'        the risk-free interest rate r)
##' @param vola d-vector of volatilities sigma
##' @param type character vector indicating the type of stochastic process
##'        sampled (Brownian motion, geometric Brownian motion or Brownian
##'        bridge)
##' @param init initial values used if type = "GBM"
##' @return (N, n+1, d)-array
##' @author Marius Hofert
rBrownian <- function(N, t, d = 1, U = matrix(runif(N * n * d), ncol = d),
                      drift = rep(0, d), vola = rep(1, d),
                      type = c("BM", "GBM", "BB"), init = rep(1, d))
{
    ## Checks and setup
    stopifnot(N >= 1, (n <- length(t)) >= 1, t > 0, (diff(t. <- c(0, t))) > 0,
              d >= 1, (dim.U <- dim(U)) == c(N * n, d), length(drift) == d,
              length(vola) == d, vola > 0, length(init) == d, init > 0)

    ## Convert U to (N, n, d)-array of increments
    dmnms <- list("path" = NULL, "time" = NULL, "component" = NULL) # dimension names for returned arrays
    t. <- c(0, t)
    s <- outer(sqrt(diff(t.)), vola) # (n, d)-matrix scaling factors of N(0,1)
    Y <- array(, dim = c(N, n, d), dimnames = dmnms)
    for(i in 1:N)
        Y[i,,] <- s * qnorm(U[n * (i-1) + 1:n, ])

    ## Main
    type <- match.arg(type)
    switch(type,
           "BM" = {
               ## X_{t_k,j} = \mu_j t_k + \sum_{i=1}^{k} \sigma_j * \sqrt{t_{i}-t_{i-1}} * \Phi^{-1}(U_{i,j})
               X <- array(, dim = c(N, n+1, d), dimnames = dmnms)
               X[,1,] <- 0
               Y.csum <- matrix(0, nrow = N, ncol = d) # (N, d)-matrix
               for(k in 1:n) {
                   Y.csum <- Y.csum + Y[,k,]
                   X[,k+1,] <- rep(drift * t[k], each = N) + Y.csum
               }
               X
           },
           "GBM" = {
               ## S_{t_k,j} = S_{0,j} * e^{X_{t_k,j}}
               S <- array(, dim = c(N, n+1, d), dimnames = dmnms)
               S[,1,] <- rep(init, each = N) # to be treated as an (N, d)-matrix
               Y.csum <- matrix(0, nrow = N, ncol = d) # (N, d)-matrix
               for(k in 1:n) {
                   Y.csum <- Y.csum + Y[,k,]
                   S[,k+1,] <- S[,1,] * exp(rep(drift * t[k], each = N) + Y.csum)
               }
               S
           },
           "BB" = {
               ## B_{t_k,j} = X_{t_k,j} - (t/T) X_{t_n,j}
               W <- rBrownian(N, t = t, d = d, U = U, drift = drift, vola = vola,
                              type = "BM") # (N, n+1, d)-array
               t.T <- t. / t[n] # (n+1)-vector (0, t_1/T, ..., t_n/T = 1)
               W.T <- W[,n+1,] # (N, d)-matrix
               B <- array(, dim = c(N, n+1, d), dimnames = dmnms)
               for(k in 1:(n+1)) { # by far the easiest version; both sapply() and array() require reorderings
                   B[,k,] <- W[,k,] - t.T[k] * W.T
               }
               B
           })
}
