### GEV(xi, mu, sigma) distribution ############################################

##' @title Density of the GEV(xi, mu, sigma) distribution
##' @param x evaluation points
##' @param xi parameter xi (real)
##' @param mu parameter mu (real)
##' @param sigma parameter sigma (also real here; density is 0 if sigma <= 0)
##' @param log logical indicating whether the log density is computed
##' @return density of the GEV(xi, mu, sigma) distribution
##' @author Marius Hofert
dGEV <- function(x, xi, mu = 0, sigma = 1, log = FALSE)
{
    l <- length(x)
    if(sigma <= 0)
        return(rep(if(log) -Inf else 0, l)) # for logLik_GEV()
    y <- (x - mu) / sigma # sigma > 0
    if(xi == 0) { # xi == 0
        res <- -log(sigma) - (y + exp(-y))
        ## Note: - for mu >> x, y << 0 => exp(-y) = Inf => density correctly 0 then
        ##       - y should be >= -log(.Machine$double.xmax) for exp(-y) < Inf
    } else { # xi != 0
        res <- rep(-Inf, l) # correctly extend log-density
        xiy <- xi * y
        ii <- 1 + xiy > 0 # those indices for which density is positive
        res[ii] <- -log(sigma) + (-1/xi - 1) * log1p(xiy[ii]) - (1 + xiy[ii])^(-1/xi)
    }
    if(log) res else exp(res)
}

##' @title Distribution function of the GEV(xi, mu, sigma) distribution (vectorized in q)
##' @param q quantile
##' @param xi parameter xi
##' @param mu parameter mu
##' @param sigma parameter sigma
##' @param lower.tail logical indicating whether lower/upper tail is used
##' @param log.p logical indicating whether probabilities are given as log()
##' @return distribution function of the GEV(xi, mu, sigma) distribution
##' @author Marius Hofert
pGEV <- function(q, xi, mu = 0, sigma = 1, lower.tail = TRUE, log.p = FALSE)
{
    stopifnot(sigma > 0)
    y <- (q-mu)/sigma
    if(xi == 0) { # xi == 0
        if(lower.tail)
            if(log.p) -exp(-y) else exp(-exp(-y))
        else if(log.p) log1p(-exp(-exp(-y))) else 1-exp(-exp(-y))
    } else { # xi != 0
        xiy <- pmax(xi*y, -1) # see dGEV()
        if(lower.tail) {
            if(log.p) {
                -(1+xiy)^(-1/xi) # log H
            } else {
                exp(-(1+xiy)^(-1/xi)) # H
            }
        } else {
            if(log.p) {
                log1p(-exp(-(1+xiy)^(-1/xi))) # log(bar{H})
            } else {
                1-exp(-(1+xiy)^(-1/xi)) # bar{H}
            }
        }
        ## Formerly (also works):
        ## xiy <- xi*y
        ## ii <- 1+xiy > 0
        ## if(lower.tail) {
        ##     if(log.p) {
        ##         res <- if(xi < 0) rep(0, length(q)) else rep(-Inf, length(q)) # log('see case below')
        ##         res[ii] <- -(1+xiy[ii])^(-1/xi) # log H
        ##     } else {
        ##         res <- if(xi < 0) rep(1, length(q)) else rep(0, length(q)) # correctly extend; think of case '1+xi*x = 0' in both cases
        ##         res[ii] <- exp(-(1+xiy[ii])^(-1/xi)) # H
        ##     }
        ## } else {
        ##     if(log.p) {
        ##         res <- if(xi < 0) rep(-Inf, length(q)) else rep(0, length(q)) # log('see case below')
        ##         res[ii] <- log1p(-exp(-(1+xiy[ii])^(-1/xi))) # log(bar{H})
        ##     } else {
        ##         res <- if(xi < 0) rep(0, length(q)) else rep(1, length(q)) # 1 - 'case above'
        ##         res[ii] <- 1-exp(-(1+xiy[ii])^(-1/xi)) # bar{H}
        ##     }
        ## }
        ## res
    }
}

##' @title Quantile function of the GEV(xi, mu, sigma) distribution (vectorized in p)
##' @param p probability
##' @param xi parameter xi
##' @param mu parameter mu
##' @param sigma parameter sigma
##' @param lower.tail logical indicating whether lower/upper tail is used
##' @param log.p logical indicating whether probabilities are given as log()
##' @return quantile function of the GEV(xi, mu, sigma) distribution
##' @author Marius Hofert
qGEV <- function(p, xi, mu = 0, sigma = 1, lower.tail = TRUE, log.p = FALSE)
{
    stopifnot(sigma > 0)
    p <- if(log.p) pmin(p, 0) else pmin(pmax(p, 0), 1) # correctly extend
    if(xi == 0) { # xi == 0
        res <- if(lower.tail)
                   if(log.p) -log(-p) else -log(-log(p))
               else if(log.p) -log(-log1p(-exp(p))) else -log(-log1p(-p))
    } else { # xi != 0
          res <- if(lower.tail)
                     if(log.p) ((-p)^(-xi)-1)/xi else ((-log(p))^(-xi)-1)/xi
                 else if(log.p) ((-log1p(-exp(p)))^(-xi)-1)/xi else ((-log1p(-p))^(-xi)-1)/xi
          res <- if(xi < 0) pmin(res, -1/xi) else pmax(res, -1/xi)
    }
    mu + sigma * res
}

##' @title Generating random numbers from a GEV(xi, mu, sigma) distribution
##' @param n sample size n
##' @param xi parameter xi
##' @param mu parameter mu
##' @param sigma parameter sigma
##' @return n-vector containing GEV(xi, mu, sigma) random variates
##' @author Marius Hofert
rGEV <- function(n, xi, mu = 0, sigma = 1)
    qGEV(runif(n), xi = xi, mu = mu, sigma = sigma)


