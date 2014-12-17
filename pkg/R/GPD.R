### GPD(xi, beta) distribution #################################################

## Density of the GPD(xi, beta) distribution
## vectorized in x
dGPD <- function(x, xi, beta, log=FALSE)
{
    stopifnot(beta > 0)
    res <- rep(0, length(x)) # correctly extend
    if(xi == 0) { # xi == 0
        ind <- x>=0
        if(any(ind))
            res[ind] <- if(log) -x[ind]/beta-log(beta) else exp(-x[ind]/beta)/beta
    } else { # xi != 0
        ind <- if(xi > 0) x>=0 else 0<=x & x<=-beta/xi
        if(any(ind))
            res[ind] <- if(log) -(1/xi+1)*log1p(xi*x[ind]/beta)-log(beta)
            else (1+xi*x[ind]/beta)^(-(1/xi+1))/beta
    }
    res
}

## Distribution function of the GPD(xi, beta) distribution
## vectorized in q
pGPD <- function(q, xi, beta, lower.tail=TRUE, log.p=FALSE)
{
    stopifnot(beta > 0)
    if(xi == 0) { # xi == 0
        q <- pmax(q, 0) # correctly extend (instead of stopifnot(x >= 0))
        if(lower.tail)
            if(log.p) log1p(-exp(-q/beta)) else 1-exp(-q/beta)
        else if(log.p) -q/beta else exp(-q/beta)
    } else { # xi != 0
        q <- if(xi > 0) pmax(q, 0) else pmin(pmax(q, 0), -beta/xi) # correctly extend
        if(lower.tail)
            if(log.p) log1p(-(1+xi*q/beta)^(-1/xi)) else 1-(1+xi*q/beta)^(-1/xi)
        else if(log.p) -log1p(xi*q/beta)/xi else (1+xi*q/beta)^(-1/xi)
    }
}

## Quantile function of the GPD(xi, beta) distribution
## vectorized in p
qGPD <- function(p, xi, beta, lower.tail=TRUE, log.p=FALSE)
{
    stopifnot(beta > 0)
    p <- if(log.p) pmin(p, 0) else pmin(pmax(p, 0), 1) # correctly extend
    if(xi == 0) { # xi == 0
        if(lower.tail)
            if(log.p) (-beta)*log1p(-exp(p)) else (-beta)*log1p(-p)
        else if(log.p) (-beta)*p else (-beta)*log(p)
    } else { # xi != 0
        if(lower.tail)
            if(log.p) (beta/xi)*((-expm1(p))^(-xi)-1) else (beta/xi)*((1-p)^(-xi)-1)
        else if(log.p) (beta/xi)*expm1(-xi*p) else (beta/xi)*(p^(-xi)-1)
    }
}

## Random variate generation from the GPD(xi, beta) distribution
rGPD <- function(n, xi, beta)
    qGPD(runif(n), xi=xi, beta=beta)


### Par(theta) = GPD(1/theta, 1/theta), theta > 0 distribution #################

## Density of the Par(theta) distribution
dPar <- function(x, theta, log=FALSE)
    dGPD(x, xi=1/theta, beta=1/theta, log=log)

## Distribution function of the Par(theta) distribution
pPar <- function(q, theta, lower.tail=TRUE, log.p=FALSE)
    pGPD(q, xi=1/theta, beta=1/theta, lower.tail=lower.tail, log.p=log.p)

## Quantile function of the Par(theta) distribution
qPar <- function(p, theta, lower.tail=TRUE, log.p=FALSE)
    qGPD(p, xi=1/theta, beta=1/theta, lower.tail=lower.tail, log.p=log.p)

## Random variate generation from the Par(theta) distribution
rPar <- function(n, theta)
    qGPD(runif(n), xi=1/theta, beta=1/theta)

## Primitive of the Par(theta) survival function
bar_pPar_primitive <- function(q, theta)
    if(theta==1) log1p(q) else (1+q)^(1-theta) / (1-theta)

