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

## Note: - hard-coded here to be vectorized in the main argument and theta
##       - F(x) = 1-(1+x)^{-theta}

## Density of the Par(theta) distribution
dPar <- function(x, theta, log=FALSE)
    if(log) log(theta)-(theta+1)*log1p(x) else theta*(1+x)^(-theta-1)

## Distribution function of the Par(theta) distribution
pPar <- function(q, theta, lower.tail=TRUE, log.p=FALSE)
    if(lower.tail)
        if(log.p) log(1-(1+x)^(-theta)) else 1-(1+x)^(-theta)
    else if(log.p) -theta*log1p(x) else (1+x)^(-theta)

## Quantile function of the Par(theta) distribution
qPar <- function(p, theta, lower.tail=TRUE, log.p=FALSE)
    if(lower.tail)
        if(log.p) (-expm1(p))^(-1/theta)-1 else (1-p)^(-1/theta)-1
    else if(log.p) expm1(-p/theta) else p^(-1/theta)-1

## Random variate generation from the Par(theta) distribution
rPar <- function(n, theta) qPar(runif(n), theta=theta)

## Primitive of the Par(theta) survival function
bar_pPar_primitive <- function(q, theta)
    if(theta==1) log1p(q) else (1+q)^(1-theta) / (1-theta)

