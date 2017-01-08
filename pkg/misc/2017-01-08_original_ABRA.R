
#' @title The Adaptive Block Rearrangement Algorithm (ABRA)
#'
#' @description A function that carries out an adaptive version of an optimized version of the Block
#'              Rearrangement Algorithm (BRA) from Bernard and McLeish (2014).
#'
#' @param alpha Value-at-Risk confidence level (e.g., 0.99).
#' @param qF A list of lenght \code{d} (the number of dimensions) containing the marginal quantile
#'           functions.
#' @param N.exp The exponents of the number of discretization points (a vector) over which the
#'              algorithm iterates to find the smallest number of discretization points for which the
#'              desired accuracy (specified by \code{reltol}) is attained; for each number of
#'              discretization points, at most \code{max.ra}-many column rearrangements of the
#'              underlying matrix of quantiles are considered.
#' @param M The lookback period, i.e. the number of steps (iterations) to look back for row sum
#'          variance when computing the (absolute) change in the row sum variance.
#' @param max.ra The maximal number of (considered) block rearrangements of the underlying
#'               matrix of quantiles (can be set to \code{Inf}).
#' @param method Indicates whether the bounds for the best or for the worst Value-at-Risk should be
#'               computed. These bounds are termed \underline\{s\}_N and \overline\{s\}_N in the l
#'               iterature (and below) and are theoretically not guaranteed bounds of (best or worst)
#'               Value-at-Risk; however, they are treated as such in practice and are typically in line
#'               with results from \code{VaR_bounds_hom()} in the homogeneous case, for example.
#' @param tol A vector of length two containing the individual absolute convergence and the joint
#'            relative convergence tolerances. The individual absolute tolerance is used to determine
#'            convergence of the row sum variances for \underline\{s\}_N and \overline\{s\}_N, and the
#'            joint relative tolerance (i.e., the relative difference between the computed
#'            \underline\{s\}_N and \overline\{s\}_N with respect to \overline\{s\}_N) is used to determine
#'            the convergence of \underline\{s\}_N and \overline\{s\}_N.
#' @param sample A \code{\link[base]{logical}} indicating whether each column of the two underlying
#'               matrices of quantiles (see Step 3 of the Rearrangement Algorithm in Embrechts et al.
#'               (2013)) are randomly permuted before the rearrangements begin.
#'
#' @return Returns a \code{\link[base]{list}} of the following elements:
#'         \item{bounds}{A bivariate vector containing the computed \underline\{s\}_N and \overline\{s\}_N
#'         (the so-called rearrangement range) which are typically treated as bounds for (the worst or
#'         the best) Value-at-Risk.}
#'         \item{rel.ra.gap}{The reached relative tolerance (also known as relative rearrangement gap)
#'         between \underline\{s\}_N and \overline\{s\}_N computed with respect to \overline\{s\}_N.}
#'         \item{tol}{A trivariate \code{\link[base]{vector}} containing the reached individual
#'         absolute tolerances and the reached joint relative tolerance.}
#'         \item{converged}{A trivariate logical vector indicating convergence of the computed
#'         \underline\{s\}_N and \overline\{s\}_N as well as their joint relative covergence.}
#'         \item{N.used}{The actual \code{}} used for computing the (final) \underline\{s\}_N and
#'         \overline\{s\}_N.}
#'         \item{X.rearranged}{The rearranged matrices for \underline{s}_N and \overline{s}_N with
#'         which the \code{ABRA()} is terminated.}
#'
#' @details The most important difference to the \code{ARA()} aside from the fact that
#'          \code{block_rearrange()} is used to rearrange the initial matrices of quantiles is the
#'          implicit selection of the value for the lookback period \code{M} using an externally fitted
#'          Poisson GLM (see \code{\link[stats]{glm}}). See also the details for
#'          \code{\link{block_rearrange}}.
#'
#' @references Bernard, C. and McLeish, D. (2014). Algorithms for Finding Copulas Minimizing Convex
#'             Functions of Sums. See \url{http://arxiv.org/abs/1502.02130v3}.
#'
#'             Embrechts, P., Puccetti, G., Rueschendorf, L. (2013). Model Uncertainty and VaR
#'             Aggregation. \emph{Journal of Banking & Finance} \strong{37}, 2750-2764.
#'
#'             Hofert, M., Memartoluie, A., Saunders, D. and Wirjanto, T. (2015). Improved Algorithms for
#'             Computing Worst Value-at-Risk: Numerical Challenges and the Adaptive Rearrangement Algorithm.
#'             See \url{http://arxiv.org/abs/1505.02281}.
#' @examples 
#' # Initial setup
#' alpha <- 0.99
#' d <- 5
#' qf <- function(x) qnorm(x, mean = 10, sd = 3)
#' qfs <- rep(list(qf), d)
#' ABRA(alpha = alpha, qf = qfs, M = NULL, method = "worst", tol =c(0, 0.01), max.ra = Inf,
#'      sample = TRUE)
#' 
#' @author Martin Stefanik
#' @export
ABRA <- function(alpha, qF, N.exp = seq(8, 19, by = 1), M = NULL, method = c("worst", "best"), 
                 tol = c(0, 0.01), max.ra = Inf, sample = TRUE)
{
  ltol <- length(tol)
  stopifnot(0 < alpha, alpha < 1, ltol == 1 || ltol == 2, tol >= 0, 
            length(N.exp) >= 1, N.exp >= 1, is.logical(sample), is.list(qF),
            sapply(qF, is.function), (d <- length(qF)) >= 2)
  method <- match.arg(method)
  itol <- if (ltol == 2) tol[1] else 0
  jtol <- if (ltol == 2) tol[2] else tol[1]
  
  # Function that returns the optimal M
  m.opt <- function(N, d, col.var) {
    b0 <-  1.6482635640
    b1 <- -0.0013742363
    b2 <-  0.0112121293
    b3 <-  0.0001273265
    ceiling(exp(b0 + b1 * N + b2 * d + b3 * col.var))
  }
  
  for (N in 2^N.exp) {
    
    # Create the discretization matrix for the approximation from below
    p <- if (method == "worst") {
      alpha + (1 - alpha) * (0:(N - 1)) / N
    } else {
      alpha * (0:(N - 1)) / N
    }
    
    X.low <- sapply(qF, function(qF) qF(p))
    if (method == "best") 
      X.low[1, ] <- sapply(1:d, function(j) if (is.infinite(X.low[1, j])) {
        qF[[j]](alpha / (2 * N))
      } else {
        X.low[1, j]
      })
    
    # Rearrange for the approxiation from below
    M.current <- if (is.null(M)) m.opt(N, d, mean(apply(X.low, 2, var))) else M
    res.low <- rearrangeABRA(X.low, M = M.current, tol = itol,
                             max.ra = max.ra, sample = sample)
    
    # Create the discretization matrix for the approximation from above
    p <- if (method == "worst") alpha + (1 - alpha) * (1:N) / N else alpha * (1:N) / N
    X.up <- sapply(qF, function(qF) qF(p))
    if (method == "worst") 
      X.up[N, ] <- sapply(1:d, function(j) if (is.infinite(X.up[N, j])) {
        qF[[j]](alpha + (1 - alpha) * (1 - 1 / (2 * N)))
      } else {
        X.up[N, j]
      })
    
    # Rearrange for the approxiation from above
    M.current <- if (is.null(M)) m.opt(N, d, mean(apply(X.up, 2, var))) else M
    res.up <- rearrangeABRA(X.up, M = M.current, tol = itol, 
                            max.ra = max.ra, sample = sample)
    
    # Check joint convergence
    joint.tol <- abs((res.low$bound - res.up$bound) / res.up$bound)
    joint.tol.reached <- joint.tol <= jtol
    if (res.low$converged && res.up$converged && joint.tol.reached) break
  }
  
  list(bounds = c(low = res.low$bound, up = res.up$bound), 
       rel.ra.gap = abs((res.up$bound - res.low$bound) / res.up$bound), 
       tol = c(low = res.low$tol, up = res.up$tol, joint = joint.tol), 
       converged = c(low = res.low$converged, up = res.up$converged, 
                     joint = joint.tol.reached),
       N.used = N, num.iter = c(low = res.low$iter, up = res.up$iter), 
       X.rearranged = list(low = res.low$X.rearranged, up = res.up$X.rearranged))
}
