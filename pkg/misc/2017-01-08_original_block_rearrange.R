
#' @title Rearrangement function for the ABRA
#'
#' @description The workhorse of the \code{\link{ABRA}()} function that carries out
#'              the rearrangement of the matrix of quantiles. It is a modified version of the block
#'              rearrangement function described in Bernard and McLeish (2014).
#'
#' @param X An \code{(N,d)}-matrix of quantiles (to be rearranged).
#' @param tol Absolute tolerance to determine (the individual) convergence. It needs to be
#'            greater or equal to 0.
#' @param M The lookback period, i.e. the number of steps (iterations) to look back for row sum
#'          variance when computing the (absolute) change in the row sum variance.
#' @param max.ra The maximal number of (considered) block rearrangements of the underlying
#'               matrix of quantiles (can be set to \code{Inf}).
#' @param method Indicates whether the bounds for the best or for the worst Value-at-Risk should be
#'               computed. These bounds are termed \underline\{s\}_N and \overline\{s\}_N in the
#'               literature (and below) and are theoretically not guaranteed bounds for (the best or the
#'               worst) Value-at-Risk; however, they are treated as such in practice and are typically in
#'               line with results from \code{VaR_bounds_hom()} in the homogeneous case, for example.
#' @param sample A \code{\link[base]{logical}} indicating whether each column of the two underlying
#'               matrices of quantiles (see Step 3 of the Rearrangement Algorithm in Embrechts et al.
#'               (2013)) are randomly permuted before the rearrangements begin.
#'
#' @return Returns a \code{\link[base]{list}} of the following elements:
#'         \item{X.rearranged}{The rearranged matrix \code{X}.}
#'         \item{bound}{The computed \underline\{s\}_N or \overline\{s\}_N (depending on the value
#'         of the \code{method} argument).}
#'         \item{tol}{The reached tolerance (i.e., the absolute change of the row sum variance
#'         after the last rearranged block of columns).}
#'         \item{converged}{A logical indicating whether the desired tolerance \code{tol} has been
#'         reached.}
#'         \item{iter}{The number of iterations (equal to the number of random partitions sampled)
#'         carried out before termination.}
#'
#' @details Unlike \code{rearrange()}, \code{block_rearrange()} uses row sum variance to determine
#'          its convergence. This has two implications. First, the target row sum variance is known
#'          in all cases and it is equal to 0, which makes the absolute tolerance preferable to the
#'          relative tolerance. Second, we do not have to change the convergence criterion when
#'          switching between the best and the worst Value-at-Risk.
#'
#'          Similarily to the \code{\link{rearrange}()}, \code{block_rearrange()} checks whether
#'          convergence has occurred after every iteration by comparing the absolute change to the
#'          row sum variance from \code{M} steps back to \code{tol}. While for the \code{rearrange()} we
#'          essentially have \code{M=d}, here for the \code{block_rearrange()} it is a free parameter and
#'          the main driver of the its runtime.
#' @references Bernard, C. and McLeish, D. (2014). Algorithms for Finding Copulas Minimizing Convex
#'             Functions of Sums. See \url{http://arxiv.org/abs/1502.02130v3}.
#'
#'             Embrechts, P., Puccetti, G., Rueschendorf, L. (2013). Model Uncertainty and VaR
#'             Aggregation. \emph{Journal of Banking & Finance} \strong{37}, 2750-2764.
#'
#'             Hofert, M., Memartoluie, A., Saunders, D. and Wirjanto, T. (2015). Improved Algorithms for
#'             Computing Worst Value-at-Risk: Numerical Challenges and the Adaptive Rearrangement Algorithm.
#'             See \url{http://arxiv.org/abs/1505.02281}.
#'          
#' @examples
#' # Set up a matrix of quantiles
#' N <- 500
#' d <- 4
#' alpha <- 0.99
#' qf <- function(x) qPar(x, theta = 4)
#' qfs <- rep(list(qf), d)
#' p <- alpha + (1 - alpha) * 0:(N - 1) / N
#' X <- sapply(qf, function(f) f(p))
#' 
#' # Rearrange a randomized initial matrix and compute a lower bound for the worst VaR
#' block_rearrange(X = X, tol = 0, M = 2, max.ra = Inf, method = "worst", sample = TRUE)
#' 
#' @author Martin Stefanik
#' @export
block_rearrange <- function(X, tol = 0, M = ncol(X), max.ra = Inf, method = c("worst", "best"), 
                            sample = TRUE)
{
  N <- nrow(X)
  d <- ncol(X)
  method <- match.arg(method)
  optim.fun <- if (method == "worst") min else max
  tol.fun <- function(x, y) abs(x - y)
  if (sample) X <- apply(X, 2, sample) else X
  X.rs <- .rowSums(X, m = N, n = d)
  rs.var <- rep_len(0, M)
  j <- 0
  iter <- 0

  while (TRUE) {
    iter <- iter + 1
    
    # Sample a random two-set partition
    n.cols <- sample(1:(d - 1), 1)
    cols <- sample(1:d, n.cols)
    rs.one <- .rowSums(X[, cols], m = N, n = n.cols)
    rs.two <- X.rs - rs.one
    
    # Rearrange one of the blocks
    X[, -cols] <- X[order(rs.two), ][order(order(rs.one, decreasing = TRUE)), -cols]
    
    # Assess convergence
    X.rs <- .rowSums(X, m = N, n = d)
    j <- if (j == M) 1 else j + 1
    if (iter == max.ra) {
      tol.reached <- FALSE
      tol. <- tol.fun(var(X.rs), rs.var.old)
      break
    } else if (iter > M) {
      rs.var.cur <- var(X.rs)
      rs.var.old <- rs.var[j]
      tol. <- tol.fun(rs.var.cur, rs.var.old)
      rs.var[j] <- rs.var.cur
      if (tol. <= tol) {
        tol.reached <- TRUE
        break
      } 
    } else {
      rs.var[j] <- var(X.rs)
    }
  }
  
  list(X.rearranged = X, bound = optim.fun(X.rs), tol = tol., converged = tol.reached, iter = iter)
}

