\name{worst_VaR}
\alias{crude_VaR_bounds}
\alias{worst_VaR_hom}
\alias{dual_bound}
\alias{rearrange}
\alias{RA}
\alias{ARA}
\title{Worst Value-at-Risk for Given Margins}
\description{
  Compute the worst Value-at-Risk for given marginal distributions under
  various setups.
}
\usage{
crude_VaR_bounds(alpha, d, qF, ...)
worst_VaR_hom(alpha, d, method=c("Wang", "Wang.Par", "Wang.Par.trafo", "dual"),
              interval=NULL, tol=NULL, ...)
dual_bound(s, d, pF, tol=.Machine$double.eps^0.25, ...)

rearrange(X, tol, tol.type=c("relative", "absolute"), maxiter=Inf,
          method=c("worst", "best"), impl=c("R", "C"))
RA(alpha, d, qF, N, abstol=NULL, maxiter=Inf, method=c("worst", "best"),
   sample=TRUE, impl=c("R", "C"))
ARA(alpha, d, qF, N=2^seq(8, 20, by=1), reltol=c(0.001, 0.01), maxiter=10,
    method=c("worst", "best"), sample=TRUE, impl=c("R", "C"))
}
\arguments{
  \item{alpha}{Value-at-Risk confidence level (e.g., 0.99).}
  \item{d}{Dimension (number of risk factors; \eqn{\ge 2}{>=2}).}
  \item{qF}{The marginal quantile function in the homogeneous case or,
    for \code{crude_VaR_bounds()} and \code{RA()}, a list (of length
    \code{d}) with the marginal quantile functions.}
  \item{method}{A \code{\link{character}} string. For
    \describe{
      \item{\code{worst_VaR_hom()}:}{\code{method="Wang"},
	\code{method="Wang.Par"} and \code{method="Wang.Par.trafo"}
	apply the approach of Embrechts et al. (2014,
	Proposition 3.1). The latter two assume Pareto margins and thus do
	not require numerical integration;
          \code{method="Wang.Par.trafo"} first transforms the
          root finding problem to a different scale.
	\code{method="dual"} applies the
	dual bound approach as in Embrechts et al. (2013, Proposition 4).}
      \item{\code{rearrange()}, \code{RA()},
        \code{ARA()}:}{\code{method} indicates whether bounds
	for the worst (i.e., largest) or the best (i.e., smallest)
	Value-at-Risk are computed.}
  }}
  \item{interval}{Initial interval (a \code{\link{numeric}(2)}).
    provided, these are the defaults chosen:
    \describe{
      \item{\code{method="Wang"}:}{The initial
	interval is \eqn{[0,(1-\alpha)/d]}{[0,(1-alpha)/d]}.}
      \item{\code{method="Wang.Par"}:}{The initial
	interval is \eqn{[c_l,(1-\alpha)/d]}{[c_l,(1-alpha)/d]},
	where \eqn{c_l} is chosen as in Hofert et al. (2015).}
      \item{\code{method="Wang.Par.trafo"}:}{The initial
	interval is \eqn{[1,x_u]},
	where \eqn{x_u} is chosen as in Hofert et al. (2015).}
      \item{\code{method="dual"}:}{In this case, no good defaults are known.
	Note that the lower endpoint of the initial interval has to be
	sufficiently large in order for the the inner root-finding algorithm
	to find a root; see Details.}
  }}
  \item{tol}{
    \describe{
      \item{\code{worst_VaR_hom()}:}{The tolerance for the involved
	root-finding algorithm as passed
	to \code{\link{uniroot}()}. This defaults to
	\eqn{2.2204*10^{-16}}{2.2204*10^{-16}}
	for \code{method="Wang"} or \code{method="Wang.Par"} (where a
	smaller tolerance is crucial) and to \code{\link{uniroot}()}'s
	default \code{.Machine$double.eps^0.25} otherwise. Note that for
	\code{method="dual"}, \code{tol} is used for both the outer
	and the inner root-finding procedure.}
      \item{\code{rearrange()}:}{Tolerance to determine the (individual)
	convergence; if \code{NULL}, the iteration is done until the
	matrix does not change anymore.}
  }}
  \item{tol.type}{\code{\link{character}} string indicating the
    type of convergence tolerance function to be used (\code{"relative"}
    for relative tolerance and \code{"absolute"} for absolute tolerance).}
  \item{s}{Dual bound evaluation point.}
  \item{pF}{The marginal loss distribution function.}
  \item{X}{An (\code{N}, \code{d})-matrix containing
    \eqn{N}-discretizations of the \eqn{d} quantiles to be rearranged.}
  \item{maxiter}{Maximal number of iterations over the columns of the
    underlying matrix of quantiles (can be set to \code{Inf}).}
  \item{impl}{\code{\link{character}} string providing the type of
    implementation used (C vs R).}
  \item{N}{Number of discretization points. For
    \describe{
      \item{\code{RA()}:}{A positive integer.}
      \item{\code{ARA()}:}{A vector of positive integers. For each
	such number, at most \code{maxiter}-many iterations through the columns
	of the underlying matrix of quantiles are conducted, checking after
	each of these iterations whether the desired relative tolerances
	(both individually and jointly; see below) has been
	achieved.}
  }}
  \item{abstol}{Absolute convergence tolerance \eqn{\epsilon}{epsilon}
    (to determine the individual convergence, i.e., the change in the computed minimal
    (for worst Value-at-Risk) or maximal (for best Value-at-Risk)
    row sums for the lower and upper bound). If \code{NULL} (the
    default) and \code{maxiter} has not been reached, the iterations
    through the columns of the underlying matrix of quantiles is done
    until the matrix does not change anymore, i.e., until each column
    is oppositely ordered to the sum of all others.}
  \item{reltol}{A vector containing two relative convergence tolerances,
    the individual relative tolerance (first component; used to determine
    convergence of the minimal (for worst Value-at-Risk) or maximal
    (for best Value-at-Risk) row sums for the lower and upper bound)
    and the joint relative tolerance (second component; relative tolerance
    between the computed lower and upper bound computed with respect
    to the upper bound). Note that \code{reltol[1]} is allowed to be
    \code{NULL}; see \code{abstol} above for how to interpret this case.}
  \item{sample}{\code{\link{logical}} indicating whether each column of
    the two internal matrices of quantiles (see Step 3 of the Rearrangement
    Algorithm in Embrechts et al. (2013))
    are sampled before the iteration begins.}
  \item{\dots}{
    \describe{
      \item{\code{worst_VaR_hom()}:}{For \code{d=2} \dots must contain the quantile
        function \code{qF()}. The same holds for \code{method="Wang"},
        but additional arguments passed via \dots are passed on to
	% the underlying objective function Wang_h() from which they
	% are passed on to
	the underlying \code{\link{integrate}()}. For
        \code{method="Wang.Par"} and \code{method="Wang.Par.trafo"},
	\dots must contain the Pareto parameter \eqn{\theta>0}{theta>0}.
	For \code{method="dual"}, \dots must contain the distribution
	function \code{pF()}, the initial interval \code{interval} for
	the outer root finding procedure; additional arguments
	are passed on to the underlying
	\code{\link{integrate}()} for computing the dual bound
	\eqn{D(s)}.}
      \item{\code{dual_bound()}:}{The \dots arguments are passed to the
	underlying \code{\link{integrate}()}.}
      \item{\code{crude_VaR_bounds()}:}{The \dots arguments are passed
        to (all provided) quantile functions.}
    }}
}
\value{
  \code{crude_VaR_bounds()} returns crude lower and upper bounds for
  Value-at-Risk at confidence level \eqn{\alpha}{alpha} for any
  \eqn{d}-dimensional model with marginal quantile functions
  specified by \code{qF}.

  \code{worst_VaR_hom()} returns the worst (i.e., largest) Value-at-Risk at
  confidence level \eqn{\alpha}{alpha} for \eqn{d} risks with equal
  distribution function specified by the ellipsis \code{...}.

  \code{dual_bound()} returns the value of the dual bound \eqn{D(s)} as
  given in Embrechts, Puccetti, \enc{Rüschendorf}{Rueschendorf}
  (2013, Eq. (12)).

  \code{rearrange()} returns a \code{\link{list}} containing
  \describe{
    \item{\code{bound}:}{The computed (lower
      or upper) bound for (worst or best) Value-at-Risk.}
    \item{\code{err}:}{The tolerance (i.e., the (absolute or
      relative) change of the minimal (for worst Value-at-Risk) or
      maximal (for best Value-at-Risk) row sum in the last iteration).}
    \item{\code{num.iter}:}{The number of iterations over the
      columns of the matrix of quantiles.}
    \item{\code{row.sums}:}{A \code{\link{matrix}} with the same number of rows as
      \code{X}, giving the row sums after each iteration of
      rearrangements of \code{X}; one iteration corresponds to looping
      over all columns and oppositely reorder each of them with respect
      to the sum of all others. The number of columns of \code{row.sums}
      thus corresponds to \code{num.iter}.}
    \item{\code{m.row.sums}:}{A \code{\link{numeric}} containing minimal
      (for worst Value-at-Risk) or maximal (for best Value-at-Risk)
      row sums for the computed bound.}
    \item{\code{num.opp.ordered}:}{The number of oppositely ordered
      columns.}
  }

  \code{RA()} returns a \code{\link{list}} containing
  \describe{
    \item{\code{bounds}:}{A bivariate vector containing the computed lower
      and upper bound for (worst or best) Value-at-Risk; see also
      \code{bound} returned by \code{rearrange()} above.}
    \item{\code{rel.DU.spread}:}{relative tolerance (computed with respect
      to the upper bound) between lower and upper bound.}
    \item{\code{individual.err}:}{A bivariate vector containing the
      individual tolerances (i.e., the absolute change of the minimal
      (for worst Value-at-Risk) or maximal (for best Value-at-Risk)
      row sums for the lower (first component) and upper bound (second
      component) in the last iteration of \code{rearrange()}); see also
      \code{err} returned by \code{rearrange()} above.}
    \item{\code{num.iter}:}{A bivariate vector containing, for the lower
      and upper bound, the number of iterations over the
      columns of the underlying matrix of quantiles; see also
      \code{num.iter} returned by \code{rearrange()} above.}
    \item{\code{row.sums}:}{A list of length two containing the computed
      matrices of row sums (one column for each iteration) for the lower
      and upper bound; see also \code{rearrange()} above; see also
      \code{row.sums} returned by \code{rearrange()} above.}
    \item{\code{m.row.sums}:}{A list of length two containing minimal
      (for worst Value-at-Risk) or maximal (for best Value-at-Risk)
      row sums for the computed lower and upper bounds; see also \code{m.row.sums}
      returned by \code{rearrange()} above.}
    \item{\code{num.opp.ordered}:}{A bivariate vector containing
      the number of oppositely ordered columns for the lower and upper
      bound; see also \code{num.opp.ordered} returend by
      \code{rearrange()} above.}
  }

  \code{ARA()} returns a \code{\link{list}} containing
  \describe{
    \item{\code{bounds}:}{See \code{RA()}.}
    \item{\code{rel.DU.spread}:}{See \code{RA()}.}
    \item{\code{individual.err}:}{See \code{RA()}.}
    \item{\code{joint.err}:}{The joint relative tolerance reached (i.e.,
      the relative tolerance between the computed lower and upper bound
      computed with respect to the upper bound).}
    \item{\code{N.used}:}{The actual \code{N} used for computing
      the (final) lower and upper bound.}
    \item{\code{num.iter}:}{See \code{RA()}.}
    \item{\code{row.sums}:}{See \code{RA()}.}
    \item{\code{m.row.sums}:}{See \code{RA()}.}
    \item{\code{num.opp.ordered}:}{See \code{RA()}.}
  }
}
\details{
  For \code{d=2}, \code{worst_VaR_hom()} uses the method of
  Embrechts et al. (2013,
  Proposition 2). For \code{method="Wang"}, \code{method="Wang.Par"} and
  \code{method="Wang.Par.trafo"},
  the method presented in Embrechts et al. (2014, Proposition 3.1) is
  implemented. This requires
  one \code{\link{uniroot}()} and, for the generic \code{method="Wang"},
  one \code{\link{integrate}()}. The critical part for the
  generic \code{method="Wang"} is the lower endpoint of the initial
  interval for \code{\link{uniroot}()}. If the (marginal)
  distribution function has finite first moment, this can be taken as
  0. However, if it has infinite first moment, the lower endpoint has to
  be positive (below the unknown root and smaller than the upper
  endpoint \eqn{(1-\alpha)/d}{(1-alpha)/d}) (which also happens to be a
  root). In the case of Pareto margins, Hofert et al. (2015) have
  derived such an initial lower endpoint (which is used by
  \code{method="Wang.Par"} and \code{method="Wang.Par.trafo"}; for the
  latter, due to the transformation of the objective function involved,
  this is actually an initial upper endpoint, see Hofert et al. (2015)
  for more details. Also note that the chosen smaller default tolerances for
  \code{\link{uniroot}()} in case of \code{method="Wang"} and
  \code{method="Wang.Par"} are crucial for obtaining reliable
  Value-at-Risk values; see Hofert et al. (2015).

  For \code{method="dual"}, the method presented of Embrechts et
  al. (2013, Proposition 4) is implemented. This requires two (nested)
  \code{\link{uniroot}()}, and an \code{\link{integrate}()}. For the
  inner root-finding procedure to find a root, the lower endpoint of the
  provided initial \code{interval} has to be \dQuote{sufficiently
    large}.

  Note that these approaches for computing the worst (i.e., largest)
  Value-at-Risk in the homogeneous case are numerically non-trivial;
  see the source code and \code{demo(worst_VaR)} for more details. As a
  rule of thumb, use \code{method="Wang"} if you have to (i.e., if the
  margins are not Pareto) and \code{method="Wang.Par"} if you can (i.e.,
  if the margins are Pareto); in the former case, it might be wise to
  follow along the lines of \code{method="Wang.Par.trafo"} to implement
  one's own method (e.g., if \code{method="Wang"} turns out to be
  numerically critical), see Hofert et al. (2015) for more details.
  Furthermore, it is not recommended to use (the more challenging)
  \code{method="dual"}.


  Concerning the inhomogeneous case, \code{rearrange()} is an auxiliary
  function which is called by \code{RA()} and \code{ARA()}. For
  performance reasons, no checking is done. It should only be used by
  experts.

  For the Rearrangement Algorithm \code{RA()}, unless
  \code{maxiter}-many iterations have been reached, convergence is
  determined for the lower and the upper bound if the minimal
  (for the worst Value-at-Risk) or maximal (for the best
  Value-at-Risk) row sums after an iteration did not change by more than
  (the absolute tolerance) \code{abstol} (so \eqn{\le\epsilon}{<= eps}). This
  is slightly different from Embrechts et al. (2013) who use \eqn{<\epsilon}{< eps}
  but this slight modification has the advantage of allowing
  \eqn{\epsilon=0}{eps=0}.

  For the Adaptive Rearrangement Algorithm \code{ARA()},
  again unless \code{maxiter}-many iterations
  have been reached, convergence is determined if the minimal
  (for the worst Value-at-Risk) or maximal (for the best
  Value-at-Risk) row sums after an iteration did not change by more than
  (the individual relative tolerance) \code{reltol[1]} \emph{and} the
  relative (joint) tolerance between both bounds is at most \code{reltol[2]}.

  Note that both \code{RA()} and \code{ARA()} need to evalute the
  0-quantile (for the lower bound for the best Value-at-Risk) and
  the 1-quantile (for the upper bound for the
  worst Value-at-Risk). As the algorithms can only handle finite values, the
  0-quantile and the 1-quantile need to be adjusted if infinite. Instead
  of the 0-quantile, the \eqn{\alpha/(2N)}{alpha/(2N)}-quantile is
  computed and instead of the 1-quantile the
  \eqn{\alpha+(1-\alpha)(1-1/(2N))}{alpha+(1-alpha)(1-1/(2N))}-quantile
  is computed for such margins (if the 0-quantile or the 1-quantile is
  finite, no adjustment is made).

  As a rule of thumb (see the examples below,
  \code{demo(worst_VaR)} and Hofert et al. (2015) for the reasons),
  it is recommended to use \code{ARA()} instead of \code{RA()}.
}
\author{Marius Hofert}
\references{
  Embrechts, P., Puccetti, G., \enc{Rüschendorf}{Rueschendorf}, L.,
  Wang, R., Beleraj, A. (2014).  An Academic Response to Basel
  3.5. \emph{Risks} \bold{2}(1), 25--48.

  Embrechts, P., Puccetti, G., \enc{Rüschendorf}{Rueschendorf}, L. (2013).
  Model uncertainty and VaR aggregation. \emph{Journal of Banking \&
    Finance} \bold{37}, 2750--2764.

  Hofert, M., Memartoluie, A., Saunders, D., Wirjanto, T. (2015).
  Improved Algorithms for Computing Worst Value-at-Risk:
  Numerical Challenges and the Adaptive Rearrangement Algorithm.
  See http://arxiv.org/abs/1505.02281.
}
\seealso{
  \code{demo(worst_VaR)} for more example calls, numerical challenges
  encoutered and a comparison of the different methods for computing
  the worst (i.e., largest) Value-at-Risk.
}
\examples{
require(qrmtools)
## Pareto setup
alpha <- 0.99 # VaR confidence level
th <- 2 # Pareto parameter theta
qF <- function(p, theta=th) qPar(p, theta=theta) # Pareto quantile function
pF <- function(q, theta=th) pPar(q, theta=theta) # Pareto distribution function

## d=2: Compute worst VaR explicitly (homogeneous case, 'method' is ignored)
resW <- worst_VaR_hom(alpha, d=2, qF=qF)
set.seed(271) # set seed (for reproducibility)
resARA <- ARA(alpha, d=2, qF=qF) # Adaptive Rearrangement Algorithm (ARA)
resRA <- RA(alpha, d=2, qF=qF, N=256, abstol=0) # RA with N and abstol chosen to match ARA
stopifnot(all.equal(resW, resARA$bounds[2], resRA$bounds[2], tolerance=0.005))

## d=8: Compare the various methods (in the homogeneous case)
d <- 8 # dimension
I <- crude_VaR_bounds(alpha, d=d, qF=qF) # crude bound
(w.VaR.Wang           <- worst_VaR_hom(alpha, d=d, method="Wang", qF=qF))
(w.VaR.Wang.Par       <- worst_VaR_hom(alpha, d=d, method="Wang.Par", theta=th))
(w.VaR.Wang.Par.trafo <- worst_VaR_hom(alpha, d=d, method="Wang.Par.trafo", theta=th))
(w.VaR.dual           <- worst_VaR_hom(alpha, d=d, method="dual", interval=I, pF=pF))
stopifnot(all.equal(c(w.VaR.Wang, w.VaR.Wang.Par, w.VaR.Wang.Par.trafo),
                    rep(w.VaR.dual, 3), tolerance=5e-6))
## Note: The same check is passed for d <- 100

## Compare with ARA (with defaults)
set.seed(271)
resARA <- ARA(alpha, d=d, qF=qF)
resARA$bounds
## Compare with RA
## 1) Choose N and abstol chosen to match ARA results
abstol <- tail(abs(diff(resARA$m.row.sums[[2]])), n=1)
resRA <- RA(alpha, d=d, qF=qF, N=512, abstol=abstol)
resRA$bounds
## 2) Iterate until each column is opp. ordered to the sum of all others
resRA. <- RA(alpha, d=d, qF=qF, N=512) # abstol=NULL
resRA.$bounds

## Using the additional results computed by RA(); for resRA.
xlim <- c(1, max(sapply(resRA.$m.row.sums, length)))
ylim <- range(resRA.$m.row.sums)
plot(resRA.$m.row.sums[[2]], type="b", xlim=xlim, ylim=ylim, xlab="Iteration",
     ylab=paste0("Minimal row sum (# of oppositely ordered columns: lower = ",
     resRA.$num.opp.ordered[1],", upper = ",
     resRA.$num.opp.ordered[2],")"),
     main=substitute("Rearrangement Algorithm for"~alpha==a.*","~d==d.*" and Par("*
                     th.*")", list(a.=alpha, d.=d, th.=th)))
lines(1:length(resRA.$m.row.sums[[1]]), resRA.$m.row.sums[[1]], type="b",
      col="blue")
legend("bottomright", bty="n", lty=rep(1,2),
       col=c("black", "blue"), legend=c("upper bound", "lower bound"))
## => One should use ARA() instead of RA()

## "Reproducing" an example from Embrechts et al. (2013)
## (seed unknown, eps unknown)
set.seed(271)
(resRA <- RA(alpha, d=3, qF=qF, N=50))
stopifnot(all.equal(resRA$bounds[1], 45.11593, tolerance=1e-7),
          all.equal(resRA$bounds[2], 46.40413, tolerance=1e-7),
          resRA$num.opp.ordered == c(3,3))

## An example from Embrechts et al. (2013) with larger N
## (first two decimals are equal)
if(FALSE) {
    resRA <- RA(alpha, d=3, qF=qF, N=1e5)
    stopifnot(all.equal(resRA$bounds, rep(45.98, 2), tolerance=3e-4))
}
}
\keyword{distribution}