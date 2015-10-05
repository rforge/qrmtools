\name{VaR_bounds}
\alias{crude_VaR_bounds}
\alias{VaR_bounds_hom}
\alias{dual_bound}
\alias{rearrange}
\alias{RA}
\alias{ARA}
\title{Best and Worst Value-at-Risk for Given Margins}
\description{
  Compute the best and worst Value-at-Risk for given marginal distributions.
}
\usage{
crude_VaR_bounds(alpha, d, qF, ...)
VaR_bounds_hom(alpha, d, method=c("Wang", "Wang.Par", "Wang.Par.trafo", "dual"),
               interval=NULL, tol=NULL, ...)
dual_bound(s, d, pF, tol=.Machine$double.eps^0.25, ...)

rearrange(X, tol, tol.type=c("relative", "absolute"), maxiter=Inf,
          method=c("worst", "best"))
RA(alpha, d, qF, N, abstol=NULL, maxiter=Inf, method=c("worst", "best"),
   sample=TRUE)
ARA(alpha, d, qF, N.exp=seq(8, 20, by=1), reltol=c(0.001, 0.01), maxiter=12,
    method=c("worst", "best"), sample=TRUE)
}
\arguments{
  \item{alpha}{Value-at-Risk confidence level (e.g., 0.99).}
  \item{d}{Dimension (number of risk factors; \eqn{\ge 2}{>=2}).}
  \item{qF}{The marginal quantile function in the homogeneous case or,
    for \code{crude_VaR_bounds()}, \code{RA()} and \code{ARA()},
    a list (of length \code{d}) with the marginal quantile functions.}
  \item{method}{A \code{\link{character}} string. For
    \describe{
      \item{\code{VaR_bounds_hom()}:}{\code{method="Wang"},
	\code{method="Wang.Par"} and \code{method="Wang.Par.trafo"}
         apply the approaches of Embrechts et al. (2014,
	 Proposition 3.1) for computing best (i.e., smallest) and
	 worst (i.e., largest) Value-at-Risk. The
         latter two methods assume Pareto margins and thus do
	 not require numerical integration;
         \code{method="Wang.Par.trafo"} first transforms the
         root finding problem for computing worst Value-at-Risk to
	 a different scale.
	\code{method="dual"} applies the dual bound approach as
	 in Embrechts et al. (2013, Proposition 4)
         for computing worst Value-at-Risk (no value for the best
         Value-at-Risk can be obtained with this approach and thus
          \code{\link{NA}} is returned for the best Value-at-Risk).}
      \item{\code{rearrange()}, \code{RA()}, \code{ARA()}:}{
	\code{method} indicates whether bounds
	for the best or for the worst Value-at-Risk should be computed.
	These bounds are termed
          \eqn{\underline{s}_N} and \eqn{\overline{s}_N} in the
          literature (and below) and are theoretically not guaranteed
          bounds of (best or worst) Value-at-Risk
	  (however treated as such in practice and typically in line
	  with results from \code{VaR_bounds_hom()} in the homogeneous
	  case, for example).}
  }}
  \item{interval}{Initial interval (a \code{\link{numeric}(2)}) for
    computing worst Value-at-Risk. If not provided, these are the defaults chosen:
    \describe{
      \item{\code{method="Wang"}:}{The initial
	interval is \eqn{[0,(1-\alpha)/d]}{[0,(1-alpha)/d]}.}
      \item{\code{method="Wang.Par"}:}{The initial
	interval is \eqn{[c_l,c_u]}, where \eqn{c_l} and \eqn{c_u}
	are chosen as in Hofert et al. (2015).}
      \item{\code{method="Wang.Par.trafo"}:}{The initial
	interval is \eqn{[x_l,x_u]},
	where \eqn{x_l} and \eqn{x_u} are chosen as in Hofert et al. (2015).}
      \item{\code{method="dual"}:}{In this case, no good defaults are known.
	Note that the lower endpoint of the initial interval has to be
	sufficiently large in order for the the inner root-finding algorithm
	to find a root; see Details.}
  }}
  \item{tol}{
    \describe{
      \item{\code{VaR_bounds_hom()}:}{The tolerance for the
	root-finding algorithm for computing worst Value-at-Risk as passed
	to \code{\link{uniroot}()}. This defaults to
	\eqn{2.2204*10^{-16}}{2.2204*10^{-16}}
	for \code{method="Wang"} or \code{method="Wang.Par"} (where a
	smaller tolerance is crucial) and to \code{\link{uniroot}()}'s
	default \code{.Machine$double.eps^0.25} otherwise. Note that for
	\code{method="dual"}, \code{tol} is used for both the outer
	and the inner root-finding procedure.}
      \item{\code{rearrange()}:}{(Absolute or relative) tolerance to
	determine the (individual) convergence; if \code{NULL},
	the iteration is done until the matrix does not change anymore.}
  }}
  \item{tol.type}{\code{\link{character}} string indicating the
    type of convergence tolerance function to be used (\code{"relative"}
    for relative tolerance and \code{"absolute"} for absolute tolerance).}
  \item{s}{Dual bound evaluation point.}
  \item{pF}{The marginal loss distribution function.}
  \item{X}{An (\code{N}, \code{d})-matrix containing
    \eqn{N}-discretizations of the \eqn{d} quantiles to be rearranged.}
  \item{maxiter}{Maximal number of iterations over all columns of the
    underlying matrix of quantiles (can be set to \code{Inf}).}
  \item{N}{The number of discretization points.}
  \item{N.exp}{The exponents of the number of discretization points (a vector);
    the iteration then takes place over all discretization points
    \code{2^N.exp}. For each such number, at most \code{maxiter}-many
    iterations through all columns of the underlying matrix of quantiles
    are conducted, checking after each such iteration whether the
    desired relative tolerances	(both individually and jointly; see below)
    have been reached.}
  \item{abstol}{Absolute convergence tolerance \eqn{\epsilon}{epsilon}
    to determine the individual convergence, i.e., the change in the computed minimal
    (for worst Value-at-Risk) or maximal (for best Value-at-Risk)
    row sums for the lower bound \eqn{\underline{s}_N} and the upper
    bound \eqn{\overline{s}_N}. If \code{NULL} (the
    default) and \code{maxiter} has not been reached, the iterations
    through the columns of the underlying matrix of quantiles is done
    until the matrix does not change anymore, i.e., until each column
    is oppositely ordered to the sum of all others.}
  \item{reltol}{A vector containing two relative convergence tolerances,
    the individual relative tolerance (first component; used to determine
    convergence of the minimal (for worst Value-at-Risk) or maximal
    (for best Value-at-Risk) row sums for
    \eqn{\underline{s}_N} and \eqn{\overline{s}_N})
    and the joint relative tolerance (second component; relative tolerance
    between the computed \eqn{\underline{s}_N} and
    \eqn{\overline{s}_N} with respect to \eqn{\overline{s}_N}).
    Note that \code{reltol[1]} is allowed to be
    \code{NULL}; see \code{abstol} above for how to interpret this case.}
  \item{sample}{\code{\link{logical}} indicating whether each column of
    the two internal matrices of quantiles (see Step 3 of the Rearrangement
    Algorithm in Embrechts et al. (2013))
    are sampled before the iteration begins.}
  \item{\dots}{
    \describe{
      \item{\code{VaR_bounds_hom()}:}{For \code{d=2} \dots must contain the quantile
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

  \code{VaR_bounds_hom()} returns the best and worst Value-at-Risk at
  confidence level \eqn{\alpha}{alpha} for \eqn{d} risks with equal
  distribution function specified by the ellipsis \code{...}.

  \code{dual_bound()} returns the value of the dual bound \eqn{D(s)} as
  given in Embrechts, Puccetti, \enc{Rüschendorf}{Rueschendorf}
  (2013, Eq. (12)).

  \code{rearrange()} returns a \code{\link{list}} containing
  \describe{
    \item{\code{bound}:}{The computed \eqn{\underline{s}_N}
      or \eqn{\overline{s}_N}.}
    \item{\code{tol}:}{The tolerance (i.e., the (absolute or
      relative) change of the minimal (for worst Value-at-Risk) or
      maximal (for best Value-at-Risk) row sum in the last iteration).}
    \item{\code{X.rearranged}:}{A \code{link{matrix}} containing the
      rearranged input matrix \code{X}.}
    \item{\code{row.sums}:}{A \code{\link{matrix}} with the same number of rows as
      \code{X}, giving the row sums after each iteration of
      rearrangements of \code{X}; one iteration corresponds to looping
      over all columns and oppositely reorder each of them with respect
      to the sum of all others. The number of columns of \code{row.sums}
      thus corresponds to \code{num.iter}.}
  }

  \code{RA()} returns a \code{\link{list}} containing
  \describe{
    \item{\code{bounds}:}{A bivariate vector containing the computed
      \eqn{\underline{s}_N} and \eqn{\overline{s}_N} (the so-called
      rearrangement range) which are typically treated as bounds for
      (worst or best) Value-at-Risk; see also above.}
    \item{\code{rel.ra.gap}:}{relative tolerance (also known as
      relative rearrangement gap) between
      \eqn{\underline{s}_N} and \eqn{\overline{s}_N} computed with
      respect to \eqn{\overline{s}_N}.}
    \item{\code{individual.tol}:}{A bivariate vector containing the
      individual tolerances (i.e., the absolute change of the minimal
      (for worst Value-at-Risk) or maximal (for best Value-at-Risk)
      row sums for computing \eqn{\underline{s}_N} and \eqn{\overline{s}_N}
      in the last iteration of \code{rearrange()}); see also
      \code{tol} returned by \code{rearrange()} above.}
    \item{\code{X}:}{The constructed (\code{N}, \code{d})-matrices
      for computing \eqn{\underline{s}_N} and \eqn{\overline{s}_N}.}
    \item{\code{X.rearranged}:}{A \code{\link{list}} containing the rearranged
      matrices \code{X} for the computed \eqn{\underline{s}_N} and
      \eqn{\overline{s}_N}; see also \code{X.rearranged} returned by
      \code{rearrange()} above.}
    \item{\code{num.iter}:}{A bivariate vector containing (for
      \eqn{\underline{s}_N} and \eqn{\overline{s}_N}), the number
      of iterations over all columns of the underlying matrix \code{X}
      of quantiles; see also \code{num.iter}.}
    \item{\code{row.sums}:}{A list of length two containing the computed
      matrices of row sums (one column for each iteration over all
      columns of \code{X}) for computing \eqn{\underline{s}_N} and
      \eqn{\overline{s}_N}; see also \code{row.sums} returned by
      \code{rearrange()} above.}
    \item{\code{m.row.sums}:}{A list of length two containing minimal
      (for the worst Value-at-Risk) or maximal (for the best Value-at-Risk)
      row sums for computing \eqn{\underline{s}_N} and \eqn{\overline{s}_N}.}
  }

  \code{ARA()} returns a \code{\link{list}} containing
  \describe{
    \item{\code{bounds}:}{See \code{RA()}.}
    \item{\code{rel.ra.gap}:}{See \code{RA()}.}
    \item{\code{individual.tol}:}{See \code{RA()}.}
    \item{\code{joint.tol}:}{The joint relative tolerance reached (i.e.,
      the relative tolerance between
      \eqn{\underline{s}_N} and \eqn{\overline{s}_N}
      computed with respect to the upper bound).}
    \item{\code{X}:}{See \code{RA()}.}
    \item{\code{X.rearranged}:}{See \code{RA()}.}
    \item{\code{N.used}:}{The actual \code{N} used for computing
      the (final) \eqn{\underline{s}_N} and \eqn{\overline{s}_N}.}
    \item{\code{num.iter}:}{See \code{RA()}.}
    \item{\code{row.sums}:}{See \code{RA()}.}
    \item{\code{m.row.sums}:}{See \code{RA()}.}
  }
}
\details{
  For \code{d=2}, \code{VaR_bounds_hom()} uses the method of
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
  be positive (but must lie below the unknown root). Note that the upper
  endpoint \eqn{(1-\alpha)/d}{(1-alpha)/d} also happens to be a
  root and thus one needs a proper initional interval containing the
  root and being stricticly contained in
  \eqn{(0,(1-\alpha)/d}{(1-alpha)/d)}.
  In the case of Pareto margins, Hofert et al. (2015) have
  derived such an initial (which is used by
  \code{method="Wang.Par"} and \code{method="Wang.Par.trafo"};
  see Hofert et al. (2015) for more details.
  Also note that the chosen smaller default tolerances for
  \code{\link{uniroot}()} in case of \code{method="Wang"} and
  \code{method="Wang.Par"} are crucial for obtaining reliable
  Value-at-Risk values; see Hofert et al. (2015).

  For \code{method="dual"} for computing worst Value-at-Risk, the method
  presented of Embrechts et al. (2013, Proposition 4) is implemented.
  This requires two (nested) \code{\link{uniroot}()}, and an
  \code{\link{integrate}()}. For the inner root-finding procedure to
  find a root, the lower endpoint of the provided initial
  \code{interval} has to be \dQuote{sufficiently large}.

  Note that these approaches for computing the
  Value-at-Risk bounds in the homogeneous case are numerically non-trivial;
  see the source code and \code{demo(VaR_bounds)} for more details. As a
  rule of thumb, use \code{method="Wang"} if you have to (i.e., if the
  margins are not Pareto) and \code{method="Wang.Par"} if you can (i.e.,
  if the margins are Pareto); in the former case, it might be wise to
  follow along the lines of \code{method="Wang.Par.trafo"} to implement
  one's own method (e.g., if our \code{method="Wang"} turns out to be
  numerically critical), see Hofert et al. (2015) for more details.
  Furthermore, it is not recommended to use (the numerically more challenging)
  \code{method="dual"}.


  Concerning the inhomogeneous case, \code{rearrange()} is an auxiliary
  function which is called by \code{RA()} and \code{ARA()}. For
  performance reasons, no checking is done. It should only be used by
  experts.

  For the Rearrangement Algorithm \code{RA()}, unless
  \code{maxiter}-many iterations have been reached, convergence of
  \eqn{\underline{s}_N} and \eqn{\overline{s}_N} is determined if the
  minimal (for the worst Value-at-Risk) or maximal (for the best
  Value-at-Risk) row sums after an iteration did not change by more than
  (the absolute tolerance) \code{abstol} (so \eqn{\le\epsilon}{<=
  eps}). This is slightly different from Embrechts et al. (2013) who use
  \eqn{<\epsilon}{< eps} but this slight modification has the advantage
  of allowing \eqn{\epsilon=0}{eps=0}.

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
  \code{demo(VaR_bounds)} and Hofert et al. (2015) for the reasons),
  it is recommended to use \code{ARA()} instead of \code{RA()}.

  On the theoretical side, let us again stress the following.
  \code{rearrange()}, \code{RA()} and \code{ARA()} compute
  \eqn{\underline{s}_N} and \eqn{\overline{s}_N} which are, from a
  practical point of view, treated as bounds for the worst
  (i.e., largest) or the best (i.e., smallest) Value-at-Risk
  (whatever is chosen with \code{method}), but which are not
  known to be such bounds from a theoretical point of view; see also
  above. Calling them \dQuote{bounds} for worst or best Value-at-Risk is
  thus theoretically not correct (unless proven) but \dQuote{practical}.
  The literature thus speaks of \eqn{(\underline{s}_N, \overline{s}_N)}
  as the rearrangement range (rather than an interval containing worst
  or best Value-at-Risk).
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
  \code{demo(VaR_bounds)} for more example calls, numerical challenges
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


## d=2: Compute best/worst VaR explicitly (hom. case) and compare with (A)RA ###

d <- 2 # dimension

## Explicit
VaRbounds <- VaR_bounds_hom(alpha, d=d, qF=qF) # (best VaR, worst VaR)

## Adaptive Rearrangement Algorithm (ARA)
set.seed(271) # set seed (for reproducibility)
ARAbest  <- ARA(alpha, d=d, qF=qF, method="best")
ARAworst <- ARA(alpha, d=d, qF=qF)

## Rearrangement Algorithm (RA) with N as in ARA and abstol (roughly) chosen as in ARA
## Best VaR
RAbest  <- RA(alpha, d=d, qF=qF, N=ARAbest$N.used,
              abstol=mean(tail(abs(diff(ARAbest$m.row.sums$low)), n=1),
                          tail(abs(diff(ARAbest$m.row.sums$up)), n=1)),
              method="best")
## Worst VaR
RAworst <- RA(alpha, d=d, qF=qF, N=ARAworst$N.used,
              abstol=mean(tail(abs(diff(ARAworst$m.row.sums$low)), n=1),
                          tail(abs(diff(ARAworst$m.row.sums$up)), n=1)))

## Compare (comparably large tolerances required to contain both (A)RA bounds)
## Best VaR
stopifnot(all.equal(c(ARAbest$bounds[1], ARAbest$bounds[2],
                      RAbest$bounds[1],  RAbest$bounds[2]),
                    rep(VaRbounds[1], 4), tolerance=0.004))
## Worst VaR
stopifnot(all.equal(c(ARAworst$bounds[1], ARAworst$bounds[2],
                      RAworst$bounds[1],  RAworst$bounds[2]),
                    rep(VaRbounds[2], 4), tolerance=0.003))


## d=8: Compute best/worst VaR (hom. case) and compare with (A)RA ##############

d <- 8 # dimension

## Compute VaR bounds with various methods
I <- crude_VaR_bounds(alpha, d=d, qF=qF) # crude bound
VaR.Wang           <- VaR_bounds_hom(alpha, d=d, method="Wang", qF=qF)
VaR.Wang.Par       <- VaR_bounds_hom(alpha, d=d, method="Wang.Par", theta=th)
VaR.Wang.Par.trafo <- VaR_bounds_hom(alpha, d=d, method="Wang.Par.trafo", theta=th)
VaR.dual           <- VaR_bounds_hom(alpha, d=d, method="dual", interval=I, pF=pF)

## Adaptive Rearrangement Algorithm (ARA) (with defaults)
set.seed(271) # set seed (for reproducibility)
ARAbest  <- ARA(alpha, d=d, qF=qF, method="best")
ARAworst <- ARA(alpha, d=d, qF=qF)

## Rearrangement Algorithm (RA) with N as in ARA and abstol (roughly) chosen as in ARA
## Best VaR
RAbest  <- RA(alpha, d=d, qF=qF, N=ARAbest$N.used,
              abstol=mean(tail(abs(diff(ARAbest$m.row.sums$low)), n=1),
                          tail(abs(diff(ARAbest$m.row.sums$up)), n=1)),
              method="best")
## Worst VaR
RAworst <- RA(alpha, d=d, qF=qF, N=ARAworst$N.used,
              abstol=mean(tail(abs(diff(ARAworst$m.row.sums$low)), n=1),
                          tail(abs(diff(ARAworst$m.row.sums$up)), n=1)))

## Compare (comparably large tolerances required to contain both (A)RA bounds)
## Best VaR
stopifnot(all.equal(c(VaR.Wang[1], VaR.Wang.Par.trafo[1],
                      ARAbest$bounds[1], ARAbest$bounds[2],
                      RAbest$bounds[1],  RAbest$bounds[2]),
                    rep(VaR.Wang.Par[1], 6), tolerance=0.004))
## Worst VaR
stopifnot(all.equal(c(VaR.Wang[2], VaR.Wang.Par.trafo[2], VaR.dual[2],
                      ARAworst$bounds[1], ARAworst$bounds[2],
                      RAworst$bounds[1],  RAworst$bounds[2]),
                    rep(VaR.Wang.Par[2], 7), tolerance=0.003))


## Using (some of) the additional results computed by (A)RA() ##################

xlim <- c(1, max(sapply(RAworst$m.row.sums, length)))
ylim <- range(RAworst$m.row.sums)
plot(RAworst$m.row.sums[[2]], type="b", xlim=xlim, ylim=ylim, xlab="Iteration",
     ylab=paste0("Minimal row sum"),
     main=substitute("RA for worst VaR with"~alpha==a.*","~d==d.*" and Par("*
                     th.*")", list(a.=alpha, d.=d, th.=th)))
lines(1:length(RAworst$m.row.sums[[1]]), RAworst$m.row.sums[[1]], type="b",
      col="blue")
legend("bottomright", bty="n", lty=rep(1,2),
       col=c("black", "blue"), legend=c("upper bound", "lower bound"))
## => One should use ARA() instead of RA()


## "Reproducing" examples from Embrechts et al. (2013) #########################

## "Reproducing" Table 1 (but seed and eps unknown)
set.seed(271)
RAworst <- RA(alpha, d=3, qF=qF, N=50)
stopifnot(all.equal(RAworst$bounds[1], 44.84026, tolerance=1e-6),
          all.equal(RAworst$bounds[2], 46.40413, tolerance=1e-6))

## "Reproducing" Table 3 for alpha=0.99 (but seed unknown)
N <- 1e5 # we use a smaller N here to save run time
eps <- 0.1 # absolute tolerance
xi <- c(1.19, 1.17, 1.01, 1.39, 1.23, 1.22, 0.85, 0.98)
beta <- c(774, 254, 233, 412, 107, 243, 314, 124)
qF.lst <- lapply(1:8, function(j) { function(p) qGPD(p, xi=xi[j], beta=beta[j]) })
set.seed(271)
RAbest  <- RA(0.99, qF=qF.lst, N=N, abstol=eps, method="best")
RAworst <- RA(0.99, qF=qF.lst, N=N, abstol=eps)
stopifnot(all.equal(RAbest$bounds[1], 1.780585e5,  tolerance=1e-7),
          all.equal(RAbest$bounds[2], 1.78304e5,   tolerance=1e-7),
          all.equal(RAworst$bounds[1], 2.556191e6, tolerance=1e-7),
          all.equal(RAworst$bounds[2], 2.556441e6, tolerance=1e-7))
}
\keyword{distribution}