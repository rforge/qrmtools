/* C header for VaR_bounds.c **************************************************/

#ifndef VaR_bounds_H
#define VaR_bounds_H

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

/* For faster rank() in rearrange() */
SEXP rank_(SEXP x);

SEXP colsplit_(SEXP x);

#endif

