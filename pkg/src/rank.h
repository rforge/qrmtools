/* C header for rank.c ********************************************************/

#ifndef rank_H
#define rank_H

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

/* For faster rank() in rearrange() */
SEXP rank_aux_(SEXP x);

#endif

