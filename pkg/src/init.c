/* Register routines with R ***************************************************/

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include "rank.h"


static const R_CallMethodDef callMethods[] = {
	{"rank_aux_", (DL_FUNC) &rank_aux_, 5},
	{NULL, NULL, 0}
};

void R_init_qrng(DllInfo *dll)
{
    R_useDynamicSymbols(dll, FALSE);
    R_registerRoutines(dll, NULL, callMethods, NULL, NULL); /* s. WRE (2015, Section 5.4) */
}
