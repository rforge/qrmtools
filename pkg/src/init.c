/* Register routines with R ***************************************************/

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include "VaR_bounds.h"


static const R_CallMethodDef callMethods[] = {
	{"rank", (DL_FUNC) &rank_, 1},
	{"colsplit", (DL_FUNC) &colsplit_, 1},
	{NULL, NULL, 0}
};

void R_init_qrmtools(DllInfo *dll)
{
    R_useDynamicSymbols(dll, FALSE);
    R_registerRoutines(dll, NULL, callMethods, NULL, NULL); /* s. WRE (2015, Section 5.4) */
}
