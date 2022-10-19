
/* 
   Native symbol registration table for polyclip package

*/

#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

SEXP Csimplify(SEXP A,
	       SEXP pft,
	       SEXP X0,
	       SEXP Y0,
	       SEXP Eps);

SEXP Cclipbool(SEXP A,
	       SEXP B,
	       SEXP pftA,
	       SEXP pftB,
	       SEXP ct,
	       SEXP X0,
	       SEXP Y0,
	       SEXP Eps,
               SEXP clo); 

SEXP Cpolyoffset(SEXP A,
		 SEXP del,
		 SEXP jt,
		 SEXP mlim,
		 SEXP atol,
		 SEXP X0,
		 SEXP Y0,
		 SEXP Eps);

SEXP Clineoffset(SEXP A,
		 SEXP del,
		 SEXP jt,
		 SEXP et,
		 SEXP mlim,
		 SEXP atol,
		 SEXP X0,
		 SEXP Y0,
		 SEXP Eps);

SEXP Cminksum(SEXP A,   
	      SEXP B,   
	      SEXP clo, 
	      SEXP X0,
	      SEXP Y0,
	      SEXP Eps);

SEXP Cpiptest(SEXP P,
	      SEXP A,
	      SEXP X0,
	      SEXP Y0,
	      SEXP Eps); 

static const R_CMethodDef CEntries[] = {
    {NULL, NULL, 0}
};
  
static const R_CallMethodDef CallEntries[] = {
    {"Csimplify",        (DL_FUNC) &Csimplify,         5},
    {"Cclipbool",        (DL_FUNC) &Cclipbool,         9},
    {"Cpolyoffset",      (DL_FUNC) &Cpolyoffset,       8},
    {"Clineoffset",      (DL_FUNC) &Clineoffset,       9},
    {"Cminksum",         (DL_FUNC) &Cminksum,          6},
    {"Cpiptest",         (DL_FUNC) &Cpiptest,          5},
    {NULL, NULL, 0}
};


void R_init_polyclip(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
