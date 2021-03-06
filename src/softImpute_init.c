// Automatically generated by SUtools, edited by TH because pcol was integer
#ifndef R_SOFTIMPUTE_H
#define R_SOFTIMPUTE_H
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("softImpute", String)
#else
#define _(String) (String)
#endif


#define FDEF(name)  {#name, (DL_FUNC) &F77_SUB(name), sizeof(name ## _t)/sizeof(name ## _t[0]), name ##_t}
void F77_SUB(plusregc)(
int *nrow,
int *ncol,
int *nrank,
double *u,
double *v,
int *irow,
int *pcol,
int *nomega,
double *x,
double *zy,
double *zz
);
 
static R_NativePrimitiveArgType plusregc_t[] = {
INTSXP,
INTSXP,
INTSXP,
REALSXP,
REALSXP,
INTSXP,
INTSXP,
INTSXP,
REALSXP,
REALSXP,
REALSXP
};
void F77_SUB(suv)(
int *nrow,
int *ncol,
int *nrank,
double *u,
double *v,
int *irow,
int *jcol,
int *nomega,
double *r
);
 
static R_NativePrimitiveArgType suv_t[] = {
INTSXP,
INTSXP,
INTSXP,
REALSXP,
REALSXP,
INTSXP,
INTSXP,
INTSXP,
REALSXP
};
void F77_SUB(suvc)(
int *nrow,
int *ncol,
int *nrank,
double *u,
double *v,
int *irow,
int *pcol,
int *nomega,
double *r
);
 
static R_NativePrimitiveArgType suvc_t[] = {
INTSXP,
INTSXP,
INTSXP,
REALSXP,
REALSXP,
INTSXP,
INTSXP,
INTSXP,
REALSXP
};

static R_FortranMethodDef fMethods[] = {
FDEF(plusregc) ,
FDEF(suv) ,
FDEF(suvc) ,
{NULL, NULL, 0}
};

void R_init_softImpute(DllInfo *dll){
  R_registerRoutines(dll, NULL, NULL, fMethods, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
#endif
