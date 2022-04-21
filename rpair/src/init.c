#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>


/* .Fortran calls */
extern void F77_NAME(pfishnet)(void *, void *, void *);
extern void F77_NAME(plognet)(void *, void *, void *);
extern void F77_NAME(psqhnet)(void *, void *, void *);
extern void F77_NAME(phuhnet)(void *, void *, void *);
extern void F77_NAME(get_int_parms)(void *, void *, void *);
extern void F77_NAME(chg_fract_dev)(void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
  {"pfishnet", (DL_FUNC) &F77_NAME(pfishnet), 29},
  {"plognet", (DL_FUNC) &F77_NAME(plognet), 31},
  {"psqhnet", (DL_FUNC) &F77_NAME(psqhnet), 26},
  {"phuhnet", (DL_FUNC) &F77_NAME(phuhnet), 27},
  {"get_int_parms", (DL_FUNC) &F77_NAME(get_int_parms), 7},
  {"chg_fract_dev", (DL_FUNC) &F77_NAME(chg_fract_dev), 1},
  {NULL, NULL, 0}
};

void R_init_rpair(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
