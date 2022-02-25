#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>


/* .Fortran calls */
  extern void F77_NAME(pfishnet)(void *, void *, void *);
extern void F77_NAME(plognet)(void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
  {"pfishnet", (DL_FUNC) &F77_NAME(pfishnet), 29},
  {"plognet", (DL_FUNC) &F77_NAME(plognet), 31},
  {NULL, NULL, 0}
};

void R_init_rpair(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
