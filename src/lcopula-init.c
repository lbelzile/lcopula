#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>


/* .Call calls */
extern SEXP _lcopula_ctspecdens(SEXP, SEXP, SEXP);
extern SEXP _lcopula_dirspecdens(SEXP, SEXP, SEXP, SEXP);
extern SEXP _lcopula_marginCombo(SEXP, SEXP);
extern SEXP _lcopula_negdirspecdens(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_lcopula_ctspecdens",     (DL_FUNC) &_lcopula_ctspecdens,     3},
    {"_lcopula_dirspecdens",    (DL_FUNC) &_lcopula_dirspecdens,    4},
    {"_lcopula_marginCombo",    (DL_FUNC) &_lcopula_marginCombo,    2},
    {"_lcopula_negdirspecdens", (DL_FUNC) &_lcopula_negdirspecdens, 4},
    {NULL, NULL, 0}
};

void R_init_lcopula(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
