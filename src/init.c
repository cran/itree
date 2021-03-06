/*
 * Same as init.c in rpart, but the last function
 * is called R_init_itree instead of R_init_rpart.
 */
#include <R.h>
#include <Rinternals.h>
#include "R_ext/Rdynload.h"
#include "rpartS.h"
#include "node.h"
#include "rpartproto.h"

SEXP init_rpcallback(SEXP rhox, SEXP ny, SEXP nr, SEXP expr1x, SEXP expr2x);
void rpartexp2(Sint *n2, double *y, double *eps, long *keep);


static const R_CMethodDef CEntries[] = {
    {"pred_rpart", (DL_FUNC) &pred_rpart, 13},
    {"rpartexp2", (DL_FUNC) &rpartexp2, 4},
    {"s_to_rp", (DL_FUNC) &s_to_rp, 16},
    {"s_to_rp2", (DL_FUNC) &s_to_rp2, 14},
    {"s_xpred", (DL_FUNC) &s_xpred, 19},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"init_rpcallback", (DL_FUNC) &init_rpcallback, 5},
    {NULL, NULL, 0}
};

//ALG: 6/13/2013
//new name, same function.
void R_init_itree(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
