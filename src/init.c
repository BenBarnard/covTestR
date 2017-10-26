#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _covTestR_Ahmad2017Stat(SEXP);
extern SEXP _covTestR_bilinearcube(SEXP);
extern SEXP _covTestR_bilinearoff(SEXP);
extern SEXP _covTestR_bilinearquad(SEXP);
extern SEXP _covTestR_bilinearsquare(SEXP);
extern SEXP _covTestR_BoxesMStat(SEXP);
extern SEXP _covTestR_c1(SEXP);
extern SEXP _covTestR_c2(SEXP);
extern SEXP _covTestR_c3(SEXP);
extern SEXP _covTestR_Chaipitak2013poolStat(SEXP);
extern SEXP _covTestR_Chaipitak2013Stat(SEXP);
extern SEXP _covTestR_Ishii2016Stat(SEXP);
extern SEXP _covTestR_quadra(SEXP);
extern SEXP _covTestR_Schott2001Stat(SEXP);
extern SEXP _covTestR_Schott2007pooledStat(SEXP);
extern SEXP _covTestR_Schott2007Stat(SEXP);
extern SEXP _covTestR_Srivastava2007Stat(SEXP);
extern SEXP _covTestR_Srivastava2014poolStat(SEXP);
extern SEXP _covTestR_Srivastava2014Stat(SEXP);
extern SEXP _covTestR_SrivastavaYanagihara2010Stat(SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_covTestR_Ahmad2017Stat",                (DL_FUNC) &_covTestR_Ahmad2017Stat,                1},
    {"_covTestR_bilinearcube",                 (DL_FUNC) &_covTestR_bilinearcube,                 1},
    {"_covTestR_bilinearoff",                  (DL_FUNC) &_covTestR_bilinearoff,                  1},
    {"_covTestR_bilinearquad",                 (DL_FUNC) &_covTestR_bilinearquad,                 1},
    {"_covTestR_bilinearsquare",               (DL_FUNC) &_covTestR_bilinearsquare,               1},
    {"_covTestR_BoxesMStat",                   (DL_FUNC) &_covTestR_BoxesMStat,                   1},
    {"_covTestR_c1",                           (DL_FUNC) &_covTestR_c1,                           1},
    {"_covTestR_c2",                           (DL_FUNC) &_covTestR_c2,                           1},
    {"_covTestR_c3",                           (DL_FUNC) &_covTestR_c3,                           1},
    {"_covTestR_Chaipitak2013poolStat",        (DL_FUNC) &_covTestR_Chaipitak2013poolStat,        1},
    {"_covTestR_Chaipitak2013Stat",            (DL_FUNC) &_covTestR_Chaipitak2013Stat,            1},
    {"_covTestR_Ishii2016Stat",                (DL_FUNC) &_covTestR_Ishii2016Stat,                1},
    {"_covTestR_quadra",                       (DL_FUNC) &_covTestR_quadra,                       1},
    {"_covTestR_Schott2001Stat",               (DL_FUNC) &_covTestR_Schott2001Stat,               1},
    {"_covTestR_Schott2007pooledStat",         (DL_FUNC) &_covTestR_Schott2007pooledStat,         1},
    {"_covTestR_Schott2007Stat",               (DL_FUNC) &_covTestR_Schott2007Stat,               1},
    {"_covTestR_Srivastava2007Stat",           (DL_FUNC) &_covTestR_Srivastava2007Stat,           1},
    {"_covTestR_Srivastava2014poolStat",       (DL_FUNC) &_covTestR_Srivastava2014poolStat,       1},
    {"_covTestR_Srivastava2014Stat",           (DL_FUNC) &_covTestR_Srivastava2014Stat,           1},
    {"_covTestR_SrivastavaYanagihara2010Stat", (DL_FUNC) &_covTestR_SrivastavaYanagihara2010Stat, 1},
    {NULL, NULL, 0}
};

void R_init_covTestR(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
