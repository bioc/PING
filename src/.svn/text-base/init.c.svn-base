#include <R.h>
#include <R_ext/Rdynload.h>
#include "ping.h"


static const R_CallMethodDef CallEntries[] = {
    {"fitPING", (DL_FUNC)&fitPING, 8},
    {"segReadsAll", (DL_FUNC)&segReadsAll, 8},
    {"getSegL", (DL_FUNC)&getSegL, 1},
    {"getScore", (DL_FUNC)&getScore, 1},
    {"getMin", (DL_FUNC)&getMin, 1},
    {"getMax", (DL_FUNC)&getMax, 1},
    {"getScoreR", (DL_FUNC)&getScoreR, 1},
    {"getScoreF", (DL_FUNC)&getScoreF, 1},
    {"getChr", (DL_FUNC)&getChr, 1},
    {"getMap", (DL_FUNC)&getMap, 1},
    {"getK", (DL_FUNC)&getK, 1},
    {"getVector", (DL_FUNC)&getVector, 2},
    {"getDensity", (DL_FUNC)&getDensity, 6},
    {"getDensityList", (DL_FUNC)&getDensityList, 6},

    {NULL, NULL, 0}
};

static const R_CMethodDef CEntries[] = {
    {"background", (DL_FUNC) &background, 11},
    {NULL, NULL, 0}
};

void R_init_PING(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);

}





