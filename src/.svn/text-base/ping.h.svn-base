#include <Rinternals.h>
#include <Rdefines.h>
#include <Rmath.h>

void background(double *dataF, double *dataR, int *nReadsF, int *nReadsR, int *mapS, int *mapE, int *lMap, int *gapS, int *gapE, int *lGap, double *pRetain);
SEXP fitPING(SEXP segReadsList, SEXP paraEM, SEXP paraPrior, SEXP minReads, SEXP detailS, SEXP rescaleS, SEXP calphaS, SEXP PES);


SEXP getVector(SEXP list, SEXP ind);
SEXP getK(SEXP list);
SEXP getScore(SEXP list);
SEXP getScoreF(SEXP list);
SEXP getScoreR(SEXP list);
SEXP getSegL(SEXP list);

SEXP getMax(SEXP list);
SEXP getMin(SEXP list);

SEXP getChr(SEXP list);
SEXP getMap(SEXP list);
void wThreCounts(int *step, int *dataF, int *dataR, int *nReadsF, int *nReadsR, int *width, int *scoreF, int *scoreR);
void callRegions(int *center, int *lengthCenter, int *width, int *scoreF, int *scoreR, int *scoreRegionF, int *scoreRegionR, int *cutoff, int *StartRegion, int *EndRegion, int *nRegions);
void callRegionsL(int *center, int *nProbes, int *width, int *scoreF, int *scoreR, int *scoreRegionF, int *scoreRegionR, int *cutoff, int *StartRegion, int *EndRegion, int *nRegions, int maxStep, int kStep, int minL);
SEXP segReads(SEXP chr, SEXP dataF, SEXP dataR, SEXP contF, SEXP contR, SEXP StartMap, SEXP EndMap, SEXP jitter, SEXP width, SEXP cutoff, SEXP step, SEXP maxStep, SEXP minLength);
SEXP segR(SEXP chr, SEXP dataF, SEXP dataR, SEXP contF, SEXP contR, SEXP StartRegion, SEXP EndRegion, SEXP StartMap, SEXP EndMap, SEXP jitter, int nRegions);
SEXP segR2(SEXP chr, SEXP dataF, SEXP dataR, SEXP contF, SEXP contR, SEXP StartRegion, SEXP EndRegion, SEXP StartMap, SEXP EndMap, SEXP jitter, int nRegions);
SEXP segR3(SEXP chr, SEXP dataF, SEXP dataR, SEXP contF, SEXP contR, SEXP StartRegion, SEXP EndRegion, SEXP StartMap, SEXP EndMap, SEXP jitter, int nRegions);
SEXP segReadsAll(SEXP data, SEXP dataC, SEXP StartMap, SEXP EndMap, SEXP jitter, SEXP paraSW, SEXP maxStep, SEXP minLength);

SEXP getDensity(SEXP ping, SEXP strand, SEXP step, SEXP filter, SEXP sum, SEXP scale);
SEXP getDensityList(SEXP pingList, SEXP strand, SEXP step, SEXP filter, SEXP sum, SEXP scale);
SEXP getListElement(SEXP list, const char *str);
