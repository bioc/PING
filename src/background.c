#include <R.h>
#include <Rmath.h>
#include <float.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

double ran_flat(double m, double M);
void background(double *dataF, double *dataR, int *nReadsF, int *nReadsR, int *mapS, int *mapE, int *lMap, int *gapS, int *gapE, int *lGap, double *pRetain)
{
  int FLAG=1,i=0,j=0,l=0,countR=0;
  double data=0,mF=dataF[0],MF=dataF[*nReadsF-1],mR=dataR[0],MR=dataR[*nReadsR-1];
  double length=0;
  /* Random number parameters used for the imputation stuff */
  // const gsl_rng_type * T;
  // gsl_rng * r;
  // 
  // gsl_rng_env_setup();
  // 
  // /* Allocate memory for the random number stuff */
  // T = gsl_rng_default;
  // r = gsl_rng_alloc(T);
	
  GetRNGstate();
	
  for(i=0;i<*nReadsF;i++)
  {
    FLAG=0;
    while(FLAG==0)
    {
      FLAG=1;
      data=(int)ran_flat(mF, MF);

      for(l=0;l<*lMap;l++)
      {
        /* It is in a non mappable region */
        if((data>=mapS[l]) && (data<=mapE[l]))
        {
          FLAG=0;
        }
      }
      /* Even if it is in a non mappable region we can still retain the read with pro pRetain */
      if((FLAG==0) && (ran_flat(0, 1)<*pRetain))
      {
        FLAG==1;
      }
      for(l=0;l<*lGap;l++)
      {
        /* It is in a gap region */
        if(data>=gapS[l] && data<=gapE[l])
        {
          FLAG=0;
        }
      }
    }
    dataF[i]=data;
  }
    
    // length=175+gsl_ran_gaussian(r,30);
    // /* End the fragment length */
    // data+=length;
    // FLAG=1;
    // for(l=0;l<*lMap;l++)
    // {
    //   /* It is a non mappable regions */
    //   if(data>=mapS[l] && data<=mapE[l])
    //   {
    //     FLAG=0;
    //   }
    // }
    // /* Check that the end of the fragment is within a mappable region, otherwise we don't retain it */
    // if(FLAG==1)
    // {
    //   dataR[countR]=data;
    //   countR++;
    // }
  
  // *nReadsR=countR;

  for(j=0;j<*nReadsR;j++)
  {
    FLAG=0;
    while(FLAG==0)
    {
      FLAG=1;
      data=ran_flat(mR, MR);
      for(l=0;l<*lMap;l++)
      {
        /* It is in a non mappable region */
        if((data>=mapS[l]) && (data<=mapE[l]))
        {
          FLAG=0;
        }
      }
      /* Even if it is in a non mappable region we can still retain the read with pro pRetain */
      if((FLAG==0) && (ran_flat(0, 1)<*pRetain))
      {
        FLAG==1;
      }
      for(l=0;l<*lGap;l++)
      {
        /* It is in a gap region */
        if((data>=gapS[l]) && (data<=gapE[l]))
        {
          FLAG=0;
        }
      }
    }
    dataR[j]=data;      
  }
  // gsl_rng_free(r);
  PutRNGstate();
}

double ran_flat(double m, double M)
{
  double u=unif_rand();
  return((M-m)*u);
}
