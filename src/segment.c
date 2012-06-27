#include <R.h>
#include <Rmath.h>
#include <float.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <Rdefines.h>
//GSL stuff
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

SEXP getDensity(SEXP ping, SEXP strand, SEXP step, SEXP filter, SEXP sum, SEXP scale);
SEXP getDensityList(SEXP pingList, SEXP strand, SEXP step, SEXP filter, SEXP sum, SEXP scale);
SEXP getListElement(SEXP list, const char *str);

SEXP getDensityList(SEXP pingList, SEXP strand, SEXP step, SEXP filter, SEXP sum, SEXP scale)
{
  int l=0,i=0;
  SEXP ping, List, ans, ansTmp;
  SEXP x, y, chr, chromosome,names;
  int nProtected=0,totalLength=0;
  double *range;

  PROTECT(List=GET_SLOT(pingList,install("List"))); nProtected++;
  PROTECT(ans=NEW_LIST(3)); nProtected++;

  for(l=0;l<length(List);l++)
  {
    ping=VECTOR_ELT(List,l);
    if(strcmp(CHAR(STRING_ELT(GET_CLASS(ping), 0)),"ping")==0)
    {
      range=REAL(GET_SLOT(ping,install("range")));
      totalLength+=(int)((range[1]-range[0])/REAL(step)[0]);
    }
  }

  PROTECT(x=NEW_NUMERIC(totalLength));nProtected++;
  PROTECT(y=NEW_NUMERIC(totalLength));nProtected++;
  PROTECT(chr=NEW_STRING(totalLength));nProtected++;
  
  totalLength=0;
  for(l=0;l<length(List);l++)
  { 
    ping=VECTOR_ELT(List,l);
    if(strcmp(CHAR(STRING_ELT(GET_CLASS(ping), 0)),"ping")==0)
    {
    
      chromosome=GET_SLOT(ping,install("chr"));
      range=REAL(GET_SLOT(ping,install("range")));

      PROTECT(ansTmp=getDensity(ping, strand, step, filter, sum, scale));nProtected++;
      for(i=0;i<(int)((range[1]-range[0])/REAL(step)[0]);i++)
      {
        REAL(x)[i+totalLength]=REAL(VECTOR_ELT(ansTmp,0))[i];
        REAL(y)[i+totalLength]=REAL(VECTOR_ELT(ansTmp,1))[i];
        STRING_PTR(chr)[i+totalLength]=STRING_PTR(chromosome)[0];
      }    
      totalLength+=(int)((range[1]-range[0])/REAL(step)[0]);
      UNPROTECT(1);nProtected--;
    }
  }
    // I have added the option to interrupt R
    // R_CheckUserInterrupt();
  SET_VECTOR_ELT(ans,0,x);
  SET_VECTOR_ELT(ans,1,y);
  SET_VECTOR_ELT(ans,2,chr);
  PROTECT(names = allocVector(STRSXP, 3)); nProtected++; 
  SET_STRING_ELT(names, 0, mkChar("x"));
  SET_STRING_ELT(names, 1, mkChar("density"));
  SET_STRING_ELT(names, 2, mkChar("chr"));
  setAttrib(ans, R_NamesSymbol, names);  
  UNPROTECT(nProtected);
  return(ans);
}

SEXP getDensity(SEXP ping, SEXP strand, SEXP step, SEXP filter, SEXP sum, SEXP scale)
{
  int k=0,K=1,i=0;
  double *w, *mu, *delta, *sF, *sR, *range, *se, sumW=0;
  double *muFilter, *deltaFilter, *sFFilter, *sRFilter, *seFilter, *scoreFilter, *score;
  SEXP filter1;
  gsl_vector *One;  
  gsl_matrix *Density;
  int length, nProtected=0;
  SEXP ans, x, y, chr, names;
  gsl_vector_view view;
  chr=GET_SLOT(ping,install("chr"));
  
  if(strcmp(CHAR(STRING_ELT(GET_CLASS(ping), 0)),"pingError")==0)
  {
    return(R_NilValue);
  }
  else
  {
    K=length(VECTOR_ELT(GET_SLOT(ping,install("estimates")), 0));
    w=REAL(VECTOR_ELT(GET_SLOT(ping,install("estimates")), 0));
    mu=REAL(VECTOR_ELT(GET_SLOT(ping,install("estimates")), 1));
    delta=REAL(VECTOR_ELT(GET_SLOT(ping,install("estimates")),2));
    sF=REAL(VECTOR_ELT(GET_SLOT(ping,install("estimates")), 3));
    sR=REAL(VECTOR_ELT(GET_SLOT(ping,install("estimates")), 4));
    score=REAL(GET_SLOT(ping,install("score")));
    se=REAL(VECTOR_ELT(GET_SLOT(ping,install("estimates")), 5));

    deltaFilter=REAL(getListElement(filter, "delta"));
    sFFilter=REAL(getListElement(filter, "sigmaSqF"));
    sRFilter=REAL(getListElement(filter, "sigmaSqR"));
    seFilter=REAL(getListElement(filter, "se"));
    scoreFilter=REAL(getListElement(filter, "score"));

    range=REAL(GET_SLOT(ping,install("range")));

    PROTECT(ans=NEW_LIST(2));nProtected++;
    
    // Compute the number of steps I want
    length=(int)((range[1]-range[0])/REAL(step)[0]);
    // Allocate the memory for x and y
    PROTECT(x=NEW_NUMERIC(length));nProtected++;
    PROTECT(y=NEW_NUMERIC(length));nProtected++;
    Density=gsl_matrix_calloc(length,K);
    One=gsl_vector_alloc(K);
    gsl_vector_set_all(One, 1.0);

    for(i=0;i<length;i++)
    {
      //Grid
      REAL(x)[i]=range[0]+i*REAL(step)[0];

      for(k=0;k<K;k++)
      {
      //Check if this is a valid binding event
        if((delta[k]>deltaFilter[0] & delta[k]<deltaFilter[1]) & 
          (sF[k]>sFFilter[0] & sF[k]<sFFilter[1]) & 
          (sR[k]>sRFilter[0] & sR[k]<sRFilter[1]) & 
          (se[k]>seFilter[0] & se[k]<seFilter[1]) &
          (score[k]>scoreFilter[0] & score[k]<scoreFilter[1]))
        {
        //Keep track of the sum of the weights to renormalize the density
          sumW+=w[k];
          if(REAL(strand)[0]==1)
          {
            gsl_matrix_set(Density,i,k,w[k]*gsl_ran_tdist_pdf((REAL(x)[i]-mu[k]+delta[k]/2.0)/sqrt(sF[k]), 4.0)/sqrt(sF[k]));
          }
          else if(REAL(strand)[0]==-1)
          {
            gsl_matrix_set(Density,i,k,w[k]*gsl_ran_tdist_pdf((REAL(x)[i]-mu[k]-delta[k]/2.0)/sqrt(sR[k]), 4.0)/sqrt(sR[k]));
          }
          else if(REAL(strand)[0]==0)
          {
            //Here I compute the shift density
            gsl_matrix_set(Density,i,k,0.5*w[k]*gsl_ran_tdist_pdf((REAL(x)[i]-mu[k])/sqrt(sR[k]), 4.0)/sqrt(sR[k]));
            gsl_matrix_set(Density,i,k,gsl_matrix_get(Density,i,k)+0.5*w[k]*gsl_ran_tdist_pdf((REAL(x)[i]-mu[k])/sqrt(sF[k]), 4.0)/sqrt(sF[k]));
          }
          // Here I scale by the score
          if(LOGICAL(scale)[0])
          {
            gsl_matrix_set(Density,i,k,gsl_matrix_get(Density,i,k)*score[k]);
          }
        }
      }
    }
      // We compute the overall density
    if(LOGICAL(sum)[0])
    {
          //Sum up over all binding events
      view=gsl_vector_view_array(REAL(y),length);
      gsl_blas_dgemv(CblasNoTrans, 1.0, Density, One, 0, &view.vector);
          //Rescale
      if(sumW>0)
      {
        gsl_vector_scale(&view.vector, 1./sumW);
      }          
    }//Otherwise we take the most likely event at each position
    else
    {
      for(i=0;i<length;i++)
      {
        view=gsl_matrix_row(Density,i);
        REAL(y)[i]=gsl_vector_max(&view.vector);
      }
    }

    SET_VECTOR_ELT(ans,0,x);
    SET_VECTOR_ELT(ans,1,y);
    PROTECT(names = allocVector(STRSXP, 2)); nProtected++; 
    SET_STRING_ELT(names, 0, mkChar("x"));
    SET_STRING_ELT(names, 1, mkChar("density"));
    setAttrib(ans, R_NamesSymbol, names);  

    gsl_vector_free(One);
    gsl_matrix_free(Density);
    UNPROTECT(nProtected++);
    return(ans);
  }
}

SEXP getListElement(SEXP list, const char *str)
     {
       SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
       int i;
     
       for (i = 0; i < length(list); i++)
         if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
           elmt = VECTOR_ELT(list, i);
           break;
         }
       return elmt;
     }
