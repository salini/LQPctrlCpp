#include "QPdef.h"
    
void iZero(int n,
           int *x,
           int *s)
{
  int i;

  if (s) {
    for(i=0; i<n; ++i)
      x[s[i]]=0;
  }
  else
    memset(x,0,n*sizeof(n));
} /* iZero */

void iSet(int n,
          int a,
          int *x,
          int *id)
{
  int i;

  if (!id)
    for(i=0; i<n; ++i)
      x[i]=a;
  else
    for(i=0; i<n; ++i)
      x[id[i]]=a;
} /* iSet */

void iSwap(int i1,
           int i2,
           int *v)
{
  int temp;

  temp  = v[i1];
  v[i1] = v[i2];
  v[i2] = temp;
} /* iSwap */

void iCopy(int n,
           int *s,
           int *d)
{
  memcpy(d,s,n*sizeof(int));
} /* iCopy */

int iSum(int n,
         int *x)
{
  int i,sum=0;
  
  for (i=0; i<n; i++)
    sum+=x[i];
    
  return sum;
} /* iSum */

void dCopy(int    n,
           double *s,
           double *d)
{
  if (n) memcpy(d,s,n*sizeof(double));
} /* dCopy */

void dCat(int    n,
          int    *ix,
          double *s,
          double *d)
{
  int i;
  
  for (i=0; i<n; i++) {
    d[i]=s[ix[i]];
    s[ix[i]]=0.0;
  }
} /* dCat */

double dDot(double *x,
            double *y,
            int n)
{
  int    i;
  double r=0.0;
  
  for (i=0; i<n; i++)
    r+=x[i]*y[i];
    
  return r;
} /* dDot */

double svDot(array*  x,
             double* y)
{
  int    i;
  double r=0.0;
  
  for(i=0; i<x->nn0; ++i)
    r+=x->an[i]*y[x->ja[i]];
    
  return (r);
} /* svDot */

void setArray(double  alf,
              array*  x,
              double* y)
{
  int k;
  
  if (x->nn0<=0||alf==0.0)
    return;
  
  if (alf==1.0) {
    for (k=0; k<x->nn0; k++)
      y[x->ja[k]]+=x->an[k];
  }
  else if (alf==-1.0) {
    for (k=0; k<x->nn0; k++)
      y[x->ja[k]]-=x->an[k];
  }
  else {
    for (k=0; k<x->nn0; k++)
      y[x->ja[k]]+=alf*x->an[k];
  }
} /* setArray */
