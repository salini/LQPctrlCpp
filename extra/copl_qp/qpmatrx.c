#include "QPdef.h"

static symatx* qAlloc(int n,int nz)
{
  symatx* r=NULL;
  
  r=(symatx*)calloc(1,sizeof(symatx));
  if (!r)
  {
    ShutDown(NOT_MEMSPC);
    exit(0);
  }
  if (n>0)
  {
    r->ncol=n;
    r->diag=dAlloc(n);
    r->cols=(array*)calloc(n,sizeof(array));
    if (!r->cols)
    {
      ShutDown(NOT_MEMSPC);
      exit(0);
    }
    if (nz)
    {
      r->nnzo    =nz;
      r->cols->ja=iAlloc(nz);
      r->cols->an=dAlloc(nz);
    }
  }

  return r;
} /* sAlloc */

static void makeQmatx(symatx* q,
                      int     nnzo)
{
  int i;
  
  if (q->cols&&q->nnzo)
  {
    iFree(&q->cols->ja);
    dFree(&q->cols->an);
    q->cols->ja=iAlloc(nnzo);
    q->cols->an=dAlloc(nnzo);
    q->nnzo    =nnzo;
  }
  else if (q->cols)
  {
    for (i=0; i<q->ncol; i++)
    {
      if (q->cols+i&&q->cols[i].nn0)
      {
        iFree(&q->cols[i].ja);
        dFree(&q->cols[i].an);
        q->cols[i].nn0=0;
      }
      q->cols->ja=iAlloc(nnzo);
      q->cols->an=dAlloc(nnzo);
      q->nnzo    =nnzo;
    }
  }
} /* makeQmatx */

static void InitQmatx(symatx *q)
{
  int i;
  
  for (i=0; i<q->ncol; i++)
    q->cols[i].nn0=0;
} /* InitQmatx */

static void ScanQmatx(symatx* q,
                      symatx* qt)
{
  int i,j,k,ncol=q->ncol,nn0;
  
  for (i=0; i<ncol; i++) {
    nn0=q->cols[i].nn0;
    for (k=0; k<nn0; k++) {
      j=q->cols[i].ja[k];
      qt->cols[j].nn0++;
    }
  }
} /* ScanQmatx */

int qTrans(int      ncol,
           symatx*  q,
           int      pid,
           int*     p,
           int      vid,
           symatx** qt)
{
  int    i,j,k,t,nn0;
  symatx *r;
  array  *qi,*qj;
  
  nn0=0;
  for (i=0; i<q->ncol; i++)
    nn0+=q->cols[i].nn0;
  
  r=*qt;
  
  if (r) {
    if (!vid&&nn0>r->nnzo)
    {
      ShutDown(NOT_MEMSPC);
      exit(0);
    }
    if (ncol>r->ncol)
    {
      ShutDown(NOT_MEMSPC);
      exit(0);
    }
    r->ncol=ncol;
  }

  else
  {
    r=qAlloc(ncol,0);
    makeQmatx(r,nn0);
  }
  
  if (vid)
    for (j=0; j<ncol; j++)
      r->cols[j].nn0=0;
  
  else
  {
    InitQmatx(r);
    ScanQmatx(q,r);
    
    nn0=0;
    for (j=0; j<ncol; j++) {
      r->cols[j].ja=r->cols->ja+nn0;
      r->cols[j].an=r->cols->an+nn0;
      nn0+=r->cols[j].nn0;
      r->cols[j].nn0=0;
    }
    
    if (nn0>r->nnzo)
    {
      ShutDown(NOT_MEMSPC);
      exit(0);
    }
    r->nnzo=nn0;
  }
  
  for (k=0; k<q->ncol; k++)
  {
    i=k;
    if (p)
    {
      qi=q->cols+p[k];
      if (!pid) i=p[k];
    }
    else qi=q->cols+k;
    
    for (t=0; t<qi->nn0; t++) {
      j  =qi->ja[t];
      qj =r->cols+j;
      nn0=qj->nn0;
      
      qj->ja[nn0]=i;
      qj->an[nn0]=qi->an[t];
      qj->nn0++;
    }
  }
  
  *qt=r;
  
  return true;
} /* qTrans */

