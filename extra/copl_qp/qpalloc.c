#include "QPdef.h"

static int PrimeSet[_N_PRIM]={
                 29,  229,    883,    1671,   2791,
               4801, 8629,  15289,   25303,  34843,
              65269, 99709,129403,  147673, 166669,
             201403,222163,242729,  261431, 303491,
             320237,402761,501131,  602309, 701507,
             800999,900551,1000619,1100837,1200359
           };
           
char* cAlloc(int n)
{
  char* r=NULL;
  
  if (n<=0)
    return r;

  r=(char*)calloc(n,sizeof(char));
  if (!r) {
    ShutDown(NOT_MEMSPC);
    exit(0);
  }
  return r;
} /* cAlloc */

void cFree(char** x)
{
  char* r=*x;

  if (!r)
    return;

  free(r);
  *x=NULL;
} /* cFree */

int* iAlloc(int n)
{
  int* r=NULL;
  
  if (n) {
    r=(int*)calloc(n,sizeof(int));
    if (!r) {
      ShutDown(NOT_MEMSPC);
      exit(0);
    }
  }
  return r;
} /* iAlloc */

void iFree(int** x)
{
  int* r=*x;
  
  if (r) {
    free(r);
    *x=NULL;
  }
} /* iFree */

double* dAlloc(int n)
{
  double* r=NULL;
  
  if (n) {
    r=(double*)calloc(n,sizeof(double));
    if (!r) {
      ShutDown(NOT_MEMSPC);
      exit(0);
    }
  }
  return r;
} /* dAlloc */

void dFree(double** x)
{
  double* r=*x;
  
  if (r) {
    free(r);
    *x=NULL;
  }
} /* dFree */

char** cPallo(int m,int n)
{
  int    i;
  char** r=NULL;
  
  if (m>0) {
    r=(char**)calloc(m,sizeof(char*));
    if (!r) {
      ShutDown(NOT_MEMSPC);
      exit(0);
    }
    if (n>0) {
      r[0]=cAlloc(n*m);
      for (i=1; i<m; i++)
        r[i]=r[i-1]+n;
    }
  }
  return r;
} /* cPallo */

void cPfre(char ***x,int n)
{
  int    i;
  char** r=*x;
  
  if (r) {
    for (i=0; i<n; i++)
      cFree(&r[i]);
    free(r);
    *x=NULL;
  }
} /* cPfre */

int** iPallo(int m,int n)
{
  int   i;
  int** r=NULL;

  if (m>0) {
    r=(int**)calloc(m,sizeof(int*));
    if (!r) {
      ShutDown(NOT_MEMSPC);
      exit(0);
    }
    if (n>0) {
      r[0]=iAlloc(n*m);
      for (i=1; i<m; i++)
        r[i]=r[i-1]+n;
    }
  }
  return r;
} /* iPallo */

void iPfre(int*** x)
{
  int** r=*x;

  if (r) {
    iFree(&r[0]);
    free(r);
    *x=NULL;
  }
} /* iPfre */

double** dPallo(int m,int n)
{
  int      i;
  double** r=NULL;

  if (m>0) {
    r=(double**)calloc(m,sizeof(double*));
    if (!r) {
      ShutDown(NOT_MEMSPC);
      exit(0);
    }
    if (n>0) {
      r[0]=dAlloc(n*m);
      for (i=1; i<m; i++)
        r[i]=r[i-1]+n;
    }
  }
  return r;
} /* dPallo */

void dPfre(double*** x)
{
  double** r=*x;

  if (r) {
    dFree(&r[0]);
    free(r);
    *x=NULL;
  }
} /* dPfre */

hashtab* hAlloc(int lssz)
{
  int      i;
  hashtab* tab;
  
  tab=(hashtab*)malloc(sizeof(hashtab));
  if (!tab) {
    ShutDown(NOT_MEMSPC);
    exit(0);
  }
  
  if (lssz) {
    for (i=0; i<_N_PRIM; i++) {
      if (lssz<PrimeSet[i]) {
        lssz=PrimeSet[i];
        break;
      }
    }
    
    if (i>_N_PRIM)
    {
      ShutDown(NOT_PRIM);
      exit(0);
    }
    
    tab->lssz=lssz;
    tab->list=(plst*)calloc(lssz,sizeof(plst));
    if (!tab->list) {
      ShutDown(NOT_MEMSPC);
      exit(0);
    }
    
    for (i=0; i<lssz; i++)
      tab->list[i]=NULL;
  }
  return tab;
} /* hAlloc */

void hFree(hashtab **tab)
{
  int      i;
  plst     ptr,next;
  hashtab* r=*tab;
  
  if (!r) return;
  
  if (r->list) {
    for (i=0; i<r->lssz; i++) {
      ptr=r->list[i];
      while(ptr) {
        next=ptr->next;
        cFree(&ptr->key);
        cFree((char**)&ptr);
        ptr=next;
      }
    }    
    free(r->list);
    r->list=NULL;
  }
  free(r);
  *tab=NULL;
} /* hFree */

matrix *mAlloc(int  m,
               int  nz)
{
  matrix* r;
  
  r=(matrix*)calloc(1,sizeof(matrix));
  if (!r) {
    ShutDown(NOT_MEMSPC);
    exit(0);
  }
  
  if (m) {
    r->ia=(array*)calloc(m,sizeof(array));
    if (!r) {
      ShutDown(NOT_MEMSPC);
      exit(0);
    }
    
    if (nz) {
      r->ia->ja=iAlloc(nz);
      r->ia->an=dAlloc(nz);
    }
  }
  r->mnrs=m;
  r->nrow=m;
  r->mems=nz;
  return r;
} /* mAlloc */

void mFree(matrix** x)
{
  matrix* r=*x;
  
  if (r) {
    if (r->ia) {
      iFree(&r->ia->ja);
      dFree(&r->ia->an);
      free(r->ia);
      r->ia=NULL;
    }
    r->mnrs=0;
    r->nrow=0;
    free(r);
    *x=NULL;
  }
} /* mFree */

symatx* sAlloc(int n,int nz)
{
  int     i;
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
    else
    {
      for (i=0; i<n; i++)
      {
        r->cols[i].nn0=0;
        r->cols[i].ja=NULL;
        r->cols[i].an=NULL;
      }
    }
  }

  return r;
} /* sAlloc */

void sFree(symatx** x)
{
  int     i;
  symatx* r=*x;

  if (r)
  {
    if (r->cols)
    {
      if (r->nnzo)
      {
        iFree(&r->cols->ja);
        dFree(&r->cols->an);
      }
      else
      {
        for (i=0; i<r->ncol; i++)
        {
          iFree(&(r->cols[i].ja));
          dFree(&(r->cols[i].an));
        }
      }
      free(r->cols);
      r->cols=NULL;
    }
    r->ncol=0;
    r->nnzo=0;
    dFree(&r->diag);
    free(r);
    *x=NULL;
  }
} /* sFree */

optrow* rAlloc(int n)
{
  optrow* r=NULL;

  if (n>0) {
    r=(optrow*)calloc(n,sizeof(optrow));
    if (!r) {
      ShutDown(NOT_MEMSPC);
      exit(0);
    }
  }
  return r;
} /* rAlloc */

void rFree(optrow** x)
{
  int     i;
  optrow* r=*x;

  if (r) 
  {
    for (i=0; i<qt->m_intCol; i++)
    {
      iFree(&(r[i].ja));
      dFree(&(r[i].an));
    }
    free(r);
    *x=NULL;
  }
} /* rFree */

chol* chAllo(int n)
{
  chol* r=NULL;
  
  if (n)
  {
    r=(chol*)calloc(1,sizeof(chol));
    if (!r)
    {
      ShutDown(NOT_MEMSPC);
      exit(0);
    }
    
    r->mrow =n;
    r->nrow =n;
      
    r->snnz =0;
    r->shead=iAlloc(n);
    r->ssize=iAlloc(n);
    r->ssub =NULL;
    r->diag =dAlloc(n);
    r->unnz =0;
    r->ujnz =0;
    r->ujbeg=iAlloc(n);
    r->uhead=iAlloc(n);
    r->ujsze=iAlloc(n);
    r->usub =NULL;
    r->uval =NULL;
    r->perm =iAlloc(n);
    r->invp =iAlloc(n);
    r->nsnds=0;
    r->subg =iAlloc(n+1);      
  }
  return r;
} /* chAllo */

void chFre(chol** x)
{
  chol* r=*x;
  
  if (r)
  {
    iFree(&r->shead);
    iFree(&r->ssize);
    iFree(&r->ssub);
    dFree(&r->diag);
    iFree(&r->ujbeg);
    iFree(&r->uhead);
    iFree(&r->ujsze);
    iFree(&r->usub);
    dFree(&r->uval);
    iFree(&r->perm);
    iFree(&r->invp);
    iFree(&r->subg);
    iFree(&r->dhead);
    iFree(&r->dbeg);
    iFree(&r->dsub);
    free(r);
    *x=NULL;
  }
} /* chFre */

order *oAlloc(int  nnod,
              int  nn0)
{
  order* r;
  
  r=(order*)calloc(1,sizeof(order));
  if (!r)
  {
    ShutDown(NOT_MEMSPC);
    exit(0);
  }
  
  r->nnod=nnod;
  r->nn0 =nn0;
  
  r->adjn=iAlloc(nn0);
  r->rbeg=iAlloc(nnod);
  r->rexs=iAlloc(nnod);
  r->rlen=iAlloc(nnod);
  r->rend=iAlloc(nnod);
  r->pres=iAlloc(nnod);
  r->succ=iAlloc(nnod);
  
  return (r);
} /* oAlloc */

void oFree(order** x)
{
  order* r=*x;
  
  if (r) {
    iFree(&r->adjn);
    iFree(&r->rbeg);
    iFree(&r->rexs);
    iFree(&r->rlen);
    iFree(&r->rend);
    iFree(&r->pres);
    iFree(&r->succ);
    free(r);
    *x=NULL;  
  }
} /* OdFree */

static void xClear(xlist* xt)
{
  int i,sze;
  
  sze     =xt->last;
  xt->idep=xt->most+1;
  xt->lowp=xt->idep;
  xt->cure=sze;
  xt->ntot=0;
  
  for (i=0; i<xt->idep; i++)
    xt->head[i]=xt->last;
    
  for (i=0; i<sze; i++) {
    xt->port[i]=xt->idep;
    xt->fwrd[i]=xt->last;
    xt->bwrd[i]=xt->last;
  }
} /* xClear */

xlist* xAlloc(int  last,
              int  most)
{
  xlist* r;
  
  r=(xlist*)calloc(1,sizeof(xlist));
  if (!r) 
  {
    ShutDown(NOT_MEMSPC);
    exit(0);
  }
    
  r->loca=true;
  r->last=last;
  r->most=most;
  r->ntot=0;
    
  r->head=iAlloc(most+1);
  r->port=iAlloc(last);
  r->fwrd=iAlloc(last);
  r->bwrd=iAlloc(last);
  
  xClear(r);
  
  return (r);
} /* xAlloc */

void xFree(xlist** x)
{
  xlist* r=*x;
  
  if (r) {
    if (r->loca) {
      iFree(&r->head);
      iFree(&r->port);
      iFree(&r->fwrd);
      iFree(&r->bwrd);
    }
    free(r);
    *x=NULL;
  }
} /* XtFree */

void itAllo(int  m,
            int  n,
            int* x[])
{
  int i;
  
  if (n) {
    for (i=0; i<m; i++) {
      x[i]=(int*)calloc(n,sizeof(int));
      if (!x[i])
      {
        ShutDown(NOT_MEMSPC);
        exit(0);
      }
    }
  }
} /* itAllo */

void itFre(int  m,
           int* x[])
{
  int i;
  
  for (i=0; i<m; i++)
    iFree(&x[i]);
} /* itFre */

