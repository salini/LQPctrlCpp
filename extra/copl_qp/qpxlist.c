#include "QPdef.h"

static int ChkXlist(xlist *xt,
                    int   ind)
{
  return (xt->port[ind]==xt->idep);
} /* ChkXlist */

int XtSucc(xlist *xt)
{
  int t,last=xt->last,most=xt->most,
      *head;

  if (xt->cure==last)
    return (false);
    
  if (xt->fwrd[xt->cure]!=last)
    xt->cure=xt->fwrd[xt->cure];
      
  else {
    head=xt->head;
     
    for(t=xt->port[xt->cure]+1; t<=most && head[t]==last; ++t);
    
    if (t>most)
      xt->cure=last;
    else
      xt->cure=xt->head[t];
  }

  return (true);
} /* XtSucc */

void XtDel(xlist *xt,
           int    e)
{
  int p;

  if (!ChkXlist(xt,e)) {
    
    if (xt->ntot<=0)
    {
      ShutDown(SYS_ERROR);
      exit(0);
    }
     
    xt->ntot--;
    
    if (xt->cure==e) {
      if (xt->ntot)
        XtSucc(xt);
      else
        xt->cure=xt->last;
    }
    
    p          =xt->port[e];
    xt->port[e]=xt->idep;
    
    if (xt->bwrd[e]!=xt->last)
      xt->fwrd[xt->bwrd[e]]=xt->fwrd[e];
    else
      xt->head[p]=xt->fwrd[e];
    
    if (xt->fwrd[e]!=xt->last)
      xt->bwrd[xt->fwrd[e]]=xt->bwrd[e];

    if (xt->head[p]==xt->last &&
         xt->lowp==p) {
      xt->lowp=xt->idep;
      if (xt->ntot) {
        for(++p; p<=xt->most; ++p){
          if (xt->head[p]!=xt->last){
            xt->lowp=p;
            break;
          }
        }
      }
    }
  }
} /* XtDel */

void XtPut(xlist *xt,
           int   e,
           int   p)
{
  if (0<=e && e<xt->last && 0<=p && p<=xt->most) {
    XtDel(xt,e);
    xt->ntot++;
    xt->port[e] =p;
    xt->fwrd[e] =xt->head[p];
    xt->bwrd[e]=xt->last;
    
    if (xt->head[p]!=xt->last)
      xt->bwrd[xt->head[p]]=e;
    
    xt->head[p]=e;
    xt->lowp    =min(p,xt->lowp);
  }
  
  else
  {
    ShutDown(SYS_ERROR);
    exit(0);
  }
} /* XtPut */

int XtLeast(xlist *xt)
{
  if (xt->lowp==xt->idep) {
    if (xt->ntot!=0)
    {
      ShutDown(SYS_ERROR);
      exit(0);
    }
    
    xt->cure=xt->last;
    return false;
  }
  
  else {
    if (xt->ntot<=0)
    {
      ShutDown(SYS_ERROR);
      exit(0);
    }
    
    xt->cure=xt->head[xt->lowp];
    return true;
  }
} /* XtLeast */

int XtGet(xlist *xt,
          int   *e,
          int   *p)
{
  if (xt->cure>xt->last)
  {
    ShutDown(SYS_ERROR);
    exit(0);
  }
  
  if (xt->cure==xt->last)
    return false;
    
  *e=xt->cure;
  *p=xt->port[*e];
  
  return true;
} /* XtGet */
