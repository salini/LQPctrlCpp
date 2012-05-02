#include "QPdef.h"

static void SetChol(qpdat* opt)
{
  int   h,j,k,m,n,nn,nm;
  array *aj,*qj;
  chol* cl;
  
  m =opt->m_intRow;
  n =opt->m_colUse;
  nn=opt->m_intCol;
  nm=m+nn;
  cl=opt->cl;
  
  h=0;
  for (j=0; j<n; j++)
  {
    qj=opt->qh->cols+j;
    for (k=0; k<qj->nn0; k++)
    {
      cl->ssub[h]=qj->ja[k];
      h++;
    }
    aj=opt->at->ia+j;
    for (k=0; k<aj->nn0; k++)
    {
      cl->ssub[h]=nn+aj->ja[k];
      h++;
    }
  }
  for (j=n; j<nn; j++)
  {
    aj=opt->at->ia+j;
    for (k=0; k<aj->nn0; k++)
    {
      cl->ssub[h]=nn+aj->ja[k];
      h++;
    }
  }
} /* SetChol */

static void plusXs(int  n,
                   int* x,
                   int* s)
{
  int i;
  
  if (!s) {
    for (i=0; i<n; i++)
      x[i]++;
  }
  else {
    for (i=0; i<n; i++)
      x[s[i]]++;
  }
} /* plusXs */

static int LocIntPos(int n,
                     int i,
                     int *v)
{
  int j;
  
  for(j=0; j<n && i!=v[j]; ++j);
  return (j);
} /* LocIntPos */

static void InsertDplRow(int i,
                         int ifir,
                         int *ifirst,
                         int *ilink)
{
  int temp;
  
  temp=ifirst[ifir];
  ifirst[ifir]=i;
  ilink[i]=temp;
} /* InsertDplRow */

static void PermTransSym(int nrow,
                         int *fir,
                         int *nnz,
                         int *sub,
                         int *p,
                         int rowwise,
                         int *firt,
                         int *nnzt,
                         int *subt)
{
  int i,j,s,t,stopt;

  iZero(nrow,nnzt,NULL);
  
  if (rowwise) {
    if (p) {
      for(s=0; s<nrow; ++s) {
        j =p[s];
        for(t=fir[s], stopt=t+nnz[s]; t<stopt; ++t) {
          i=p[sub[t]];
          nnzt[max(i,j)]++;
        }
      }
    }
    else {
      for(j=0; j<nrow; ++j) {
        for(t=fir[j], stopt=t+nnz[j]; t<stopt; ++t) {
          i=sub[t];
          nnzt[max(i,j)]++;
        }
      }
    }
  }
  
  else {
    if (p) {
      for(s=0; s<nrow; ++s) {
        j =p[s];
        for(t=fir[s], stopt=t+nnz[s]; t<stopt; ++t) {
          i=p[sub[t]];
          nnzt[min(i,j)]++;
        }
      }
    }
    else {
      for(j=0; j<nrow; ++j) {
        for(t=fir[j], stopt=t+nnz[j]; t<stopt; ++t) {
          i=sub[t];
          nnzt[min(i,j)]++;
        }
      }
    }
  }

  firt[0]=0;
  for(i=1; i<nrow; ++i) {
    firt[i]  =firt[i-1] + nnzt[i-1];
    nnzt[i-1]=0;
  }
  nnzt[nrow-1]=0;

  if (rowwise) {
    if (p) {
      for(s=0; s<nrow; ++s) {
        j=p[s];
        for(t=fir[s], stopt=t+nnz[s]; t<stopt; ++t) {
          i=p[sub[t]];
          if (i>j) {
            subt[firt[i]+nnzt[i]]=j;
            nnzt[i]++;
          }
          else {
            subt[firt[j]+nnzt[j]]=i;
            nnzt[j]++;
          }
        }
      }
    }
    else {
      for(j=0; j<nrow; ++j) {
        for(t=fir[j], stopt=t+nnz[j]; t<stopt; ++t) {
          i=sub[t];
          if (i>j) {
            subt[firt[i]+nnzt[i]]=j;
            nnzt[i]++;
          }
          else {
            subt[firt[j]+nnzt[j]]=i;
            nnzt[j]++;
          }
        }
      }
    }
  }
  
  else {
    if (p) {
      for(s=0; s<nrow; ++s) {
        j=p[s];
        for(t=fir[s], stopt=t+nnz[s]; t<stopt; ++t) {
          i=p[sub[t]];
          if (i<j) {
            subt[firt[i]+nnzt[i]]=j;
            nnzt[i]++;
          }
          else {
            subt[firt[j]+nnzt[j]]=i;
            nnzt[j]++;
          }
        }
      }
    }
    else {
      for(j=0; j<nrow; ++j) {
        for(t=fir[j], stopt=t+nnz[j]; t<stopt; ++t) {
          i=sub[t];
          if (i<j) {
            subt[firt[i]+nnzt[i]]=j;
            nnzt[i]++;
          }
          else {
            subt[firt[j]+nnzt[j]]=i;
            nnzt[j]++;
          }
        }
      }
    }
  }
} /* PermTransSym */

static void LocDplRow(int nrow,
                      int ncol,
                      int fir[],
                      int sze[],
                      int *sub,
                      int map[],
                      int ifirst[],
                      int ilink[],
                      int ilist[],
                      int *icount,
                      int i1nrow[])
{
  int i,new,n,oisze,isze,s,count,
      temp,k,nexti,*cur=i1nrow;
  
  n      =nrow;
  *icount=0;
  
  for (i=0; i<nrow; i++) {
    cur[i]  =0;
    ilink[i]=n;
    ilist[i]=n;
  }
  iSet(ncol,n,ifirst,NULL);
  
  isze =0;
  count=0;
  oisze=isze;
  new  =n;
  for(i=0; i<nrow; ++i) {
    if (map)
      for(; cur[i]<sze[i] && !map[sub[fir[i]+cur[i]]]; ++cur[i]);
    
    if ( cur[i]<sze[i] ) {
      s=sub[fir[i]+cur[i]];
      if ( ifirst[s]==n )
        ilist[isze++]=s;
      
      InsertDplRow(i,s,ifirst,ilink);
      
      cur[i]++;
    }
    
    else {
      temp    =new;
      new     =i;
      ilink[i]=temp;
    }
  }
  
  for(k=oisze; k<isze; ++k) {
    temp            =ifirst[ilist[k]];
    ifirst[ilist[k]]=n;
    ilist[k]        =temp;
  }
  
  if (new!=n) {
    count++;
    ilist[nrow-count]=new;
  }
  
  while(isze) {
    isze--;
    oisze      =isze;
    
    i          =ilist[isze];
    ilist[isze]=n;
    
    if (i==n)
      exit(0);
    
    new=n;
    if (ilink[i]==n)
      new=i;
    else {
      for(; i!=n; i=nexti) {
        nexti   =ilink[i];
        ilink[i]=n;
        
        if ( map )
          for(; cur[i]<sze[i] && !map[sub[fir[i]+cur[i]]]; ++cur[i]);
          
        if (cur[i]<sze[i]) {
          s =sub[fir[i]+cur[i]];
          cur[i]++;
           
          if (ifirst[s]==n)
            ilist[isze++]=s;
           
          temp     =ifirst[s];
          ifirst[s]=i;
          ilink[i] =temp;
        }
        
        else {
          temp    =new;
          new     =i;
          ilink[i]=temp;
        }
      }
    }
    
    for(k=oisze; k<isze; ++k) {
      temp            =ifirst[ilist[k]];
      ifirst[ilist[k]]=n;
      ilist[k]        =temp;
    }
    
    if (new!=n) {
      count++;
      ilist[nrow-count]=new;
    }
  }
  
  *icount=count;
  for(k=0; k<count; ++k)
    ilist[k]=ilist[nrow-count+k];
} /* LocDplRow */

static int CompIntElem(const void* e1,
                       const void* e2)
{
  int *i1,*i2;

  i1=(int *) e1;
  i2=(int *) e2;

  if (*i1<*i2)
    return (-1);
  else if(*i1>*i2)
    return (1);
  return (0);
} /* CompIntElem */

static void iSort(int  n,
                  int* x)
{
  qsort((void *)x,n,sizeof(int),CompIntElem);
} /* iSort */

static void DetectDenseNodes(chol* sf,
                             int*  i1nrow,
                             int*  i2nrow,
                             int*  i3nrow,
                             int*  i4nrow,
                             int*  i5nrow,
                             int*  i6nrow)
{
  int j,k,l,t,ndens,nil=sf->nrow,
      *subg=sf->subg,*ujbeg=sf->ujbeg,
      *ujsze=sf->ujsze,*usub=sf->usub,
      *fir,*sze,*ilist,*ilink;
  
  if (!sf->nsnds||
      !i1nrow  || !i2nrow || !i3nrow ||
      !i4nrow  || !i5nrow || !i6nrow) {
    sf->sdens=false;
    return;
  }
  
  sf->sdens  =true;
  fir        =i1nrow;
  sze        =i2nrow;
  ilist      =i3nrow;
  ilink      =i4nrow;
  
  sf->nsndn=0;
  
  l=subg[sf->nsnds-1];
  for(k=0; k+1<sf->nsnds; ++k) {
    j=subg[k];
    for(t=0; t<ujsze[j] && usub[ujbeg[j]+t]<l; ++t);
    
    fir[k] =ujbeg[j]+t;
    sze[k] =ujsze[j]-t;
  }
  
  LocDplRow(sf->nsnds-1,sf->nrow,fir,sze,usub,
            NULL,
            i6nrow,ilink,ilist,&ndens,i5nrow);
  
  sf->dhead=iAlloc(ndens+1);
  sf->dsub=iAlloc(sf->nsnds);
  sf->dbeg=iAlloc(sf->nsnds);
  
  nil        =sf->nsnds-1;
  sf->ndens   =0;
  sf->nsndn   =0;
  sf->dhead[0]=0;
  
  for(k=0; k<ndens; ++k) {
    j=ilist[k];
    if ( sze[j] ) {
      sf->dhead[sf->ndens+1]=sf->dhead[sf->ndens];
      sf->ndens++;
      for(; j!=nil; j=ilink[j]) {
        sf->nsndn                    +=  sf->subg[j+1]-sf->subg[j];
        sf->dsub[sf->dhead[sf->ndens]]=j;
        sf->dbeg[sf->dhead[sf->ndens]]=fir[j]-ujbeg[subg[j]];
        sf->dhead[sf->ndens]++;
      }
      iSort(sf->dhead[sf->ndens]-sf->dhead[sf->ndens-1],
            sf->dsub+sf->dhead[sf->ndens-1]);
      
      for(t=sf->dhead[sf->ndens-1]; t<sf->dhead[sf->ndens]; ++t)
        sf->dbeg[t]=fir[sf->dsub[t]]-ujbeg[subg[sf->dsub[t]]];
    }
  }
} /* DetectDenseNodes */

static int ChlSymb(chol* sf,
                   int   ulnnz)
{
  int chksn,i,j,t,stopt,sze,first,cur,k,
      ffree=0,ipos,nrow=sf->nrow,nil=nrow,
      *nnz,*fir,*pssub,*link,*buf,*mask,
      *usub,*tusub,*i1nrow,*i2nrow,*i3nrow,
      *i4nrow,*p=sf->perm,*invp=sf->invp,
      *ujbeg=sf->ujbeg,*ujsze=sf->ujsze,
      *subg=sf->subg;
  
  pssub=iAlloc(sf->snnz);
  
  for(i=0; i<nrow; ++i)
    invp[p[i]]=i;
  
  nnz=sf->uhead;
  fir=sf->subg;
  
  PermTransSym(nrow,sf->shead,sf->ssize,sf->ssub,
               invp,true,fir,nnz,pssub);
  
  PermTransSym(nrow,fir,nnz,pssub,NULL,false,
               sf->shead,sf->ssize,sf->ssub);
  
  iFree(&pssub);
  
  k       =ulnnz+nrow;
  usub    =iAlloc(k);
  buf     =usub+ulnnz;
  
  mask=sf->uhead;
  
  link=invp;
  
  iZero(nrow,mask,NULL);
  iSet(nrow,nil,link,NULL);
  
  ffree   =0;
  sf->nsnds=0;
  subg[0] =0;
  for(i=0; i<nrow; ++i) {
    sze  =sf->ssize[i];
    first=nil;
    cur  =link[i];
    chksn=false;
    
    if (cur==nil) {
      
      subg[sf->nsnds+1] =1 + subg[sf->nsnds];
      ujsze[i]          =sze;
      ujbeg[i]         =ffree;
      ffree            += sze;
      
      iCopy(sze,sf->ssub+sf->shead[i],usub+ujbeg[i]);
      if (sze) {
        first=usub[ujbeg[i]];
        for(cur=first; link[cur]!=nil; cur=link[cur]);
        link[cur]     =sf->nsnds;
        link[sf->nsnds]=nil;
      }
      sf->nsnds++;
    }
    
    else {
      mask[i]=1;
      
      iCopy(sze,sf->ssub+sf->shead[i],buf);
      iSet(sze,1,mask,buf);
        
      for(; cur!=nil; cur=link[cur]) {
        chksn |= (1+cur==sf->nsnds);
        k     =subg[cur];
        
        for(t=ujbeg[k], stopt=t+ujsze[k]; t<stopt; ++t) {
          j=usub[t];
          if ( j>i && !mask[j] ) {
            buf[sze]=j;
            mask[j] =1;
            sze++;
          }
        }
      }
      
      if (chksn) {
        k    =subg[sf->nsnds-1];
        chksn=sze==( ujsze[k]-(subg[sf->nsnds]-subg[sf->nsnds-1]) );
      }
      
      first  =nrow;
      mask[i]=0;
      for(t=0; t<sze; ++t) {
        j     =buf[t];
        mask[j]=0;
        first  =min(j,first);
      }
      
      if (chksn) {
        ipos=LocIntPos(ujsze[i-1],i,usub+ujbeg[i-1]);
         
        if (ipos==ujsze[i-1])
        {
          ShutDown(SYS_ERROR);
          exit(0);
        }
        
        iSwap(ujbeg[i-1],ipos+ujbeg[i-1],usub);

        subg[sf->nsnds]++;
        ujbeg[i] =ujbeg[i-1]+1;
        ujsze[i]  =ujsze[i-1]-1;
        
        if (usub[ujbeg[i]-1]!=i)
        {
          ShutDown(SYS_ERROR);
          exit(0);
        }
 
        if ( first!=nil ) {
          for(cur=first; link[cur]!=nil; cur=link[cur]);
          link[cur]       =sf->nsnds-1;
          link[sf->nsnds-1]=nil;
        }
      }
      
      else {
        subg[sf->nsnds+1] =1 + subg[sf->nsnds];
        ujbeg[i]         =ffree;
        ujsze[i]          =sze;
        ffree            += sze;
        
        if (ffree>ulnnz)
        {
          ShutDown(SYS_ERROR);
          exit(0);
        }
        
        iCopy(sze,buf,usub+ujbeg[i]);
        
        if ( first!=nil ) {
          for(cur=first; link[cur]!=nil; cur=link[cur]);
          link[cur]     =sf->nsnds;
          link[sf->nsnds]=nil;
        }
        sf->nsnds++;
      }
    }
     
    if (ujsze[i]+1==nrow-i)
      break;
  }
  
  for(++i; i<nrow; ++i) {
    ujsze[i] =ujsze[i-1]-1;
    ujbeg[i]=ujbeg[i-1]+1;
     
    subg[sf->nsnds]++;
  }
   
  tusub=iAlloc(ffree);
    
  fir=buf;
  nnz=sf->uhead;
  
  iZero(nrow,nnz,NULL);
  
  for(k=0; k<sf->nsnds; ++k) {
    j=subg[k];
    plusXs(ujsze[j],nnz,usub+ujbeg[j]);
  }
  
  fir[0]=0;
  for(k=1; k<nrow; ++k)
    fir[k]=fir[k-1] + nnz[k-1];
    
  iZero(nrow,nnz,NULL);
  
  for(k=0; k<sf->nsnds; ++k) {
    j=subg[k];
    for(t=ujbeg[j], stopt=t+ujsze[j]; t<stopt; ++t) {
      i                    =usub[t];
      tusub[fir[i]+nnz[i]] =j;
      nnz[i]++;
    }
    ujsze[j]=0;
  }
  
  for(i=0; i<nrow; ++i) {
    for(t=fir[i], stopt=t+nnz[i]; t<stopt; ++t) {
      j                      =tusub[t];
      usub[ujbeg[j]+ujsze[j]] =i;
      ujsze[j]++;
    }
  }
  
  iFree(&tusub);

  if (ffree<=sf->ujnz) {
    iCopy(ffree,usub,sf->usub);
    iFree(&usub);
  }
  
  else {
    sf->ujnz=0;
    iFree(&sf->usub);
    
    sf->usub=iAlloc(ffree);
    iCopy(ffree,usub,sf->usub);
    
    sf->ujnz=ffree;
    iFree(&usub);
  }
  
  i1nrow=iAlloc(4*nrow);
  i2nrow=i1nrow+nrow;
  i3nrow=i2nrow+nrow;
  i4nrow=i3nrow+nrow;
  
  DetectDenseNodes(sf,sf->uhead,sf->invp,
                   i1nrow,i2nrow,i3nrow,i4nrow);
  
  iFree(&i1nrow);
  
  sf->uhead[0]=0;
  for(i=1; i<nrow; ++i)
    sf->uhead[i]=sf->uhead[i-1]+sf->ujsze[i-1];
  
  for(i=0; i<nrow; ++i)
    invp[p[i]]=i;
  
  for(k=0; k<sf->nsnds; ++k)
    if ( subg[k]+1!=subg[k+1] )
      break;
      
  return true;
} /* ChlSymb */

static void PrintDim(qpdat* opt)
{
  char sl[_L_SIZE],sr[_L_SIZE];
  
  fprintf(fout,"KKT Matrix Statistics\n");
  fprintf(fout,"---------------------\n");
  
  sprintf(sl,"   number of rows in KKT matrix");
  sprintf(sr,"%d",opt->cl->nrow);
  LeftDots(sl,sr);
   
  sprintf(sl,"   number of nonzeros in KKT matrix");
  sprintf(sr,"%d",2*opt->cl->snnz+opt->ndia);
  LeftDots(sl,sr);
   
  sprintf(sl,"   density of KKT matrix");
  sprintf(sr,"%.3f%s",
             100.0*(((2.0*opt->cl->snnz+opt->ndia)/
             (double)opt->cl->nrow)/(double)opt->cl->nrow),
             "%");
  LeftDots(sl,sr);
   
  sprintf(sl,"   number of nonzeros in Cholesky factor");
  sprintf(sr,"%d",opt->cl->unnz+opt->cl->nrow);
  LeftDots(sl,sr);
   
  sprintf(sl,"   density of Cholesky factor");
  sprintf(sr,"%.3f%s",
             100.0*(((opt->cl->unnz+opt->cl->nrow)/
             (double)opt->cl->nrow)/(double)opt->cl->nrow),
             "%");
  LeftDots(sl,sr);
  
  fprintf(fout,"\n");

  printf("\n");

} /* PrintDim */

int SymboProc(qpdat* opt,
              qppar* par)
{
  int    idKey,i,k,m,n,t,nm,nn,
         snnz,lnnz,*nnzi;
  order* od;
  chol*  cl;
  
  printf(" BEGIN symbolic computation...\n");
  
  m =opt->m_intRow;
  n =opt->m_colUse;
  nn=opt->m_intCol;
  nm=nn+m;
  opt->cl=chAllo(nm);
  cl=opt->cl;
  
  snnz=opt->m_intQnz+opt->m_intAnz;
  if (!snnz)
    return SYS_ERROR;
  
  cl->ssub=iAlloc(snnz);
  cl->snnz=snnz;
  
  k=0;
  for (i=0; i<n; i++)
  {
    snnz=opt->qh->cols[i].nn0+opt->at->ia[i].nn0;
    cl->shead[i]=k;
    cl->ssize[i]=snnz;
    k+=snnz;
  }
  for (i=n; i<nn; i++)
  {
    snnz=opt->at->ia[i].nn0;
    cl->shead[i]=k;
    cl->ssize[i]=snnz;
    k+=snnz;
  }
  for (i=nn; i<nm; i++)
  {
    cl->shead[i]=k;
    cl->ssize[i]=0;
  }
  
  SetChol(opt);
    
  nnzi=cl->perm;
  memset(nnzi,0,nm*sizeof(int));
  for (i=0; i<nm; i++)
  {
    nnzi[i]+=cl->ssize[i];
    plusXs(cl->ssize[i],nnzi,cl->ssub+cl->shead[i]);
  }
  
  od=oAlloc(nm,2*cl->snnz);
  idKey=OdInit(od,nnzi);
  if (idKey!=OPT_OK)
    return idKey;
    
  for (i=0; i<nm; i++)
    for (t=0; t<cl->ssize[i]; ++t)
      OdIndex(od,i,cl->ssub[cl->shead[i]+t]);
  
  GetOrder(od,cl->perm);
  lnnz=od->ntot; 
  oFree(&od);
  
  idKey=ChlSymb(cl,lnnz);
  
  lnnz=iSum(cl->nrow,cl->ujsze);
  if (lnnz>cl->unnz)
  {
    cl->unnz=0;
    if (cl->uval)
      dFree(&cl->uval);
    cl->uval=dAlloc(lnnz);
    cl->unnz=lnnz;
  }
  printf(" END symbolic computation\n");  
  PrintDim(opt);
  
  return OPT_OK;
} /* SymboProc */
