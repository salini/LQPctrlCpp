#include "QPdef.h"

static void DataSwap(qpdat* opt)
{
  int     h,i,j,k,n;
  optrow* col;
  array*  aj;
  
  h  =0;
  n  =opt->m_colUse;
  col=opt->col;
  
  memset(opt->imem,0,opt->m_intRow*sizeof(int));

  for (j=0; j<n; j++)
  {
    aj=opt->at->ia+j;
    aj->ja=opt->at->ia->ja+h;
    aj->an=opt->at->ia->an+h;
    
    k        =col[j].nn0;
    aj->nn0  =k;
    opt->c[j]=col[j].rhs;
    opt->u[j]=col[j].upp;
    opt->l[j]=col[j].low;
    
    for (i=0; i<k; i++)
    {
      opt->imem[col[j].ja[i]]++;
      aj->ja[i]=col[j].ja[i];
      aj->an[i]=col[j].an[i];
    } 
    h+=k;
    iFree(&(col[j].ja));
    dFree(&(col[j].an));
  }
  rFree(&opt->col);
} /* DataSwap */

static void TreatBound(double bound,
                       int    indx,
                       qppar* par,
                       qpdat* opt,
                       int    type)
{
  int    i,j;
  double rtmp;
  array* aj;
  
  aj=opt->at->ia+indx;

  switch(type)
  {
  case LOW:
    for (j=0; j<aj->nn0; j++)
    {
      i         = aj->ja[j];
      rtmp      =-aj->an[j]*bound;
      opt->b[i]+=rtmp;
      opt->r[i]+=rtmp;
    }
    opt->obj +=opt->c[indx]*bound;
    if (opt->u[indx]<=par->q_dblMax)
      opt->u[indx]-=bound;
    break;

  case MIS:
    opt->l[indx]=-bound;
    opt->u[indx]=_P_MAXM;
    opt->c[indx]=-opt->c[indx];
      
    if (fabs(bound)>=par->q_dblAij)
    {
      for (j=0; j<aj->nn0; j++)
      {
        i         = aj->ja[j];
        aj->an[j] =-aj->an[j];
        rtmp      = aj->an[j]*bound;
        opt->b[i]+= rtmp;
        opt->r[i]+= rtmp;
      }
      opt->obj+=opt->c[indx]*bound;
    }
    else
      for (j=0; j<aj->nn0; j++)
        aj->an[j]=-aj->an[j];
    break;
      
  default:
    break;
  }
} /* TreatBound */

double xtQx(symatx* q,
            double* x)
{
  int    i,j,k,n;
  double qval,rtmp,diag;
  array* qj;

  n=q->ncol;
  qval=0.0;
  diag=0.0;

  for (j=0; j<n; j++)
  {
    qj=q->cols+j;
    diag+=x[j]*q->diag[j]*x[j];
    rtmp =0.0;
    for (k=0; k<qj->nn0; k++)
    {
      i=qj->ja[k];
      rtmp+=qj->an[k]*x[i];
    }
    qval+=x[j]*rtmp;
  }
  qval+=qval;
  qval+=diag;
  
  return qval;
} /* xtQx */

void QxPlusy(qpdat*  qt,
             symatx* q,
             double* x,
             double* y)
{
  int    i,j,k,n;
  array* qj;
  
  n=q->ncol;
  for (j=0; j<n; j++)
  {
    qj=q->cols+j;
    y[j]+=q->diag[j]*x[j];
    for (k=0; k<qj->nn0; k++)
    {
      i=qj->ja[k];
      y[j]+=qj->an[k]*x[i];
      y[i]+=qj->an[k]*x[j];
    }
  }
} /* QxPlusy */

static void PrintDim(qpdat* opt,
                     qppar* par)
{
  int  j,n,m,nz,*ibl,*ibu,*ibf;
  char sl[_L_SIZE],sr[_L_SIZE];
    
  m =opt->m_intRow;
  n =opt->m_intCol;
  nz=opt->m_intAnz;
  
  opt->nbup=0;
  opt->nblo=0;
  opt->nbfr=0;
  
  for (j=0; j<n; j++)
  {
    if (opt->u[j]<=par->q_dblMax)
      opt->nbup++;
    if (opt->l[j]>=-par->q_dblMax)
      opt->nblo++;
    if (opt->l[j]<-par->q_dblMax&&opt->u[j]>par->q_dblMax)
      opt->nbfr++;
  }
  
  ibl=opt->bid;
  ibf=opt->bid+opt->nblo;
  ibu=opt->bid+(opt->nblo-opt->nbup);
  
  for (j=0; j<n; j++)
  {
    if (opt->l[j]<-par->q_dblMax)
    {
      (*ibf)=j;
      ibf   ++;
    }
    else if (opt->u[j]<=par->q_dblMax)
    {
      (*ibu)=j;
      ibu   ++;
    }
    else
    {
      (*ibl)=j;
      ibl   ++;
    }
  }
  
  dFree(&opt->r);
  dFree(&opt->l);
  
  opt->l=opt->u;
  opt->u=dAlloc(opt->nbup);
  ibu   =opt->bid+(opt->nblo-opt->nbup);
  for (j=0; j<opt->nbup; j++)
    opt->u[j]=opt->l[ibu[j]];
  dFree(&opt->l);
  
  fprintf(fout,"Final Problem Dimension\n");
  fprintf(fout,"-----------------------\n");
  
  sprintf(sl,"   number of constraints");
  sprintf(sr,"%d",opt->m_intRow);
  LeftDots(sl,sr);
  
  sprintf(sl,"   number of variables");
  sprintf(sr,"%d",opt->m_intCol);
  LeftDots(sl,sr);
  
  sprintf(sl,"   number of upper bounds");
  sprintf(sr,"%d",opt->nbup);
  LeftDots(sl,sr);
  
  sprintf(sl,"   number of free variables");
  sprintf(sr,"%d",opt->nbfr);
  LeftDots(sl,sr);
    
  sprintf(sl,"   number of nonzeros in matrix A");
  sprintf(sr,"%d",opt->m_intAnz);
  LeftDots(sl,sr);
    
  sprintf(sl,"   number of nonzeros in matrix Q");
  sprintf(sr,"%d",2*opt->m_intQnz+opt->ndia);
  LeftDots(sl,sr);
  
  sprintf(sl,"   density of matrix A");
  sprintf(sr,"%.3f%s",
             100.0*((double)opt->m_intAnz/(double)n/(double)m),
             "%");
  LeftDots(sl,sr);
  
  sprintf(sl,"   density of matrix Q");
  sprintf(sr,"%.3f%s",
             100.0*((double)(2.0*opt->m_intQnz+opt->ndia)
             /(double)n/(double)n),"%");
  LeftDots(sl,sr);

  printf("\n");
  fprintf(fout,"\n");
} /* PrintDim */

static void InitMtx(matrix *a)
{
  int i;
  
  for (i=0; i<a->nrow; i++)
    a->ia[i].nn0=0;
} /* InitMtx */

static void ScanMtx(matrix *a,
                    matrix *at)
{
  int i,j,k,nrow=a->nrow,nn0;
  
  for (i=0; i<nrow; i++) {
    nn0=a->ia[i].nn0;
    for (k=0; k<nn0; k++) {
      j=a->ia[i].ja[k];
      at->ia[j].nn0++;
    }
  }
} /* ScanMtx */

static void makeMtxrw(matrix *a,
                      int    mems)
{
  if (a->ia) {
    iFree(&a->ia->ja);
    dFree(&a->ia->an);
    a->ia->ja=iAlloc(mems);
    a->ia->an=dAlloc(mems);
  }
} /* makeMtxrw */

int mTrans(int      ncol,
           matrix*  a,
           int      pid,
           int*     p,
           int      vid,
           matrix** at)
{
  int    i,j,k,t,nn0;
  matrix *r;
  array  *ai,*aj;
  
  nn0=0;
  for (i=0; i<a->nrow; i++)
    nn0+=a->ia[i].nn0;
  
  r=*at;
  
  if (r)
  {
    if (!vid&&nn0>r->mems)
      return false;
      
    if (ncol>r->mnrs)
      return false;
      
    r->nrow=ncol;
  }

  else
  {
    r=mAlloc(ncol,0);
    makeMtxrw(r,nn0);
  }
  
  if (vid)
    for (j=0; j<ncol; j++)
      r->ia[j].nn0=0;
  
  else
  {
    InitMtx(r);
    ScanMtx(a,r);
    
    nn0=0;
    for (j=0; j<ncol; j++)
    {
      r->ia[j].ja=r->ia->ja+nn0;
      r->ia[j].an=r->ia->an+nn0;
      nn0+=r->ia[j].nn0;
      r->ia[j].nn0=0;
    }
    
    r->mems=nn0;
    
    if (nn0>r->mems)
    {
      printf("\n\n Exit -- 140: out of internal space.\n\n");
      exit(0);
    }
  }
  
  for (k=0; k<a->nrow; k++)
  {
    i=k;
    if (p)
    {
      ai=a->ia+p[k];
      if (!pid) i=p[k];
    }
    else ai=a->ia+k;
    
    for (t=0; t<ai->nn0; t++)
	{
      j  =ai->ja[t];
      aj =r->ia+j;
      nn0=aj->nn0;
      
      aj->ja[nn0]=i;
      aj->an[nn0]=ai->an[t];
      aj->nn0++;
    }
  }
  
  *at=r;
  
  return true;
} /* mTrans */

int PreProc(qppar* par,
            qpdat* opt)
{
  int     i,j,k,n,nz;
  double  lj,uj,*r;
  char    key;
  array   *aj,*qj;
  symatx* qr;
  FILE    *fp;
  
  printf(" BEGIN preprocess...\n");
  
  n =opt->m_colUse+opt->m_intNsk;
  nz=opt->m_intAnz+opt->m_intNsk;
  
  opt->at  =mAlloc(n,nz);
  opt->u   =dAlloc(n);
  opt->l   =dAlloc(n);
  opt->c   =dAlloc(n);
  opt->bid =iAlloc(n);
  opt->imem=iAlloc(opt->m_intRow);
  opt->rmem=dAlloc(n);
  r        =opt->rmem;
  
  DataSwap(opt);
  
  n  =opt->m_colUse;
  
  fp=fopen("coplfile.bin","ab");
  if (!fp)
  {
    ShutDown(NOT_DSKSPC);
    exit(0);
  }
  
  fwrite(opt->c,sizeof(double),n,fp);

  for (j=0; j<n; j++)
  {
    lj=opt->l[j];
    uj=opt->u[j];
    aj=opt->at->ia+j;
    
    fwrite(&aj->nn0,sizeof(int),1,fp);
    fwrite(aj->ja,sizeof(int),aj->nn0,fp);
    fwrite(aj->an,sizeof(double),aj->nn0,fp);

    if (lj>uj)
    {
      ShutDown(INF_PROB);
      return (INF_PROB);
    }
    
    if (fabs(lj)<par->q_dblAij)
    {
      k=0;
      fwrite(&k,sizeof(int),1,fp);
      fwrite(&lj,sizeof(double),1,fp);
    }
    else if (lj>=-par->q_dblMax)
    {
      k=1;
      fwrite(&k,sizeof(int),1,fp);
      fwrite(&lj,sizeof(double),1,fp);
      TreatBound(lj,j,par,opt,LOW);
    }
    else if (uj<=par->q_dblMax)
    {
      k=2;
      fwrite(&k,sizeof(int),1,fp);
      fwrite(&uj,sizeof(double),1,fp);
      TreatBound(uj,j,par,opt,MIS);
    }
    else {
      k=3;
      fwrite(&k,sizeof(int),1,fp);
      fwrite(&lj,sizeof(double),1,fp);
    }
  }
  fclose(fp);

  for (j=0; j<n; j++)
  {
    if (opt->l[j]<-par->q_dblMax)
      r[j]=0.0;
    else
      r[j]=opt->l[j];
  }
  opt->obj+=0.5*xtQx(opt->qh,r);
  QxPlusy(opt,opt->qh,r,opt->c);
  dFree(&opt->rmem);

  nz=opt->m_intAnz;
  for (i=0; i<opt->m_intRow; i++)
  {
    key=opt->rowtype[i];
    if (opt->imem[i]==0)
    {
      if (fabs(opt->b[i])>=par->q_dblAij)
      {
        if (key=='E')
        {
          ShutDown(INF_PROB);
          return INF_PROB;
        }
        else if (key=='L'&&opt->b[i]<0)
        {
          ShutDown(INF_PROB);
          return INF_PROB;
        }
        else if (key=='G'&&opt->b[i]>0)
        {
          ShutDown(INF_PROB);
          return INF_PROB;
        }
        else if (key=='R'&&(opt->b[i]<0||opt->r[i]>0))
        {
          ShutDown(INF_PROB);
          return INF_PROB;
        }
      }
    }
    
    if (key=='E')
      continue;

    j           =n;
    aj          =opt->at->ia+j;
    aj->nn0     =1;
    aj->ja      =opt->at->ia->ja+nz;
    aj->an      =opt->at->ia->an+nz;
    opt->c[j]   =0.0;
    opt->l[j]   =0.0;
    opt->u[j]   =_P_MAXM;
    aj->ja[0]   =i;
    n           ++;
    nz          ++;
    
    switch(opt->rowtype[i])
    {
    case 'L':
      aj->an[0]=1.0;
      break;
    case 'G':
      aj->an[0]=-1.0;
      break;
    case 'R':
      aj->an[0]=-1.0;
      opt->u[j]=opt->b[i]-opt->r[i];
      opt->b[i]=opt->r[i];
      break;
    default:
      break;
    }
  }
  printf(" END preprocess\n");
  
  opt->m_intAnz=nz;
  opt->m_intCol=n;
  PrintDim(opt,par);
  iFree(&opt->imem);
  
  nz=0;
  for (j=0; j<opt->qh->ncol; j++)
    nz+=opt->qh->cols[j].nn0;
    
  qr=sAlloc(n,nz);
  nz=0;
  for (j=0; j<opt->qh->ncol; j++)
  {
    aj         =qr->cols+j;
    qj         =opt->qh->cols+j;
    aj->ja     =qr->cols->ja+nz;
    aj->an     =qr->cols->an+nz;
    aj->nn0    =qj->nn0;
    qr->diag[j]=opt->qh->diag[j];
    
    for (k=0; k<qj->nn0; k++)
    {
      aj->ja[k]=qj->ja[k];
      aj->an[k]=qj->an[k];
    }
    nz+=k;
  }
  for (j=opt->qh->ncol; j<n; j++)
  {
    aj     =qr->cols+j;
    aj->ja =qr->cols->ja+nz;
    aj->an =qr->cols->an+nz;
    aj->nn0=0;
    qr->diag[j]=0.0;
  }
  sFree(&(opt->qh));
  opt->qh=qr;
  
  mTrans(opt->m_intRow,opt->at,false,NULL,false,&opt->a);
  mTrans(opt->m_intCol,opt->a, false,NULL,false,&opt->at);
  mFree(&opt->a);
  
  return OPT_OK;
} /* PreProc */
