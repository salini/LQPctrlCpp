#include "QPdef.h"

static void SetQpdat(qpdat* opt)
{
  int m,nn,nblo,nbup,nbfr;

  m   =opt->m_intRow;
  nn  =opt->m_intCol;
  nblo=opt->nblo;
  nbup=opt->nbup;
  nbfr=opt->nbfr;

  opt->rp  =dAlloc(m);
  opt->rd  =dAlloc(nn);
  opt->rl  =dAlloc(nblo);
  opt->ru  =dAlloc(nbup);
  opt->rr  =dAlloc(m);
  opt->ry  =dAlloc(m);
  opt->rf  =dAlloc(nbfr);
  opt->rh  =dAlloc(nbfr);
  
  opt->rz  =dAlloc(nn);
  opt->rw  =dAlloc(m);
  opt->rs  =dAlloc(nbup);
  opt->rq  =dAlloc(m);
  opt->rk  =dAlloc(nbfr);
  opt->qx  =dAlloc(nn);

  opt->diag=dAlloc(nn);
  opt->diau=dAlloc(nbup);
  opt->diaf=dAlloc(nbfr);
  opt->dial=dAlloc(m);

  opt->x   =dAlloc(nn);
  opt->y   =dAlloc(m);
  opt->g   =dAlloc(nblo);
  opt->z   =dAlloc(nblo);
  opt->t   =dAlloc(nbup);
  opt->s   =dAlloc(nbup);
  opt->v   =dAlloc(m);
  opt->w   =dAlloc(m);
  opt->p   =dAlloc(m);
  opt->q   =dAlloc(m);
  opt->f   =dAlloc(nbfr);
  opt->h   =dAlloc(nbfr);
  opt->k   =dAlloc(nbfr);
  opt->l   =dAlloc(nbfr);

  opt->dx  =dAlloc(nn);
  opt->dy  =dAlloc(m);
  opt->dg  =dAlloc(nblo);
  opt->dz  =dAlloc(nblo);
  opt->dt  =dAlloc(nbup);
  opt->ds  =dAlloc(nbup);
  opt->dv  =dAlloc(m);
  opt->dw  =dAlloc(m);
  opt->dp  =dAlloc(m);
  opt->dq  =dAlloc(m);
  opt->df  =dAlloc(nbfr);
  opt->dh  =dAlloc(nbfr);
  opt->dk  =dAlloc(nbfr);
  opt->dl  =dAlloc(nbfr);
  
  opt->imem=iAlloc(2*(m+nn));
  opt->rmem=dAlloc(2*(m+nn));
} /* SetQpdat */

static double dDots(int     n,
                    double* x,
                    double* y,
                    int*    s)
{
  int    i;
  double r;
  
  r=0.0;
  if (s)
  {
    for (i=0; i<n; i++)
      r+=x[s[i]]*y[s[i]];
  }
  else
  {
    for (i=0; i<n; i++)
      r+=x[i]*y[i];
  }
  return r;
} /* dDots */

static void ScalVect(int     n,
                     double  alfa,
                     double  *x)
{
  int     i;

  if (n<=0 || alfa==1.0)
    return;

  if (alfa==0.0) {
    memset(x,0,sizeof(double)*n);
    return;
  }
  
  else if (alfa==-1.0) {
    for(i=0; i<n; i++)
      x[i]=-x[i];
  }

  else {
    for(i=0; i<n; i++)
      x[i]*=alfa;
  }

} /* ScalVect */

static void mTimesv(int    tran,
                    int    nrow,
                    int    ncol,
                    double alfa,
                    matrix *a,
                    double *x,
                    double beta,
                    double *y)
{
  int      xsz,ysz,i;
  
  if (nrow<0||ncol<0) {
    ShutDown(SZE_ERROR);
    exit(0);
  }
  
  if (!tran) {
    xsz=ncol;
    ysz=nrow;

    ScalVect(ysz,beta,y);

    if (alfa==0.0)
      return;

    for(i=0; i<xsz; i++)
      if (x[i])
        setArray(alfa*x[i],a->ia+i,y);
  }
  else {
    xsz=nrow;
    ysz=ncol;

    ScalVect(ysz,beta,y);

    if (alfa==0.0)
      return;

    for(i=0; i<ysz; i++)
      y[i]+=alfa*svDot(a->ia+i,x);
  }
} /* mTimesv */

static void setXYind(int     nnz,
                     double* y,
                     double* x,
                     int*    s)
{
  int i;
  
  for(i=0; i<nnz; ++i) {
    x[i]=y[s[i]];
    y[s[i]]=0.0;
  }
} /* setXYind */

static void setColi(chol*   cl,
                    int     i,
                    double* ai)
{
  setXYind(cl->ujsze[i],ai,
           cl->uval+cl->uhead[i],
           cl->usub+cl->ujbeg[i]);
} /* setColi */

static void ClSetup(chol*   cl,
                    symatx* q,
                    matrix* at)
{
  int    h,j,k,i,m,n,nn,nnz,nrow,nil,next,
         *ifir,*inex,*icur,*invp;
  double *rw;
  array  *qi,*ai;
  
  n   =q->ncol;
  nn  =at->nrow;
  nrow=cl->nrow;
  m   =nrow-nn;
  invp=cl->invp;
  
  rw  =dAlloc(nrow);
  ifir=iAlloc(3*nrow);
  icur=ifir+nrow;
  inex=icur+nrow;
  
  nil=nrow;
  for (i=0; i<nrow; i++)
    ifir[i]=nil;
  
  for (i=0; i<n; i++)
  {
    icur[i]=0;
    qi  =q->cols+i;
    nnz =qi->nn0;
    
    if (nnz>0)
    {
      j=qi->ja[0];
      inex[i]=ifir[j];
      ifir[j]=i;
    }
    else
      inex[i]=nil;
  }
  
  for (i=0; i<n; i++)
  {
    memset(rw,0,nrow*sizeof(double));
    qi =q->cols+i;
    for (k=0; k<qi->nn0; k++)
      rw[invp[qi->ja[k]]]=-qi->an[k];
    
    for (j=ifir[i]; j!=nil; j=next)
    {
      next=inex[j];
      qi=q->cols+j;
      
      rw[invp[j]]=-qi->an[icur[j]];
      icur[j]++;

      if (icur[j]<qi->nn0)
      {
        k=qi->ja[icur[j]];
        inex[j]=ifir[k];
        ifir[k]=j;
      }
    }

    ai=at->ia+i;
    for (k=0; k<ai->nn0; k++)
      rw[invp[ai->ja[k]+nn]]=ai->an[k];
    
    setColi(cl,invp[i],rw);
  }
  
  for (i=n; i<nn; i++)
  {
    memset(rw,0,nrow*sizeof(double));
    ai=at->ia+i;
    for (k=0; k<ai->nn0; k++)
      rw[invp[ai->ja[k]+nn]]=ai->an[k];
    setColi(cl,invp[i],rw);

    setColi(cl,invp[i],rw);
  }
  
  for (i=0; i<nrow; i++)
    ifir[i]=nil;
  
  for (i=0; i<nn; i++)
  {
    icur[i]=0;
    ai  =at->ia+i;
    nnz =ai->nn0;
    
    if (nnz>0)
    {
      j=ai->ja[0];
      inex[i]=ifir[j];
      ifir[j]=i;
    }
    else
      inex[i]=nil;
  }
  
  for (i=0; i<m; i++)
  {
    memset(rw,0,nrow*sizeof(double));
    for (j=ifir[i]; j!=nil; j=next)
    {
      next=inex[j];
      ai  =at->ia+j;
      h   =icur[j];
      rw[invp[j]]=ai->an[h];
      icur[j]++;

      if (icur[j]<ai->nn0)
      {
        k=ai->ja[icur[j]];
        inex[j]=ifir[k];
        ifir[k]=j;
      }
    }

    setColi(cl,invp[nn+i],rw);
  }
  
  dFree(&rw);
  iFree(&ifir);

} /* ClSetup */

static void InitQpdat(qpdat* opt)
{
  int j,k,m,nl,nu,nf,nn,nm,*ibl,*ibu,*ibf,*invp;
  double rtmp,*rhs,*sol;
  
  m =opt->m_intRow;
  nl=opt->nblo;
  nu=opt->nbup;
  nf=opt->nbfr;
  nn=nl+nf;
  nm=nn+m;
  
  ibl =opt->bid;
  ibu =ibl+(nl-nu);
  ibf =ibl+nl;
  invp=opt->cl->invp;
  
  rhs=opt->rmem;
  sol=opt->rmem+(nm);
  
  ClSetup(opt->cl,opt->qh,opt->at);
  for (j=0; j<nn; j++)
    opt->cl->diag[invp[j]]=-opt->qh->diag[j]-1.0;
  for (j=nn; j<nm; j++)
    opt->cl->diag[invp[j]]=1.0;
  
  if (OPT_OK!=CholFact(opt->cl,opt->imem,opt->rmem,false))
  {
    ShutDown(FAC_ERROR);
    exit(0);
  }
  
  for (j=0; j<nn; j++)
    rhs[j]=opt->c[j];
  for (j=0; j<m; j++)
    rhs[j+nn]=opt->b[j];
    
  TriSolv(opt->cl,rhs,sol);
  
  for (k=0; k<nl; k++)
  {
    j        =ibl[k];
    opt->x[j]=sol[j];
    opt->g[k]=max(fabs(opt->x[j]),pa->q_dblIni);
    opt->z[k]=opt->g[k];
  }
  for (k=0; k<nu; k++)
  {
    j        =ibu[k];
    opt->t[k]=max(fabs(opt->u[k]-opt->x[j]),pa->q_dblIni);
    opt->s[k]=max(fabs(opt->x[j]),pa->q_dblIni);
  }
  for (k=0; k<nf; k++)
  {
    j=ibf[k];
    rtmp=sol[j];
    opt->x[j]=rtmp;
    rtmp=max(fabs(rtmp),pa->q_dblIni);
    opt->f[k]=rtmp;
    opt->h[k]=rtmp;
    opt->k[k]=rtmp;
    opt->l[k]=rtmp;
  }
  
  for (j=0; j<m; j++)
  {
    opt->y[j]=sol[j+nn];
    opt->v[j]=max(fabs(opt->y[j]),pa->q_dblIni);
    opt->w[j]=opt->v[j];
    opt->p[j]=opt->v[j];
    opt->q[j]=opt->v[j];
  }    
} /* InitQpdat */

static void ShowItrTitle(void)
{
#ifdef UNIXMACHINE
  fprintf(fout,"%-6s  %-16s  %-16s  %-9s  %-9s  %-9s\n",
               "ITER","     POBJ","     DOBJ","  RGAP","  PINF","  DIFN");
  printf(" %-5s  %-16s  %-16s  %-9s  %-9s  %-9s\n",
         "ITER","       POBJ","       DOBJ","   RGAP","   PINF","   DIFN");
#else
  fprintf(fout,"%-6s  %-16s  %-16s  %-9s  %-9s  %-9s\n",
               "ITER"," POBJ"," DOBJ","RGAP","PINF","DIFN");
  printf(" %-5s  %-16s  %-16s  %-9s  %-9s  %-9s\n",
         "ITER"," POBJ"," DOBJ","RGAP","PINF","DIFN");
#endif
  fprintf(fout,"--------------------------------------"
               "-------------------------------------\n");
  printf(" -------------------------------------"
         "-------------------------------------\n");
} /* ShowItrTitle */

static void ShowItrMsg(qpdat* opt)
{
  fprintf(fout,"%-6d  %+16.8e  %+16.8e  %9.2e  %9.2e  %9.2e\n",
               opt->iter,opt->pobj,opt->dobj,
               opt->rgap,opt->pinf,opt->dinf);
  printf(" %-5d  %+16.8e  %+16.8e  %9.2e  %9.2e  %9.2e\n",
          opt->iter,opt->pobj,opt->dobj,
          opt->rgap,opt->pinf,opt->dinf);
} /* ShowItrTitle */

static void GetVals(qpdat* opt)
{
  int    j,k,m,nl,nu,nf,nn,nm,*ibl,*ibu,*ibf;
  double rtmp;
  
  m =opt->m_intRow;
  nl=opt->nblo;
  nu=opt->nbup;
  nf=opt->nbfr;
  nn=nl+nf;
  nm=nn+m;
  
  ibl=opt->bid;
  ibu=ibl+(nl-nu);
  ibf=ibl+nl;
  
  rtmp     =0.0;
  opt->pobj=0.0;
  opt->dobj=0.0;
  opt->cgap=0.0;
  
  for (j=0; j<m; j++)
  {
    opt->rp[j]=opt->b[j]+opt->w[j];
    opt->rr[j]=-opt->w[j]-opt->p[j];
    opt->ry[j]=opt->y[j]+opt->q[j]-opt->v[j];
    
    opt->dial[j]=1.0/(opt->v[j]/opt->w[j]+opt->q[j]/opt->p[j]);
    
    opt->dobj+=opt->b[j]*opt->y[j];
    opt->cgap+=opt->v[j]*opt->w[j]+opt->p[j]*opt->q[j];
  }
  mTimesv(false,m,nn,-1.0,opt->at,opt->x,1.0,opt->rp);
    
  memset(opt->qx,0,nn*sizeof(double));
  QxPlusy(opt,opt->qh,opt->x,opt->qx);
  
  for (k=0; k<nl; k++)
  {
    j         =ibl[k];
    opt->rl[k]=-opt->x[j]+opt->g[k];
    opt->rd[j]=opt->c[j]-opt->z[k]+opt->qx[j];
    
    opt->diag[j]=opt->z[k]/opt->g[k];
    
    rtmp     +=opt->x[j]*opt->qx[j];
    opt->pobj+=opt->c[j]*opt->x[j];
    opt->cgap+=opt->z[k]*opt->g[k];
  }  
  
  for (k=0; k<nu; k++)
  {
    j         =ibu[k];
    opt->ru[k] =opt->u[k]-opt->x[j]-opt->t[k];
    opt->rd[j]+=opt->s[k];
    
    opt->diau[k] =opt->s[k]/opt->t[k];
    opt->diag[j]+=opt->diau[k];
    
    opt->dobj -=opt->u[k]*opt->s[k];
    opt->cgap +=opt->s[k]*opt->t[k];
  }

  for (k=0; k<nf; k++)
  {
    j            =ibf[k];
    opt->rd[j]   =opt->c[j]+opt->qx[j]-opt->k[k];
    opt->rf[k]   =-opt->k[k]-opt->l[k];
    opt->rh[k]   =opt->f[k]-opt->x[j]-opt->h[k];
    
    opt->diaf[k] =1.0/(opt->f[k]/opt->k[k]+opt->h[k]/opt->l[k]);
    opt->diag[j] =opt->diaf[k];
    
    rtmp        +=opt->x[j]*opt->qx[j];
    opt->pobj   +=opt->c[j]*opt->x[j];
    opt->cgap   +=opt->f[k]*opt->k[k]+opt->h[k]*opt->l[k];
  }
  
  mTimesv(true,m,nn,-1.0,opt->at,opt->y,1.0,opt->rd);
  
  opt->pinf = sqrt(dDots(m, opt->rp,opt->rp,NULL)
                  +dDots(nu,opt->ru,opt->ru,NULL)
                  +dDots(m, opt->rr,opt->rr,NULL)
                  +dDots(nl,opt->rl,opt->rl,NULL)
                  +dDots(nf,opt->rh,opt->rh,NULL))/
             (sqrt(dDots(m, opt->b,opt->b,NULL))+1.0);
  opt->dinf = sqrt(dDots(nn,opt->rd,opt->rd,NULL)
                  +dDots(nf,opt->rf,opt->rf,NULL)
                  +dDots(m, opt->ry,opt->ry,NULL))/
             (sqrt(dDots(nn,opt->c,opt->c,NULL))+1.0);
  opt->pobj+=opt->obj+0.5*rtmp;
  opt->dobj+=opt->obj-0.5*rtmp;
  opt->rgap =fabs(opt->pobj-opt->dobj)/(1.0+fabs(opt->dobj));
  opt->mu   =opt->cgap/(double)(2.0*(nf+m)+nl+nu);  
} /* GetVals */

static void GetDirect(qpdat*  opt,
                      chol*   cl,
                      int*    ibuf,
                      double* rbuf,
                      double* dz,
                      double* dw,
                      double* dq,
                      double* ds,
                      double* dx,
                      double* dy,
                      double* dg,
                      double* dv,
                      double* dp,
                      double* dt,
                      double* df,
                      double* dh,
                      double* dk,
                      double* dl)
{
  int    j,k,m,nn,nm,nl,nu,nf,*ibl,*ibu,*ibf;
  double rtmp,rho,*rhs,*sol;
  
  m =opt->m_intRow;
  nl=opt->nblo;
  nu=opt->nbup;
  nf=opt->nbfr;
  nn=nl+nf;
  nm=nn+m;
  
  ibl=opt->bid;
  ibu=ibl+(nl-nu);
  ibf=ibl+nl;
  
  rhs=rbuf;
  sol=rhs+nm;
  
  for (k=0; k<nl; k++)
  {
    j      =ibl[k];
    rhs[j] =opt->rd[j]-(opt->z[k]/opt->g[k])*dz[k];
  }
  for (k=0; k<nu; k++)
  {
    j      =ibu[k];
    rhs[j]-=(opt->s[k]/opt->t[k])*ds[k];
  }
  for (k=0; k<nf; k++)
  {
    j      =ibf[k];
    rtmp   =(opt->dk[k]+(opt->h[k]/opt->l[k])*opt->dl[k])
            /(opt->f[k]/opt->k[k]+opt->h[k]/opt->l[k]);
    rhs[j] =opt->rd[j]-rtmp;
  }
  for (j=0; j<m; j++)
  {
    rtmp=opt->dial[j]*(dw[j]-(opt->q[j]/opt->p[j])*dq[j]);
    rhs[nn+j]=opt->rp[j]-rtmp;
    dw[j]    =-rtmp;
  }
    
  TriSolv(opt->cl,rhs,sol);
  
  for (j=0; j<nn; j++)
    dx[j]=sol[j];
  for (j=0; j<m; j++)
  {
    dy[j] =sol[j+nn];
    dw[j]-=opt->dial[j]*dy[j];
    dq[j] =(opt->q[j]/opt->p[j])*(dw[j]-dq[j]);
    dv[j] =(opt->v[j]/opt->w[j])*(opt->rw[j]-dw[j]);
    dp[j] =(opt->p[j]/opt->q[j])*(opt->rq[j]-dq[j]);
  }
  for (k=0; k<nl; k++)
  {
    j    =ibl[k];
    dz[k]=(opt->z[k]/opt->g[k])*(dz[k]-dx[j]);
    dg[k]=(opt->g[k]/opt->z[k])*(opt->rz[k]-dz[k]);
  }
  for (k=0; k<nu; k++)
  {
    j    =ibu[k];
    ds[k]=(opt->s[k]/opt->t[k])*(dx[j]-ds[k]);
    dt[k]=(opt->t[k]/opt->s[k])*(opt->rs[k]-ds[k]);
  }
  for (k=0; k<nf; k++)
  {
    j    =ibf[k];
    rtmp =opt->f[k]/opt->k[k];
    rho  =opt->h[k]/opt->l[k];
    dh[k]=-(rtmp*dl[k]-dk[k]+dx[j])*rho/(rtmp+rho);
    dk[k]=(dk[k]-dx[j]-dh[k])/rtmp;
    df[k]=(opt->rk[k]-dk[k])*rtmp;
    dl[k]=opt->rf[k]-dh[k]/rho;
  }
} /* GetDirect */

static double GetStep(qpdat*  opt,
                      double* g,
                      double* dg,
                      double* w,
                      double* dw,
                      double* t,
                      double* dt,
                      double* q,
                      double* dq,
                      double* f,
                      double* df,
                      double* h,
                      double* dh)
{
  int    k,m,nl,nu,nf;
  double rtmp,xdx;
  
  m =opt->m_intRow;
  nl=opt->nblo;
  nu=opt->nbup;
  nf=opt->nbfr;
  
  rtmp=1.0;
  for (k=0; k<nl; k++)
  {
    if (dg[k]<-1.0e-20)
    {
      xdx =-dg[k]/g[k];
      rtmp=max(rtmp,xdx/pa->q_dblStp);
    }
  }
  for (k=0; k<nu; k++)
  {
    if (dt[k]<-1.0e-20)
    {
      xdx =-dt[k]/t[k];
      rtmp=max(rtmp,xdx/pa->q_dblStp);
    }
  }
  for (k=0; k<nf; k++)
  {
    if (df[k]<-1.0e-20)
    {
      xdx =-df[k]/f[k];
      rtmp=max(rtmp,xdx/pa->q_dblStp);
    }
    if (dh[k]<-1.0e-20)
    {
      xdx =-dh[k]/h[k];
      rtmp=max(rtmp,xdx/pa->q_dblStp);
    }
  }
  for (k=0; k<m; k++)
  {
    if (dw[k]<-1.0e-20)
    {
      xdx =-dw[k]/w[k];
      rtmp=max(rtmp,xdx/pa->q_dblStp);
    }
    if (dq[k]<-1.0e-20)
    {
      xdx =-dq[k]/q[k];
      rtmp=max(rtmp,xdx/pa->q_dblStp);
    }
  }
  
  return rtmp;
} /* GetStep */

static int CheckNeighbor(qpdat* opt,
                         qppar* par,
                         double apd)
{
  int    k,m,nl,nu,nf;
  double xzj,cgap,rtmp;
  
  m   =opt->m_intRow;
  nl  =opt->nblo;
  nu  =opt->nbup;
  nf  =opt->nbfr;

  cgap=0.0;
  xzj =_P_MAXM;
  
  for (k=0; k<m; k++)
  {
    rtmp =(opt->v[k]+opt->dv[k]/apd)*(opt->w[k]+opt->dw[k]/apd);
    cgap+=rtmp;
    xzj  =min(xzj,rtmp);
    rtmp =(opt->p[k]+opt->dp[k]/apd)*(opt->q[k]+opt->dq[k]/apd);
    cgap+=rtmp;
    xzj  =min(xzj,rtmp);
  }
  for (k=0; k<nl; k++)
  {
    rtmp =(opt->g[k]+opt->dg[k]/apd)*(opt->z[k]+opt->dz[k]/apd);
    cgap+=rtmp;
    xzj  =min(xzj,rtmp);
  }
  for (k=0; k<nu; k++)
  {
    rtmp =(opt->s[k]+opt->ds[k]/apd)*(opt->t[k]+opt->dt[k]/apd);
    cgap+=rtmp;
    xzj  =min(xzj,rtmp);
  }
  for (k=0; k<nf; k++)
  {
    rtmp =(opt->f[k]+opt->df[k]/apd)*(opt->k[k]+opt->dk[k]/apd);
    cgap+=rtmp;
    xzj  =min(xzj,rtmp);
    rtmp =(opt->h[k]+opt->dh[k]/apd)*(opt->l[k]+opt->dl[k]/apd);
    cgap+=rtmp;
    xzj  =min(xzj,rtmp);
  }
  
  return (xzj>par->q_dblRad*cgap);
} /* CheckNeighbor */

static void SaveSol(qpdat* opt)
{
  int    h,i,j,k,m,n,nl,nu,nf,nz,*ib,*ibl,*ibu,*ibf;
  double *l,*r,*u,*z,rtmp;
  FILE   *fp;
  char   cp,lb[20],ub[20],vl[20],
         **rowname,**colname;
  array  *aj;

  nl=opt->nblo;
  nu=opt->nbup;
  nf=opt->nbfr;

  ibl=opt->bid;
  ibu=ibl+(nl-nu);
  ibf=ibl+nl;

  ib =opt->imem;
  l  =opt->dx;
  u  =opt->rd;
  z  =opt->rmem;
  r  =z+nl+nf;

  mFree(&opt->at);

  fp=fopen("coplfile.bin","rb");
  if (!fp)
  {
    ShutDown(NOT_SCHFIL);
    exit(0);
  }
  fread(&m,sizeof(int),1,fp);
  fread(&n,sizeof(int),1,fp);
  fread(&nz,sizeof(int),1,fp);

  memset(ib,0,n*sizeof(int));
  for (j=0; j<n; j++)
  {
    l[j]=0.0;
    u[j]=_P_MAXM;
    z[j]=0.0;
  }

  rowname=cPallo(m,0);
  colname=cPallo(n,0);
  opt->at=mAlloc(n,nz);

  for (i=0; i<m; i++)
  {
    fread(&k,sizeof(int),1,fp);
    rowname[i]=cAlloc(k+1);
    fread(rowname[i],sizeof(char),k,fp);
  }
  for (j=0; j<n; j++)
  {
    fread(&k,sizeof(int),1,fp);
    colname[j]=cAlloc(k+1);
    fread(colname[j],sizeof(char),k,fp);
  }
  fread(r,sizeof(double),m,fp);
  fread(opt->b,sizeof(double),m,fp);
  fread(opt->c,sizeof(double),n,fp);

  h=0;
  for (j=0; j<n; j++)
  {
    aj    =opt->at->ia+j;
    aj->ja=opt->at->ia->ja+h;
    aj->an=opt->at->ia->an+h;

    fread(&k,sizeof(int),1,fp);
    fread(aj->ja,sizeof(int),k,fp);
    fread(aj->an,sizeof(double),k,fp);
    aj->nn0 =k;
    h      +=k;

    fread(&i,sizeof(int),1,fp);
    fread(&rtmp,sizeof(double),1,fp);

    ib[j]=i;
    if (i==1)
    {
      l[j]      =rtmp;
      opt->x[j]+=rtmp;
    }
    else if (i==2)
    {
      l[j]     =rtmp;
      opt->x[j]=rtmp-opt->x[j];
    }
  }
  fclose(fp);

  for (k=0; k<nu; k++)
  {
    j   =ibu[k];
    u[j]=opt->u[k]+l[j];
  }
  for (k=0; k<nf; k++)
  {
    j   =ibf[k];
    l[j]=-_P_MAXM;
  }

  memset(opt->qx,0,(nl+nf)*sizeof(double));
  QxPlusy(opt,opt->qh,opt->x,opt->qx);

  /*
   * Ax -> rp, c+Q*x-A^T*y -> z
   */
  memset(opt->rp,0,m*sizeof(double));
  for (j=0; j<n; j++)
  {
    aj  =opt->at->ia+j;
    z[j]=opt->c[j]+opt->qx[j];
    for (k=0; k<aj->nn0; k++)
    {
      i          =aj->ja[k];
      opt->rp[i]+=aj->an[k]*opt->x[j];
      z[j]      -=aj->an[k]*opt->y[i];
    }
  }

  fp=fopen("coplfile.sol","w");
  if (!fp)
  {
    ShutDown(NOT_DSKSPC);
    exit(0);
  }

  fprintf(fp,"COPL_QP: Solution File\n");
  fprintf(fp,"======================\n\n");
  fprintf(fp,"ROWS Section\n");
  fprintf(fp,"------------\n");

  fprintf(fp,"%-5s %-8s %4s  %-12s  %-12s  %-12s  %-12s\n",
             "order","name","type","lower","active","upper","dual value");
  for (i=0; i<m; i++)
  {
    if (opt->rowtype[i]=='E')
    {
      sprintf(lb,"%+E",r[i]);
      sprintf(ub,"%+E",opt->b[i]);
    }
    else if (opt->rowtype[i]=='L')
    {
      sprintf(lb,"-infinity");
      sprintf(ub,"%+E",opt->b[i]);
    }
    else if (opt->rowtype[i]=='G')
    {
      sprintf(lb,"%+E",opt->b[i]);
      sprintf(ub,"+infinity");
    }
    else if (opt->rowtype[i]=='R')
    {
      sprintf(lb,"%+E",r[i]);
      sprintf(ub,"%+E",opt->b[i]);
    }
    sprintf(vl,"%+E",opt->rp[i]);

    fprintf(fp,"%-5d %-8s  %c   %-13s %-13s %-13s %+E\n",
               i+1,rowname[i],opt->rowtype[i],
               lb,vl,ub,opt->y[i]);
  }

  fprintf(fp,"\nCOLUMNS Section\n");
  fprintf(fp,"---------------\n");

  fprintf(fp,"%-5s %-8s %4s  %-12s  %-12s  %-12s  %-12s\n",
             "order","name","type","lower","active","upper","dual slack");
  for (j=0; j<n; j++)
  {
    i=0;
    k=0;
    if (l[j]<-0.999*_P_MAXM)
    {
      sprintf(lb,"-infinity");
      i++;
    }
    else
      sprintf(lb,"%+E",l[j]);

    if (u[j]>0.999*_P_MAXM)
    {
      sprintf(ub,"+infinity");
      k++;
    }
    else
      sprintf(ub,"%+E",u[j]);
    sprintf(vl,"%+E",opt->x[j]);

    if (!i&&!k)
      cp='B';
    else
    {
      if (i&&k)
        cp='F';
      else if (i)
        cp='U';
      else
        cp='L';
    }

    fprintf(fp,"%-5d %-8s  %c   %-13s %-13s %-13s %+E\n",
               j+1,colname[j],cp,lb,vl,ub,z[j]);
  }
  fprintf(fp,"\n----------------------------------- EOF ");
  fprintf(fp,"-----------------------------------\n");
  fclose(fp);
//----- Possible Error Found -----

  for (i=0; i<m; i++)
  {
    cFree (&(rowname[i]));
  }
  for (j=0; j<n; j++)
  {
    cFree(&(colname[j]));
  }
  cPfre(&rowname,m);
  cPfre(&colname,n);
} /* SaveSol */

int QpSolv(qpdat* opt,
           qppar* par)
{
  int    j,k,m,n,nn,nm,nl,nu,nf,np,nd,acc,
         *ibl,*ibu,*ibf,*invp;
  double rtmp,apd,ap,ad,pinf0,dinf0;

  m   =opt->m_intRow;
  n   =opt->m_colUse;
  nn  =opt->m_intCol;
  nm  =m+nn;
  nl  =opt->nblo;
  nu  =opt->nbup;
  nf  =opt->nbfr;

  ibl =opt->bid;
  ibu =ibl+(nl-nu);
  ibf =ibl+nl;
  invp=opt->cl->invp;

  if (n<100) par->q_dblIni=100.0;
  if (nf) par->q_dblStp=0.9;

  SetQpdat(opt);
  InitQpdat(opt);
  ShowItrTitle();

  np   =0;
  nd   =0;
  pinf0=_P_MAXM;
  dinf0=_P_MAXM;

  /*
   * main loop
   */
  for (opt->iter=0; opt->iter<par->q_intItr; opt->iter++)
  {
    GetVals(opt);
    ShowItrMsg(opt);
    
    if (opt->rgap<par->q_dblGap&&
        opt->pinf<par->q_dblPri&&
        opt->dinf<par->q_dblDua)
    {
      opt->m_intSta=OPT_FOUND;
      break;
    }
    else if (nd>5)
    {
      opt->m_intSta=UNB_PROB;
      break;
    }
    else if (np>5)
    {
      opt->m_intSta=INF_PROB;
      break;
    }
    else
    {
      if (opt->pinf<pinf0-1.0e-7)
      {
        pinf0=opt->pinf;
        np=0;
      }
      else
        np++;
      if (opt->dinf<dinf0-1.0e-7)
      {
        dinf0=opt->dinf;
        nd=0;
      }
      else
        nd++;
    }

    ClSetup(opt->cl,opt->qh,opt->at);
    for (j=0; j<nn; j++)
      opt->cl->diag[invp[j]]=-opt->qh->diag[j]-opt->diag[j];
    for (j=nn; j<nm; j++)
      opt->cl->diag[invp[j]]=opt->dial[j-nn];

    if (OPT_OK!=CholFact(opt->cl,opt->imem,opt->rmem,false))
    {
      opt->m_intSta=NUM_PROB;
      break;
    }

    /*
     * predictor step
     */
    for (k=0; k<nl; k++)
    {
      opt->dz[k]= opt->rl[k]-opt->g[k];
      opt->rz[k]=-opt->z[k];
    }
    for (j=0; j<m; j++)
    {
      opt->dw[j]=opt->ry[j]+opt->v[j];
      opt->dq[j]=opt->rr[j]+opt->p[j];
      opt->rw[j]=-opt->w[j];
      opt->rq[j]=-opt->q[j];
    }
    for (j=0; j<nu; j++)
    {
      opt->ds[j]=opt->ru[j]+opt->t[j];
      opt->rs[j]=-opt->s[j];
    }
    for (k=0; k<nf; k++)
    {
      opt->rk[k]=-opt->k[k];
      opt->rf[k]=-opt->l[k];
      opt->dk[k]= opt->rh[k]-opt->f[k];
      opt->dl[k]= opt->rf[k]+opt->k[k];
    }
    GetDirect(opt,opt->cl,opt->imem,opt->rmem,
              opt->dz,opt->dw,opt->dq,opt->ds,
              opt->dx,opt->dy,
              opt->dg,opt->dv,opt->dp,opt->dt,
              opt->df,opt->dh,opt->dk,opt->dl);

    ap=GetStep(opt,opt->g,opt->dg,opt->w,opt->dw,
                   opt->t,opt->dt,opt->p,opt->dp,
                   opt->f,opt->df,opt->h,opt->dh);
    ad=GetStep(opt,opt->z,opt->dz,opt->v,opt->dv,
                   opt->s,opt->ds,opt->q,opt->dq,
                   opt->k,opt->dk,opt->l,opt->dl);
    apd=max(ad,ap);

    rtmp=(apd-1.0)/(apd+10);
    apd=rtmp*rtmp*opt->mu;

    /*
     * corrector step
     */
    for (j=0; j<nl; j++)
    {
      opt->rz[j]=apd/opt->g[j]-opt->z[j]-opt->dg[j]*opt->dz[j]/opt->g[j];
      opt->dz[j]=opt->rl[j]+opt->g[j]*opt->rz[j]/opt->z[j];
    }
    for (j=0; j<nu; j++)
    {
      opt->rs[j]=apd/opt->t[j]-opt->s[j]-opt->dt[j]*opt->ds[j]/opt->t[j];
      opt->ds[j]=opt->ru[j]-opt->t[j]*opt->rs[j]/opt->s[j];
    }
    for (j=0; j<nf; j++)
    {
      opt->rk[j]=apd/opt->f[j]-opt->k[j]-(opt->df[j]/opt->f[j])*opt->dk[j];
      opt->rf[j]=apd/opt->h[j]-opt->l[j]-(opt->dh[j]/opt->h[j])*opt->dl[j];
      opt->dk[j]=opt->rh[j]+(opt->f[j]/opt->k[j])*opt->rk[j];
      opt->dl[j]=opt->rf[j]-opt->rk[j];
    }
    for (j=0; j<m; j++)
    {
      opt->rw[j]=apd/opt->v[j]-opt->w[j]-opt->dv[j]*opt->dw[j]/opt->v[j];
      opt->rq[j]=apd/opt->p[j]-opt->q[j]-opt->dp[j]*opt->dq[j]/opt->p[j];
      opt->dw[j]=opt->ry[j]-opt->v[j]*opt->rw[j]/opt->w[j];
      opt->dq[j]=opt->rr[j]-opt->p[j]*opt->rq[j]/opt->q[j];
    }

    GetDirect(opt,opt->cl,opt->imem,opt->rmem,
              opt->dz,opt->dw,opt->dq,opt->ds,
              opt->dx,opt->dy,
              opt->dg,opt->dv,opt->dp,opt->dt,
              opt->df,opt->dh,opt->dk,opt->dl);

    ap=GetStep(opt,opt->g,opt->dg,opt->w,opt->dw,
                   opt->t,opt->dt,opt->p,opt->dp,
                   opt->f,opt->df,opt->h,opt->dh);
    ad=GetStep(opt,opt->z,opt->dz,opt->v,opt->dv,
                   opt->s,opt->ds,opt->q,opt->dq,
                   opt->k,opt->dk,opt->l,opt->dl);
    apd=max(ad,ap);

    do {
      acc=CheckNeighbor(opt,par,apd);
      if (acc) break;
      apd/=0.95;
    } while (apd>=1.0e-13);

    if (apd<1.0e-13)
    {
      opt->m_intSta=NUM_PROB;
      break;
    }

    /*
     * update iterate
     */
    for (k=0; k<nl; k++)
    {
      j         =ibl[k];
      opt->x[j]+=opt->dx[j]/apd;
      opt->g[k]+=opt->dg[k]/apd;
      opt->z[k]+=opt->dz[k]/apd;
    }
    for (k=0; k<nf; k++)
    {
      j         =ibf[k];
      opt->x[j]+=opt->dx[j]/apd;
      opt->f[k]+=opt->df[k]/apd;
      opt->h[k]+=opt->dh[k]/apd;
      opt->k[k]+=opt->dk[k]/apd;
      opt->l[k]+=opt->dl[k]/apd;
    }
    for (k=0; k<nu; k++)
    {
      opt->s[k]+=opt->ds[k]/apd;
      opt->t[k]+=opt->dt[k]/apd;
    }
    for (k=0; k<m; k++)
    {
      opt->y[k]+=opt->dy[k]/apd;
      opt->w[k]+=opt->dw[k]/apd;
      opt->v[k]+=opt->dv[k]/apd;
      opt->p[k]+=opt->dp[k]/apd;
      opt->q[k]+=opt->dq[k]/apd;
    }
  } /* main loop */

//  printf(" -------------------------------------");
//  printf("-------------------------------------\n");

  fprintf(fout,"--------------------------------------");
  fprintf(fout,"-------------------------------------\n");

  if (opt->iter>=par->q_intItr)
    opt->m_intSta=EXC_ITER;

  SaveSol(opt);

  return OPT_OK;
} /* QpSolv */
