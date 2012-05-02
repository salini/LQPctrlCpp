#include "QPdef.h"

static void SolFwdSnode(chol*  sf,
                        int    snde,
                        int    f,
                        int    l,
                        double x[])
{
  int    i,t,sze,*ls,*subg=sf->subg,
         *ujbeg=sf->ujbeg,*uhead=sf->uhead,
         *usub=sf->usub;
  double xi,*l1,*diag=sf->diag,*uval=sf->uval;

  f += subg[snde];
  l += subg[snde];

  for(i=f; i<l; ++i)
  {
    x[i] /= diag[i];
    xi    = x[i];

    ls    = usub+ujbeg[i];
    l1    = uval+uhead[i];
    sze   = l-i-1;

    for(t=0; t<sze; ++t)
      x[ls[t]] -= l1[t]*xi;
  }
} /* SolFwdSnode */

static void SolBward(int    nrow,
                     double diag[],
                     double uval[],
                     int    fir[],
                     double x[])
{
  int    i,t,sze;
  double x1,x2,rtemp,
         *x0,*l1,*l2;

  for(i=nrow; i;) {
    for(; i>1; --i) {
          -- i;
      l1   = uval+fir[i-1]+1;
      l2   = uval+fir[i  ]+0;
      sze  = nrow-i-1;
      x1   = 0.0;
      x2   = 0.0;
      x0   = x+1+i;

      for(t=0; t<sze; ++t)
      {
        rtemp = x0[t];

        x1   += l1[t]*rtemp;
        x2   += l2[t]*rtemp;
      }

      x[i]   -= x2/diag[i];
      x[i-1] -= (uval[fir[i-1]]*x[i]+x1)/diag[i-1];
    }

    for(; i;) {
          -- i;
      l1   = uval+fir[i];
      sze  = nrow-i-1;
      x1   = 0.0;
      x0   = x+1+i;

      for(t=0; t<sze; ++t)
        x1 += l1[t]*x0[t];

      x[i] -= x1/diag[i];
    }
  }
} /* SolBward */

void TriSolv(chol*  sf,
             double b[],
             double x[])
{
  int    i,k,s,t,sze,f,l,itemp,*ls,
         *subg=sf->subg,*ujsze=sf->ujsze,*usub=sf->usub,
         *ujbeg=sf->ujbeg,*uhead=sf->uhead;
  double x1,x2,
         rtemp1,rtemp2,rtemp3,rtemp4,
         rtemp5,rtemp6,rtemp7,rtemp8,
         *l1,*l3,*l2,*l4,*l5,*l6,*l7,*l8,
         *diag=sf->diag,*uval=sf->uval;
   
  for(i=0; i<sf->nrow; ++i)
    x[i] = b[sf->perm[i]];

  for(s=0; s<sf->nsnds; ++s) {
    f = subg[s];
    l = subg[s+1];

    SolFwdSnode(sf,s,0,l-f,x);

    itemp = l-f-1;
    ls    = usub+ujbeg[f]+itemp;
    sze   = ujsze[f]-itemp;
    k     = f;

    itemp = l-1;
    for(; k+7<l; k+=8) {
      l1       = uval+uhead[k+0]+itemp-(k+0);
      l2       = uval+uhead[k+1]+itemp-(k+1);
      l3       = uval+uhead[k+2]+itemp-(k+2);
      l4       = uval+uhead[k+3]+itemp-(k+3);
      l5       = uval+uhead[k+4]+itemp-(k+4);
      l6       = uval+uhead[k+5]+itemp-(k+5);
      l7       = uval+uhead[k+6]+itemp-(k+6);
      l8       = uval+uhead[k+7]+itemp-(k+7);

      rtemp1   = x[k+0];
      rtemp2   = x[k+1];
      rtemp3   = x[k+2];
      rtemp4   = x[k+3];
      rtemp5   = x[k+4];
      rtemp6   = x[k+5];
      rtemp7   = x[k+6];
      rtemp8   = x[k+7];

      for(t=0; t<sze; ++t)
        x[ls[t]] -=   rtemp1*l1[t]
                    + rtemp2*l2[t]
                    + rtemp3*l3[t]
                    + rtemp4*l4[t]
                    + rtemp5*l5[t]
                    + rtemp6*l6[t]
                    + rtemp7*l7[t]
                    + rtemp8*l8[t];
    }

    for(; k+3<l; k+=4) {
      l1       = uval+uhead[k+0]+itemp-(k+0);
      l2       = uval+uhead[k+1]+itemp-(k+1);
      l3       = uval+uhead[k+2]+itemp-(k+2);
      l4       = uval+uhead[k+3]+itemp-(k+3);

      rtemp1   = x[k+0];
      rtemp2   = x[k+1];
      rtemp3   = x[k+2];
      rtemp4   = x[k+3];

      for(t=0; t<sze; ++t)
        x[ls[t]] -=   rtemp1*l1[t]
                    + rtemp2*l2[t]
                    + rtemp3*l3[t]
                    + rtemp4*l4[t];
    }

    for(; k+1<l; k+=2) {
      l1       = uval+uhead[k+0]+itemp-(k+0);
      l2       = uval+uhead[k+1]+itemp-(k+1);

      rtemp1   = x[k+0];
      rtemp2   = x[k+1];

      for(t=0; t<sze; ++t)
        x[ls[t]] -=   rtemp1*l1[t]
                    + rtemp2*l2[t];
    }


    for(; k<l; ++k) {
      l1       = uval+uhead[k+0]+itemp-(k+0);

      rtemp1   = x[k+0];

      for(t=0; t<sze; ++t)
        x[ls[t]] -=   rtemp1*l1[t];
    }
  }

  if (sf->nsnds) {
    s = sf->nsnds - 1;
    f = subg[s];
    l = subg[s+1];

    dCopy(l-f,x+f,b+f);

    SolBward(l-f,diag+f,uval,uhead+f,b+f);

    s = sf->nsnds-1;

    for(; s>=1; --s) {
      f = subg[s-1];
      l = subg[s];
      i = l;

      for(; i>1+f; --i) {
            -- i;
        ls   = usub+ujbeg[i];
        l1   = uval+uhead[i-1]+1;
        l2   = uval+uhead[i  ]+0;
        sze  = ujsze[i];
        x1   = 0.0;
        x2   = 0.0;

        for(t=0; t<sze; ++t) {
          rtemp1 = b[ls[t]];

          x1    += l1[t]*rtemp1;
          x2    += l2[t]*rtemp1;
        }

        b[i]   = x[i  ] -  x2  / diag[i];
        b[i-1] = x[i-1] - (x1 + uval[uhead[i-1]]*b[i]) / diag[i-1];
      }

      for(; i>f;) {
            -- i;
        l1   = uval+uhead[i];
        ls   = usub+ujbeg[i];
        sze  = ujsze[i];
        x1   = 0.0;

        for(t=0; t<sze; ++t)
          x1+= l1[t]*b[ls[t]];

        b[i] = x[i] - x1/diag[i];
      }
    }
  }

  for(i=0; i<sf->nrow; ++i)
    x[i] = b[sf->invp[i]];

} /* TriSolv */

