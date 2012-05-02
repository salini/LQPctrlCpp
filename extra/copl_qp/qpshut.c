#include "QPdef.h"

void ShutDown(int code)
{
  hFree(&qt->rowtab);
  hFree(&qt->coltab);   
  rFree(&qt->col);
  cFree(&qt->mpscard);
  cFree(&qt->prbname);
  cFree(&qt->objname);
  cFree(&qt->rowtype);
  cPfre(&qt->rowname,qt->m_intRow);
  cPfre(&qt->colname,qt->m_colUse);
  iFree(&qt->imem);
  iFree(&qt->bid);
  dFree(&qt->rmem);

  dFree(&qt->b);
  dFree(&qt->c);
  dFree(&qt->l);
  dFree(&qt->u);
  dFree(&qt->r);

  dFree(&qt->rp);
  dFree(&qt->rd);
  dFree(&qt->rl);
  dFree(&qt->ru);
  dFree(&qt->rr);
  dFree(&qt->ry);
  dFree(&qt->rf);
  dFree(&qt->rh);
  dFree(&qt->rz);
  dFree(&qt->rw);
  dFree(&qt->rs);
  dFree(&qt->rq);
  dFree(&qt->rk);
  dFree(&qt->ro);
  dFree(&qt->qx);

  dFree(&qt->diag);
  dFree(&qt->diau);
  dFree(&qt->diaf);
  dFree(&qt->dial);

  dFree(&qt->x);
  dFree(&qt->y);
  dFree(&qt->g);
  dFree(&qt->z);
  dFree(&qt->t);
  dFree(&qt->s);
  dFree(&qt->v);
  dFree(&qt->w);
  dFree(&qt->p);
  dFree(&qt->q);
  dFree(&qt->f);
  dFree(&qt->h);
  dFree(&qt->k);   
    
  dFree(&qt->dx);
  dFree(&qt->dy);
  dFree(&qt->dg);
  dFree(&qt->dz);
  dFree(&qt->dt);
  dFree(&qt->ds);
  dFree(&qt->dv);
  dFree(&qt->dw);
  dFree(&qt->dp);
  dFree(&qt->dq);
  dFree(&qt->df);
  dFree(&qt->dh);
  dFree(&qt->dk);
  dFree(&qt->dl);

  mFree(&qt->a);
  mFree(&qt->at);
  sFree(&qt->qh);
  chFre(&qt->cl);

  printf("\n");
  fprintf(fout,"\n");

  switch (code)
  {
    case OPT_FOUND:
      printf(" Exit -- %d: optimal solution found.\n",code);
      fprintf(fout,"Exit -- %d: optimal solution found.\n",code);
      break;
    case UNB_PROB:
      printf(" Exit -- %d: infeasible problem.\n",code);
      fprintf(fout,"Exit -- %d: infeasible problem.\n",code);
      break;
    case INF_PROB:
      printf(" Exit -- %d: infeasible problem.\n",code);
      fprintf(fout,"Exit -- %d: infeasible problem.\n",code);
      break;
    case NOT_MEMSPC:
      printf(" Exit -- %d: out of memory.\n",code);
      fprintf(fout,"Exit -- %d: out of memory.\n",code);
      break;
    case NOT_DSKSPC:
      printf(" Exit -- %d: out of disk space.\n",code);
      fprintf(fout,"Exit -- %d: out of disk space.\n",code);
      break;
    case SYS_ERROR:
      printf(" Exit -- %d: system error.\n",code);
      fprintf(fout,"Exit -- %d: system error.\n",code);
      break;
    case NOT_SCHFIL:
      printf(" Exit -- %d: can't open the specified file.\n",code);
      fprintf(fout,"Exit -- %d: can't open the specified file.\n",code);
      break;
    case NOT_PRIM:
      printf(" Exit -- %d: larger prime number needed.\n",code);
      fprintf(fout,"Exit -- %d: larger prime number needed.\n",code);
      break;
    case FAC_ERROR:
      printf(" Exit -- %d: Numerical Difficulty.\n",code);
      fprintf(fout,"Exit -- %d: Numerical Difficulty.\n",code);
      break;
    case NUM_PROB:
      printf(" Exit -- %d: Numerical Difficulty.\n",code);
      fprintf(fout,"Exit -- %d: Numerical Difficulty.\n",code);
      break;
    default:
      printf( " Exit -- %d: system error.\n",code);
      fprintf(fout,"Exit -- %d: system error.\n",code);
      break;
  }
    
  t_end=GetTime();

#ifdef UNIXMACHINE
  printf("\n the elapsed time is : %.2f seconds.\n",
        (double)(t_end-t_start)/100.0);
  fprintf(fout,"\nthe elapsed time is : %.2f seconds.\n",
        (double)(t_end-t_start)/100.0);
#else
  printf("\n the elapsed time is : %.2f\n",
        (double)(t_end-t_start)/1000.0);
  fprintf(fout,"\nthe elapsed time is : %.2f\n",
        (double)(t_end-t_start)/1000.0);
#endif

  printf("\n ******************************* COPL ENDS");
  printf(" ********************************\n\n");
  fprintf(fout,"\n******************************* COPL ENDS");
  fprintf(fout," ********************************\n\n");
} /* ShutDown */
