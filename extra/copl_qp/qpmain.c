#include "QPdef.h"

/* global variables */
FILE *fmps=NULL,
     *fout=NULL,
     *fres=NULL;
char *nmps=NULL,
     *nout=NULL,
     *nspc=NULL;

qpdat *qt =NULL;
qppar *pa =NULL;

clock_t t_start;
clock_t t_end;

clock_t GetTime(void)
{
  clock_t t;

#ifdef UNIXMACHINE  
  struct tms {
    clock_t    tms_utime;
    clock_t    tms_stime;
    clock_t    tms_cutime;
    clock_t    tms_cstime;
  } tmrec;

  clock_t times();

  t = times(&tmrec);
  t = tmrec.tms_utime+tmrec.tms_stime;
#else
  t = clock();
#endif

  return t;
} /* GetTime */

static void PrintHead(void)
{
  char Title[8][60]={
    " ==========================================\n",
    " Computational Optimization Program Library\n",
    "                                           \n",
    "     (  C  O  P  L  )     Version 1.0      \n",
    "                                           \n",
    " Computational Optimization Lab.           \n",
    " The University of Iowa (1997)             \n",
    " ==========================================\n"};
  int i;

  printf("\n\n\n");
  for (i=0; i<8; i++)
    printf("%s", Title[i]);
  printf("\n ****************************** COPL STARTS");
  printf(" *******************************\n\n");

  return;
} /* PrintHead */

static int FprintHead(void)
{
  fout=fopen(nout,"w");
  if (!fout)
    return NOT_DSKSPC;

  fprintf(fout,"====================\n");
  fprintf(fout,"C O P L  1.0  (1997)\n");
  fprintf(fout,"====================\n\n");
  return true;
} /* FprintHead */

static int CheckArgs(int   argc,
                     char* argv[],
                     char* ss)
{
  int i,id;
  
  id = 0;
  for (i=1; i<argc; i++) {
    if(strstr(argv[i],ss)) {
      argv[i]+=2;
      if (argv[i][0] == ':') {
        argv[i]++;
        id=i;
        break;
      }
      if (argv[i][0]) {
        id=i;
        break;
      }
      id = i+1;
      if (argv[id][0] == '-')
        id = 0;
      break;
    }
  }
  return id;
} /* Checkargs */

static void setParams(qppar* par)
{
  int    i,k;
  double rtmp;
  FILE*  fp=NULL;
  char   str[80],*ptr;
  char   parname[][50]=
  {
    "bound_set_name",
    "cache_memory_size",
    "dual_feasible_tolerance",
    "initial_point_threshold",
    "iteration_limit",
    "lower_bound_default",
    "number_of_columns",
    "objective_set_name",
    "primal_dual_gap_tolerance",
    "primal_feasible_tolerance",
    "problem_type",
    "range_set_name",
    "rignt_hand_side_set_name",
    "step_factor",
    "upper_bound_default",
    ""
  };

  par->q_intCol=       100000;
  par->q_intItr=          200;
  par->q_intCac=          256;
  par->q_intMin=         true;
    
  par->q_strRhs=         NULL;
  par->q_strRng=         NULL;
  par->q_strBnd=         NULL;
  par->q_strObj=         NULL;
    
  par->q_dblLow=          0.0;
  par->q_dblUpp=        1e+40;
  par->q_dblMax=0.999*_P_MAXM;
  par->q_dblAij=        1e-10;    
  par->q_dblPiv=        1e-18;
  par->q_dblGap=        1e-06;
  par->q_dblPri=        1e-06;
  par->q_dblDua=        1e-06;
  par->q_dblStp=         0.95;
  par->q_dblRad=        1e-08;
  par->q_dblIni=        1000.;
  
  if (nspc)
  {
    fp=fopen(nspc,"r");
    if (!fp)
    {
      ShutDown(NOT_SCHFIL);
      exit(NOT_SCHFIL);
    }
  }
  else
  {
    fp=fopen("COPLQP.spc","r");
    if (!fp) return;
  }
  
  fprintf(fout,"SPC FILE\n");
  fprintf(fout,"--------\n");
  while (!feof(fp))
  {
    fgets(str,80,fp);
    fputs(str,fout);
    if (*str=='*'||*str=='\n'||*str=='\0')
      continue;
    
    for (ptr=str; *ptr==' '; ptr++);

    for (i=0; i<15; i++)
    {
      k=strlen(parname[i]);
      if (strncmp(ptr,parname[i],k)==0)
        break;
    }
    
    if (i>=15)
    {
      fprintf(fout,"XXXXX invalid parameter item: %s\n",str);
      continue;
    }
    
    for (ptr=ptr+k;
         *ptr==' '||*ptr==':'||*ptr=='=';
         ptr++);
    
    switch(i)
    {
    case 0: /* bound name */
      k=strlen(ptr);
      ptr[k-1]='\0';
      if (strcmp(ptr,"NULL"))
      {
        par->q_strBnd=cAlloc(k);
        strncpy(par->q_strBnd,ptr,k);
      }
      break;
    case 1: /* cache memory */
      k=atoi(ptr);
      if (k>0) par->q_intCac=k;
      break;
    case 2: /* told */
      rtmp=atof(ptr);
      if (rtmp>0.0) par->q_dblDua=rtmp;
      break;
    case 3: /* init point threshold */
      rtmp=atof(ptr);
      if (rtmp>0.0) par->q_dblIni=rtmp;
      break;
    case 4: /* max iter */
      k=atoi(ptr);
      if (k>0) par->q_intItr=k;
      break;
    case 5: /* lower bound */
      rtmp=atof(ptr);
      par->q_dblLow=rtmp;
      break;
    case 6: /* ncol */
      k=atoi(ptr);
      if (k>0) par->q_intCol=k;
      break;
    case 7: /* obj name */
      k=strlen(ptr);
      ptr[k-1]='\0';
      if (strcmp(ptr,"NULL"))
      {
        par->q_strObj=cAlloc(k);
        strncpy(par->q_strObj,ptr,k);
      }
      break;
    case 8: /* tolg */
      rtmp=atof(ptr);
      if (rtmp>0.0) par->q_dblGap=rtmp;
      break;
    case 9: /* tolp */
      rtmp=atof(ptr);
      if (rtmp>0.0) par->q_dblPri=rtmp;
      break;
    case 10: /* probtype */
      k=strlen(ptr);
      ptr[k-1]='\0';
      if (strcmp(ptr,"min")==0||
          strcmp(ptr,"Min")==0||
          strcmp(ptr,"MIN")==0)
        par->q_intMin=true;
      if (strcmp(ptr,"max")==0||
          strcmp(ptr,"Max")==0||
          strcmp(ptr,"MAX")==0)
        par->q_intMin=false;      
      break;
    case 11: /* range name */
      k=strlen(ptr);
      ptr[k-1]='\0';
      if (strcmp(ptr,"NULL"))
      {
        par->q_strRng=cAlloc(k);
        strncpy(par->q_strRng,ptr,k);
      }
      break;
    case 12: /* rhs name */
      k=strlen(ptr);
      ptr[k-1]='\0';
      if (strcmp(ptr,"NULL"))
      {
        par->q_strRhs=cAlloc(k);
        strncpy(par->q_strRhs,ptr,k);
      }
      break;
    case 13: /* step fac */
      rtmp=atof(ptr);
      if (rtmp>0.0) par->q_dblStp=rtmp;
      break;
    case 14: /* upper bound */
      rtmp=atof(ptr);
      par->q_dblUpp=rtmp;
      break;
    default:
      break;
    }
  }
  fprintf(fout,"\n");
  fclose(fp);
} /* setParams */

static void GetFile(char fname[])
{
  int  k;
  char *ptr;

  printf(" Please type the name of your data file: ");
  fflush(stdin);
  fgets(fname,128,stdin);
  k=strlen(fname);
  ptr=fname+k-1;
  while(k>0&&(*ptr=='\n'||*ptr==' '))
  {
    *ptr='\0';
    ptr--;
    k--;
  }
} /* GetFile */

int main(int argc,char* argv[])
{  
  /* local variables */
  int   idKey; 
  qpdat qp={0};

  t_start=GetTime();

  qt=&qp;
  pa=&qp.par;
  
  /* initialization */
  PrintHead();

  idKey=CheckArgs(argc,argv,"-o");
  if (!idKey) {
    nout=cAlloc(20);
    sprintf(nout,"coplfile.out");
  }
  else nout=argv[idKey];

  if (!FprintHead())
    return NOT_DSKSPC;
  
  idKey=CheckArgs(argc,argv,"-s");
  if (idKey) nspc=argv[idKey];

  setParams(&qp.par);
  
  /* input data */
  idKey=CheckArgs(argc,argv,"-f");
  if (!idKey) {
    nmps =cAlloc(128);
    GetFile(nmps);
  }
  else nmps=argv[idKey];

  idKey=DataInput(&qp.par,&qp,nmps);
  if (OPT_OK!=idKey)
    return idKey;
  
  /* preprocess */
  idKey=PreProc(&qp.par,&qp);
  if (OPT_OK!=idKey)
    return idKey;
  
  /* symbolic computation */
  idKey=SymboProc(&qp,&qp.par);
  if (OPT_OK!=idKey)
    return idKey;
  
  /* numerical computation */
  idKey=QpSolv(&qp,&qp.par);
  
  ShutDown(qp.m_intSta);
  fclose(fout);

  return OPT_OK;
} /* main */
