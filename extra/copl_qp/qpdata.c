#include "QPdef.h"

static unsigned secfield[8]={4,4,7,3,6,6,5,6};

static char     sec[8][10]= {
                  "NAME","ROWS","COLUMNS",
                  "RHS","RANGES","BOUNDS",
                  "QUADS","ENDATA"
                },
                bnd[6][5]= {
                  " LO "," UP "," FX ",
                  " FR "," MI "," PL "
                };

static int GetCard(char* card,
                   int*  ncrd)
{
  do {
    if (feof(fmps)) {
      ShutDown(NOT_ENDCARD);
      return (NOT_ENDCARD);
    }
    fgets(card,_L_SIZE,fmps);
    (*ncrd)++;
  } while(*card=='*'||*card=='\n');

  return (((*card)==' ')? MIDCARD:HEADCARD);
} /* GetCard */

static int ParseCard(char* card,
                     char* key[])
{
  key[0]=card+4;
  key[1]=card+14;
  key[2]=card+24;
  key[3]=card+39;
  key[4]=card+49;
  return true;
} /* ParseCard */

static int GetName(char* s1,
                   char* s2)
{
  int i,j;
  
  j = strlen(s2)-1;
  if (j>8) j = 8;
  while(j>=0) {
    j--;
    if (s2[j] !=' '&&s2[j] !='\n') break;
  }
  for (i=0; i<=j; i++)
    s1[i]=s2[i];
  s1[i]='\0';
  return j+1;
} /* GetName */

static int HashValue(hashtab  *tab,
                     char     *key,
                     unsigned scal)
{
  unsigned val;
  char     *s;
  
  val=0;
  for (s=key; *s!='\0'; s++)
    val=scal*val+(*s);
    
  return (val%tab->lssz);
} /* HashValue */

static int CompString(const char *s1,
                      const char *s2)
{
  if (!s1||!s2) return true;
  return strcmp(s1,s2); 
} /* CompString */

static int InsertHash(hashtab  *tab,
                      char     *key,
                      int      addr,
                      unsigned scal)
{
  hashlst *ptr;
  int     k,len;
  
  k=HashValue(tab,key,scal);
  
  for (ptr=tab->list[k]; ptr; ptr=ptr->next)
    if (!CompString(ptr->key,key)) break;
    
  if (!ptr) {
    ptr=(hashlst*)cAlloc(sizeof(hashlst));
    
    len=strlen(key);
    ptr->key=cAlloc(len+1);
    
    strcpy(ptr->key,key);
    ptr->addr=addr;
    
    ptr->next=tab->list[k];
    tab->list[k]=ptr;
    return false;
  }
  
  return true;
} /* InsertHash */

static int FoundHash(hashtab  *tab,
                     char     *key,
                     unsigned scal)
{
  int     k,addr;
  hashlst *ptr;
  
  addr=-1;
  k=HashValue(tab,key,scal);
  
  for (ptr=tab->list[k]; ptr; ptr=ptr->next) {
    if (!CompString(ptr->key,key)) {
      addr=ptr->addr;
      break;
    }
  }
  
  return addr;
} /* FoundHash */

static int SecProc(mpssec idSec,
                   qppar* par,
                   qpdat* opt,
                   int*   ncrd,
                   int*   nerr)
{
  int             h,i,k,ifd,ityp,m,n,nnzo,*colsub,colz;
  char            *ptr,*key[5],type,name[9];
  double          cost,rtmp,*colval;
  static unsigned rowscal,colscal;
  int             toupper();
  
  switch(idSec)
  {
    case NAME:
      rewind(fmps);
      
      ptr=opt->mpscard;       
      if (GetCard(ptr,ncrd)==HEADCARD)
      {
        if (strncmp(ptr,sec[idSec],secfield[idSec]))
        {
          ShutDown(NOT_NAMECARD);
          return (NOT_NAMECARD);
        }
        key[1]=ptr+14;
        k     =strlen(key[1]);
        opt->prbname=cAlloc(k+1);
        GetName(opt->prbname,key[1]);
        printf(" NAME: %s\n",opt->prbname);
        fprintf(fout,"%8d    ",*ncrd);
        fputs(ptr,fout);
      }

      break;

    case ROWS:
      if (GetCard(opt->mpscard,ncrd)==HEADCARD)
      {
        ptr=opt->mpscard;
        if (strncmp(ptr,sec[idSec],secfield[idSec]))
        {
          ShutDown(NOT_ROWSCARD);
          return NOT_ROWSCARD;
        }
        printf(" ROWS section ...\n");
        fprintf(fout,"%8d    %s",*ncrd,ptr);
      }
      else
      {
        fprintf(fout,"%8d    %s",*ncrd,opt->mpscard);
        fprintf(fout," ROWS card is missing here.\n");
        ShutDown(NOT_ROWSCARD);
        return (NOT_ROWSCARD);
      }
      
      for (k=0; GetCard(opt->mpscard,&i)==MIDCARD; k++);
      
      opt->rowtype=cAlloc(k);
      opt->rowtab =hAlloc(k);
      opt->rowname=cPallo(k,0);
      rowscal=(unsigned)(_G_FACT*opt->rowtab->lssz);
      
      rewind(fmps);
      do
      {
        GetCard(ptr,&i);
      } while (strncmp(ptr,sec[idSec],secfield[idSec]));
      key[0]=ptr+4;
      
      while(GetCard(opt->mpscard,ncrd)==MIDCARD)
      {
        ptr=opt->mpscard+1;
        
        if (*ptr==' ') ptr++;
        type =toupper(*ptr);
        
        ityp=(type=='N'||type=='E'||
              type=='G'||type=='L');
        
        if (!ityp)
        {
          (*nerr)++;
          fprintf(fout,"%8d    warning -- %d: ",
                  *ncrd,NOT_ROW_TYPE);
          fprintf(fout,"not exist row type %s\n",opt->mpscard);
          continue;
        }
        
        if (type=='N')
        {
          if (opt->objname) continue;
          k=strlen(key[0]);
          opt->objname=cAlloc(k+1);
          k=GetName(opt->objname,key[0]);
          continue;
        }
        
        m=opt->m_intRow;
        k=GetName(name,key[0]);
        
        if (InsertHash(opt->rowtab,name,m,rowscal))
        {
          (*nerr)++;
          fprintf(fout,"%8d    warning -- %d: ",
                  *ncrd,DPL_ROW_NAME);
          fprintf(fout,"duplicate row name %s.\n",name);
          continue;
        }
        
        if (type=='G'||type=='L')
          opt->m_intUne++;
        
        opt->rowtype[m]=type;
        opt->rowname[m]=cAlloc(k+1);
        sprintf(opt->rowname[m],name);
        opt->m_intRow++;
      }
      
      break;
      
    case COLUMNS:
      ptr=opt->mpscard;
      ParseCard(ptr,key);
      if (strncmp(ptr,sec[idSec],secfield[idSec]))
      {
        ShutDown(NOT_COLSCARD);
        return (NOT_ROWSCARD);
      }
        
      printf(" COLUMNS section ...\n");
      fprintf(fout,"%8d    %s",*ncrd,ptr);

      opt->col    =rAlloc(par->q_intCol);
      opt->coltab =hAlloc(par->q_intCol);
      opt->colname=cPallo(par->q_intCol,0);
      colscal     =(unsigned)(_G_FACT*opt->coltab->lssz);
      colsub      =iAlloc(opt->m_intRow);
      colval      =dAlloc(opt->m_intRow);
      n           =0;
      colz        =0;
      nnzo        =0;
      cost        =0.0;

      while(GetCard(ptr,ncrd)==MIDCARD)
      {        
        h=GetName(name,key[0]);
        if (CompString(opt->colname[n],name))
        {
          opt->m_intCol++;
          n=opt->m_intCol-1;
          if (n>par->q_intCol) 
            continue;

          if (InsertHash(opt->coltab,name,n,colscal))
          {
            (*nerr)++;
            fprintf(fout,"%8d    warning -- %d: ",
                    *ncrd,DPL_COL_NAME);
            fprintf(fout,"duplicate column name %s.\n",name);
            opt->m_intCol--;
            continue;
          }
            
          opt->colname[n]=cAlloc(h+1);
          sprintf(opt->colname[n],"%s",name);
          opt->m_colUse++;

          if (opt->m_intCol>1)
          {
            h              =n-1;
            opt->col[h].nn0=colz;
            opt->col[h].rhs=cost;
            opt->col[h].low=par->q_dblLow;
            opt->col[h].upp=par->q_dblUpp;
            opt->col[h].ja =iAlloc(colz);
            opt->col[h].an =dAlloc(colz);
            for (k=0; k<colz; k++)
            {
              opt->col[h].ja[k]=colsub[k];
              opt->col[h].an[k]=colval[k];
            }
            opt->m_intAnz+=colz;
          }
          
          colz=0;
          cost=0.0;
        }
        
        k=strlen(ptr);
        for (ifd=VALUE1,i=1; ifd<k; i+=2,ifd+=VALUE2)
        {
          h=GetName(name,key[i]);
          rtmp=atof(key[i+1]);
          if (!CompString(opt->objname,name))
          {
            if (par->q_intMin) cost=rtmp;
            else cost=-rtmp;
            continue;
          }
          
          m=FoundHash(opt->rowtab,name,rowscal);
          if (m<0) 
          {
            (*nerr)++;
            fprintf(fout,"%8d    warning -- %d: ",
                    *ncrd,NOT_ROW_NAME);
            fprintf(fout,"nonexisted row name %s.\n",name);
            continue;
          }
          
          if (fabs(rtmp)<par->q_dblAij)
          {
            (*nerr)++;
            fprintf(fout,"%8d    warning -- %d: ",
                    *ncrd,ZERO_ELEM);
            fprintf(fout,"zero element found: %s.\n",key[i]);
            continue;
          }
          
          colval[colz]=rtmp;
          colsub[colz]=m;
          colz++;
        }
      }
      
      if (opt->m_intCol>par->q_intCol)
      {
        iFree(&colsub);
        dFree(&colval);
        ShutDown(EXC_COLLIM);
        return (EXC_COLLIM); 
      }
      opt->col[n].nn0=colz;
      opt->col[n].rhs=cost;
      opt->col[n].low=par->q_dblLow;
      opt->col[n].upp=par->q_dblUpp;
      opt->col[n].ja =iAlloc(colz);
      opt->col[n].an =dAlloc(colz);
      for (k=0; k<colz; k++)
      {
         opt->col[n].ja[k]=colsub[k];
         opt->col[n].an[k]=colval[k];
      }
      opt->m_intAnz+=colz;
      
      iFree(&colsub);
      dFree(&colval);
      
      opt->b=dAlloc(opt->m_intRow);
      opt->r=dAlloc(opt->m_intRow);
      memset(opt->b,0,opt->m_intRow*sizeof(double));
      memset(opt->r,0,opt->m_intRow*sizeof(double));

      break;
      
    case RHS:
      ptr=opt->mpscard;
      ParseCard(ptr,key);

      if (strncmp(ptr,sec[idSec],secfield[idSec]))
        printf(" Null RHS section\n");

      else 
      {
        printf(" RHS section ...\n");
        fprintf(fout,"%8d    %s",*ncrd,ptr);
        
        while(GetCard(ptr,ncrd)==MIDCARD)
        {
          h=GetName(name,key[0]);
          if (!par->q_strRhs)
          {
            par->q_strRhs=cAlloc(h+1);
            sprintf(par->q_strRhs,"%s",name);
          }
          else if (CompString(par->q_strRhs,name))
            continue;
          
          k=strlen(ptr);
          for (ifd=VALUE1,i=1; ifd<k; i+=2,ifd+=VALUE2)
          {
            h=GetName(name,key[i]);
            rtmp=atof(key[i+1]);
            
            m=FoundHash(opt->rowtab,name,rowscal);
            if (m<0) 
            {
              (*nerr)++;
              fprintf(fout,"%8d    warning -- %d: ",
                      *ncrd,NOT_ROW_NAME);
              fprintf(fout,"nonexisted row name %s.\n",name);
              continue;
            }
            
            if (fabs(rtmp)<par->q_dblAij)
            {
              (*nerr)++;
              fprintf(fout,"%8d    warning -- %d: ",
                      *ncrd,ZERO_ELEM);
              fprintf(fout,"zero element found: %s.\n",key[i]);
              continue;
            }
            
            opt->b[m]=rtmp;
          }
        }
      }
      break;
      
    case RANGES:
      for (i=0; i<opt->m_intRow; i++)
        opt->r[i]=opt->b[i];
        
      ptr=opt->mpscard;
      ParseCard(ptr,key);

      if (!strncmp(ptr,sec[idSec],secfield[idSec]))
      {
        printf(" RANGES section ...\n");
        fprintf(fout,"%8d    %s",*ncrd,ptr);
        
        while(GetCard(ptr,ncrd)==MIDCARD)
        {
          h=GetName(name,key[0]);
          if (!par->q_strRng)
          {
            par->q_strRng=cAlloc(h+1);
            sprintf(par->q_strRng,"%s",name);
          }
          else if (CompString(par->q_strRng,name))
            continue;
          
          k=strlen(ptr);
          for (ifd=VALUE1,i=1; ifd<k; i+=2,ifd+=VALUE2)
          {
            h=GetName(name,key[i]);
            rtmp=atof(key[i+1]);
            
            m=FoundHash(opt->rowtab,name,rowscal);
            if (m<0) 
            {
              (*nerr)++;
              fprintf(fout,"%8d    warning -- %d: ",
                      *ncrd,NOT_ROW_NAME);
              fprintf(fout,"nonexisted row name %s.\n",name);
              continue;
            }
            
            cost=fabs(rtmp);
            switch(opt->rowtype[m])
            {
            case 'E':
              if (rtmp<0)
                opt->r[m]-=cost;
              else
                opt->b[m]+=cost;
              break;   
            case 'L':
              opt->r[m]-=cost;
              opt->m_intUne--;
              break;
            case 'G':
              opt->b[m]+=cost;
              opt->m_intUne--;
              break;  
            case 'R':
              (*nerr)++;
              fprintf(fout,"%8d    warning -- %d: ",
                      *ncrd,DPL_RNG_ROW);
              fprintf(fout,"duplicate range setting %s.\n",name);
              break;  
            default:
              break;
            } /* switch */
            
            opt->rowtype[m]='R';
            opt->m_intNrg++;
            
          } /* for */
        } /* while */
      } /* if */
      break;
      
    case BOUNDS:
      ptr=opt->mpscard;
      ParseCard(ptr,key);
      
      if (!strncmp(ptr,sec[idSec],secfield[idSec]))
      {
        printf(" BOUNDS section ...\n");
        fprintf(fout,"%8d    %s",*ncrd,ptr);
        
        while(GetCard(ptr,ncrd)==MIDCARD)
        {
          h=GetName(name,key[0]);
          if (!par->q_strBnd)
          {
            par->q_strBnd=cAlloc(h+1);
            sprintf(par->q_strBnd,"%s",name);
          }
          else if (CompString(par->q_strBnd,name))
            continue;
          
          k=0;
          for (i=LOW; i<=PLS; i++) 
          {
            if (!strncmp(ptr,bnd[i],4)) 
            {
              k=1;
              break;
            }
          }
          
          if (!k) 
          {
            (*nerr)++;
            fprintf(fout,"%8d    warning -- %d: ",
                    *ncrd,NOT_BND_TYPE);
            fprintf(fout,"nonexisted bound type %s.\n",ptr);
            continue;
          }
          
          GetName(name,key[1]);
          n=FoundHash(opt->coltab,name,colscal);
          if (n<0) 
          {
            (*nerr)++;
            fprintf(fout,"%8d    warning -- %d: ",
                    *ncrd,NOT_COL_NAME);
            fprintf(fout,"nonexisted row name %s.\n",name);
            continue;
          }
          if ((int)strlen(ptr)<24) rtmp=0.0;
          else rtmp=atof(key[2]);
          
          switch (i)
          {
          case LOW:
            opt->col[n].low=rtmp;
            opt->m_intNlo++;
            break;  
          case UPP:
            opt->col[n].upp=rtmp;
            opt->m_intNup++;
            break;
          case FIX:
            opt->col[n].low=rtmp;
            opt->col[n].upp=rtmp;
            opt->m_intNfx++;
            break;  
          case FRE:
            opt->col[n].upp= _P_MAXM;
            opt->col[n].low=-_P_MAXM;
            opt->m_intNfr++;
            break;
          case MIS:
            opt->col[n].low=-_P_MAXM;
            opt->col[n].upp=rtmp;
            opt->m_intNmi++;
            break;  
          case PLS:
            opt->col[n].upp=_P_MAXM;
            opt->col[n].low=rtmp;
            opt->m_intNpl++;
            break;
          default:
            break;
          }     
        } /* while */
      } /* if */
      break;
      
    case QUADS:
      colsub =iAlloc(opt->m_colUse);
      colval =dAlloc(opt->m_colUse);
      opt->qh=sAlloc(opt->m_colUse,0);
      ptr    =opt->mpscard;

      if (!strncmp(ptr,sec[idSec],secfield[idSec]))
      {
        printf(" QUADS section ...\n");
        fprintf(fout,"%8d    %s",*ncrd,ptr);
        
        h   =-1;
        colz=0;
        while(GetCard(ptr,ncrd)==MIDCARD)
        {          
          k=GetName(name,key[0]);
          n=FoundHash(opt->coltab,name,colscal);
          if (n<0) 
          {
            (*nerr)++;
            fprintf(fout,"%8d    warning -- %d: ",
                    *ncrd,NOT_COL_NAME);
            fprintf(fout,"nonexisted row name %s.\n",name);
            continue;
          }
          if (n!=h)
          {
            if (colz>0)
            {
              opt->qh->cols[h].nn0=colz;
              opt->qh->cols[h].ja =iAlloc(colz);
              opt->qh->cols[h].an =dAlloc(colz);
              for (i=0; i<colz; i++)
              {
                opt->qh->cols[h].ja[i]=colsub[i];
                opt->qh->cols[h].an[i]=colval[i];
              }
              opt->m_intQnz+=colz;
            }
            colz=0;
            h   =n;
          }
          
          k=strlen(ptr);
          for (ifd=VALUE1,i=1; ifd<k; i+=2,ifd+=VALUE2)
          {
            ityp=GetName(name,key[i]);
            rtmp=atof(key[i+1]);
            
            m=FoundHash(opt->coltab,name,colscal);
            if (m<0) 
            {
              (*nerr)++;
              fprintf(fout,"%8d    warning -- %d: ",
                      *ncrd,NOT_ROW_NAME);
              fprintf(fout,"nonexisted row name %s.\n",name);
              continue;
            }
          
            if (fabs(rtmp)<par->q_dblAij)
            {
              (*nerr)++;
              fprintf(fout,"%8d    warning -- %d: ",
                      *ncrd,ZERO_ELEM);
              fprintf(fout,"zero element found: %s.\n",key[i]);
              continue;
            }
            
            if (m==h)
            {
              if (par->q_intMin) opt->qh->diag[h]=rtmp;
              else opt->qh->diag[h]=-rtmp;
              opt->ndia++;
            }
            else
            {
              colsub[colz]=m;
              if (par->q_intMin) colval[colz]=rtmp;
              else colval[colz]=-rtmp;
              colz++;
            }
          } /* for */
        } /* while */
        
        if (colz>0)
        {
          opt->qh->cols[h].nn0=colz;
          opt->qh->cols[h].ja =iAlloc(colz);
          opt->qh->cols[h].an =dAlloc(colz);
          for (i=0; i<colz; i++)
          {
            if (colsub[i]==h)
              opt->qh->diag[h]=colval[i];
            else
            {
              opt->qh->cols[h].ja[i]=colsub[i];
              opt->qh->cols[h].an[i]=colval[i];
            }
          }
          opt->m_intQnz+=colz;
        }
        iFree(&colsub);
        dFree(&colval);
      }
      
      break;
      
    case ENDATA:
      ptr=opt->mpscard;
      if (GetCard(ptr,ncrd)==HEADCARD)
      {
        if (strncmp(ptr,sec[idSec],secfield[idSec]))
        {
          ShutDown(NOT_ENDCARD);
          return (NOT_ENDCARD);
        }
        printf(" ENDATA\n");
        fprintf(fout,"%8d    ",*ncrd);
        fputs(ptr,fout);
      }
      
      break;

    default:
      break;
  }
  
  return true;
} /* ReadSection */

static void PrintDim(qpdat* opt,
                     qppar* par)
{
  int  j,n,m,nz;
  char sl[_L_SIZE],sr[_L_SIZE];
    
  m =opt->m_intRow;
  n =opt->m_colUse;
  nz=opt->m_intAnz;
  
  for (j=0; j<n; j++)
  {
    if (opt->col[j].upp<=par->q_dblMax)
      opt->nbup++;
    if (opt->col[j].low>=-par->q_dblMax)
      opt->nblo++;
    if (opt->col[j].low<-par->q_dblMax&&opt->col[j].upp>par->q_dblMax)
      opt->nbfr++;
    if (fabs(opt->col[j].upp-opt->col[j].low)<par->q_dblAij)
      opt->nbfx++;
  }
  
  fprintf(fout,"\nProblem Dimension\n");
  fprintf(fout,"-----------------\n");
  printf(" ROWS\n");
  fprintf(fout," ROWS\n");
  
  sprintf(sl,"   total");
  sprintf(sr,"%d",m);
  LeftDots(sl,sr);
  
  sprintf(sl,"   equalities");
  sprintf(sr,"%d",m-opt->m_intUne);
  LeftDots(sl,sr);
  
  sprintf(sl,"   inequalities");
  sprintf(sr,"%d",opt->m_intUne);
  LeftDots(sl,sr);
  
  sprintf(sl,"   ranged");
  sprintf(sr,"%d",opt->m_intNrg);
  LeftDots(sl,sr);

  printf(" COLUMNS\n");
  fprintf(fout," COLUMNS\n");
  
  sprintf(sl,"   total");
  sprintf(sr,"%d",n);
  LeftDots(sl,sr);
  
  sprintf(sl,"   lower bounded");
  sprintf(sr,"%d",opt->nblo);
  LeftDots(sl,sr);
  
  sprintf(sl,"   upper bounded");
  sprintf(sr,"%d",opt->nbup);
  LeftDots(sl,sr);
  
  sprintf(sl,"   fixed");
  sprintf(sr,"%d",opt->nbfx);
  LeftDots(sl,sr);
  
  sprintf(sl,"   free");
  sprintf(sr,"%d",opt->nbfr);
  LeftDots(sl,sr);
  
  printf(" BOUND CARDS\n");
  fprintf(fout," BOUND CARDS\n");
  
  sprintf(sl,"   total");
  sprintf(sr,"%d",opt->m_intNlo
                 +opt->m_intNup
                 +opt->m_intNfx
                 +opt->m_intNfr
                 +opt->m_intNmi
                 +opt->m_intNpl);
  LeftDots(sl,sr);
  
  sprintf(sl,"   lower");
  sprintf(sr,"%d",opt->m_intNlo);
  LeftDots(sl,sr);
  
  sprintf(sl,"   upper");
  sprintf(sr,"%d",opt->m_intNup);
  LeftDots(sl,sr);
  
  sprintf(sl,"   fixed");
  sprintf(sr,"%d",opt->m_intNfx);
  LeftDots(sl,sr);
  
  sprintf(sl,"   free");
  sprintf(sr,"%d",opt->m_intNfr);
  LeftDots(sl,sr);
  
  sprintf(sl,"   minus infinity");
  sprintf(sr,"%d",opt->m_intNmi);
  LeftDots(sl,sr);
    
  sprintf(sl,"   plus infinity");
  sprintf(sr,"%d",opt->m_intNpl);
  LeftDots(sl,sr);
  
  printf(" DENSITY\n");
  fprintf(fout," DENSITY\n");
  
  sprintf(sl,"   matrix A");
  sprintf(sr," %.3f%s",
             100.0*(((double)nz/(double)m)/(double)n),
             "%");
  LeftDots(sl,sr);
  
  sprintf(sl,"   matrix Q");
  sprintf(sr," %.3f%s",
             100.0*(((double)(2*opt->m_intQnz+opt->ndia)
             /(double)n)/(double)n),
             "%");
  LeftDots(sl,sr);
  
  printf("\n");
  fprintf(fout,"\n");
} /* PrintDim */

int DataInput(qppar* par,
              qpdat* opt,
              char*  nmps)
{
  int    i,k,nCard,nErr,idKey;
  mpssec idSec;
  char   snam[128];
  FILE*  fp;

  printf(" BEGIN reading data file ...\n");  
  
  nErr        =0;
  nCard       =0;
  opt->mpscard=cAlloc(_L_SIZE);
  
  fmps=fopen(nmps,"r");
  if (!fmps)
  {
    sprintf(snam,"%s.mps",nmps);
    fmps=fopen(snam,"r");
    if (!fmps)
    {
      ShutDown(NOT_SCHFIL);
      return NOT_SCHFIL;
    }
  }
  
  fprintf(fout,"INPUT FILE\n");
  fprintf(fout,"----------\n");

  for (idSec=NAME; idSec<=ENDATA; ++idSec) {
    idKey=SecProc(idSec,par,opt,&nCard,&nErr);
    if (idKey!=true) {
      fclose(fmps);
      fclose(fout);
      return idKey;
    }  
  }

  fclose(fmps);
  
  opt->m_intNsk=opt->m_intUne+opt->m_intNrg;
  
  fprintf(fout,"\nxxxxxx    Total number of errors");
  fprintf(fout," in MPS file %d\n",nErr);
  printf(" END reading\n");
  
  PrintDim(opt,par);
  
  fp = fopen("coplfile.bin","wb");
  if (!fp)
  {
    ShutDown(NOT_DSKSPC);
    return NOT_DSKSPC;
  }
  fwrite(&opt->m_intRow,sizeof(int),1,fp);
  fwrite(&opt->m_intCol,sizeof(int),1,fp);
  fwrite(&opt->m_intAnz,sizeof(int),1,fp);

  for (i=0; i<opt->m_intRow; i++)
  {
    k=strlen(opt->rowname[i]);
    fwrite(&k,sizeof(int),1,fp);
    fwrite(opt->rowname[i],sizeof(char),k,fp);
  }
  for (i=0; i<opt->m_intCol; i++)
  {
    k=strlen(opt->colname[i]);
    fwrite(&k,sizeof(int),1,fp);
    fwrite(opt->colname[i],sizeof(char),k,fp);
  }
  fwrite(opt->r,sizeof(double),opt->m_intRow,fp);
  fwrite(opt->b,sizeof(double),opt->m_intRow,fp);
  fclose(fp);

  cFree(&par->q_strRhs);
  cFree(&par->q_strRng);
  cFree(&par->q_strBnd);
  cFree(&par->q_strObj);
    
  hFree(&opt->rowtab);
  hFree(&opt->coltab);
  
  /* RecNames(nmps);*/
  cFree(&opt->mpscard);
  cFree(&opt->prbname);
  cFree(&opt->objname);
  
  cPfre(&opt->rowname,opt->m_intRow);
  cPfre(&opt->colname,opt->m_intCol);

  return OPT_OK;
} /* DataInput */
