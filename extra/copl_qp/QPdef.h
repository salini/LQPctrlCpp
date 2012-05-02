#include "QPcnt.h"

typedef struct {
    int    q_intCol;
    int    q_intItr;
    int    q_intCac;
    int    q_intMin;
    
    char*  q_strRhs;
    char*  q_strRng;
    char*  q_strBnd;
    char*  q_strObj;
    
    double q_dblLow;
    double q_dblUpp;
    double q_dblMax;
    double q_dblAij;    
    double q_dblPiv;
    double q_dblGap;
    double q_dblPri;
    double q_dblDua;
    double q_dblStp;
    double q_dblRad;
    double q_dblIni;
  } qppar;

typedef struct list* plst;

typedef struct list {
    int   addr;
    char* key;
    plst  next;
  } hashlst;

typedef struct {
    int   lssz;
    plst* list;
  } hashtab;
  
typedef struct {
    double  low;
    double  upp;
    double  rhs;
    int     nn0;
    int*    ja;
    double* an;
  } optrow;
  
typedef struct {
    int     nn0;
    int*    ja;
    double* an;
  } array;

typedef struct {
    int    mnrs;
    int    nrow;
    int    mems;
    array  *ia;
  } matrix;

typedef struct {
    int     ncol;
    int     nnzo;
    double* diag;
    array*  cols;
  } symatx;

typedef struct {
    int     mrow;     /* number of rows allocated                    */
    int     nrow;     /* number of rows used                         */
    
    int     snnz;     /* number of indices for nonzeros in S         */
    int*    shead;    /* position of first nonzero in row i of S     */
    int*    ssize;    /* number of non-zeros in row i of S below     */
                      /* the diagonal                                */
    int*    ssub;     /* column index buffer for non-zeros in S      */
    double* diag;     /* diagonal matrix D in the factorization      */
    
    int     unnz;     /* number of nonzeros in the upper factor      */
    int     ujnz;     /* number of column indices in the compressed  */
                      /* indices buffer ujsub                        */ 
    int*    ujbeg;    /* beginning position of indices in row i of U */
    int*    uhead;    /* position of first nonzero in row i of U     */
    int*    ujsze;    /* number of indices in row i of U             */ 
    int*    usub;     /* compressed column index buffer of U         */
    double* uval;     /* nonzero values in factor U                  */
    
    int*    perm;     /* permutation order                           */
    int*    invp;     /* inverse order of perm                       */
    
    int     nsnds;    /* number of supernodes                        */
    int*    subg;     /* index of the first column in supernode i    */
    int     ndens;    /* numer of dense rows                         */
    int     nsndn;    /* number supernodes in dense rows             */
    int*    dhead;    /* pointer first column in each dense row      */
    int*    dsub;     /* indices in dense rows                       */
    int*    dbeg;     /* beginning of column index                   */ 
    int     sdens;    /* separate dense row                          */
    int     upst;     /* specified index in uval for st-cut problem  */
  } chol;

typedef struct {
    int      m_intSta;
    int      m_intCol;
    int      m_intRow;
    int      m_intAnz;
    int      m_intQnz;
    int      m_intNlo;
    int      m_intNup;
    int      m_intNfr;
    int      m_intNfx;
    int      m_intNmi;
    int      m_intNpl;
    int      m_intNsk;
    int      m_intNrg;
    int      m_intUne;
    int      m_colUse;
    
    int      ndia;
    int      nblo;
    int      nbup;
    int      nbfx;
    int      nbfr;
    
    hashtab* rowtab;
    hashtab* coltab;
    
    optrow*  col;
    
    char*    mpscard;
    char*    prbname;
    char*    objname;
    char*    rowtype;
    char**   rowname;
    char**   colname;
    
    int*     imem;      /* integer working memory       */
    int*     bid;       /* bound identification vactor  */
    int      items;     /* number of reductions         */
    int      bytes;     /* size of preprocessing file   */
    int      iter;      
    
    double*  rmem;      /* double working memory        */
    double*  b;
    double*  c;
    double*  l;
    double*  u;
    double*  r;
    double   obj;
        
    double   pobj;
    double   dobj;
    double   rgap;
    double   cgap;
    double   mu;
    double   mua;
    double   mu0;
    double   inf;
    double   rval;
    double   rpnm;
    double   rdnm;
    double   pinf;
    double   dinf;

    double   cx;
    double   by;
    double*  rp;
    double*  rd;
    double*  rl;
    double*  ru;
    double*  rr;
    double*  ry;
    double*  rf;
    double*  rh;
    double*  rz;
    double*  rw;
    double*  rs;
    double*  rq;
    double*  rk;
    double*  ro;
    double*  qx;
    

    double*  diag;
    double*  diau;
    double*  diaf;
    double*  dial;
    
    double*  x;
    double*  y;
    double*  g;
    double*  z;
    double*  t;
    double*  s;
    double*  v;
    double*  w;
    double*  p;
    double*  q;
    double*  f;
    double*  h;
    double*  k;
    
    double*  dx;
    double*  dy;
    double*  dg;
    double*  dz;
    double*  dt;
    double*  ds;
    double*  dv;
    double*  dw;
    double*  dp;
    double*  dq;
    double*  df;
    double*  dh;
    double*  dk;
    double*  dl;

    matrix*  a;
    matrix*  at;
    symatx*  qh;
    chol*    cl;
    
    qppar    par;
    char     q_strMsg[128];
    char     q_strNam[128];
  } qpdat;

typedef struct {
    int idep;
    int last;
    int most;
    int cure;
    int loca;
    int lowp;
    int ntot;
    
    int *head;
    int *port;
    int *fwrd;
    int *bwrd;
  } xlist;
  
typedef struct {
    int nnod;
    int nn0;
    int raft;
    int head;
    int last;
    int ntot;
    
    int *adjn;
    int *rbeg;
    int *rexs;
    int *rlen;
    int *rend;
    int *pres;
    int *succ;
  } order;
  
#include "QPfun.h"

extern FILE *fmps,
            *fout,
            *fspc,
            *fres;

extern char *nmps,
            *nout,
            *nspc;

extern qpdat* qt;
extern qppar* pa;

extern clock_t t_start;
extern clock_t t_end;
