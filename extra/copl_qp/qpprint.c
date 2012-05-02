#include "QPdef.h"

void LeftDots(char *sl,char *sr)
{
  char sz[75];
  int  i,k,len;
  
  k=strlen(sl);
  sprintf(sz,"%s",sl);
  i=strlen(sr);
  len=75-i;
  
  for (i=k; i<len; i++)
    sz[i]='.';
  sz[len]='\0';
  
  printf("%s%s\n",sz,sr);
  fprintf(fout,"%s%s\n",sz,sr);
} /* LeftDots */
