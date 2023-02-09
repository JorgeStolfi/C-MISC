
void WriteMatrix(FILE *f, double **d, long *it, long N, byte **L, long skip)
  {
    long i, j;
    /* Write header: */
    for(i=0;i<skip;i++)
      { for (j=0;j<skip;j++) putc(' ', f);
        putc(' ', f);
        for (j=0;j<N;j++) { putc(' ', f); putc(' ', f); putc(L[it[j]][i], f); }
        putc('\n', f);
      }
    for (j=0;j<skip;j++) putc(' ', f);
    putc(' ', f);
    for (j=0;j<N;j++) { putc(' ', f); putc('-', f); putc('-', f); }
    putc('\n', f);
    
    /* Write body: */
    for(i=0;i<N;i++)
      { for (j=0;j<skip;j++) putc(L[it[i]][j], f);
        putc(' ', f);
        for (j=0;j<N;j++) 
          { fprintf(f, " %2d", FormatDist(d[it[i]][it[j]])); }
        putc('\n', f);
      }
    fflush(f);
  }
  
int FormatDist(double d)
  {
    return (int)(99.4999*(d<1?d:1.0) + 0.5);
  }
    
