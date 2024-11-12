/* 
  Last edited on 2009-02-10 09:31:49 by stolfi
  
  Uncompress files that were compressed with "mlhpack".
  
  Created by stolfi on 95-06-10. 
  
  Usage: mlhunpk [separator] < infile > outfile
  
  This program reads from stdin a compressed
  vocabulary generate by the "mlhpack" compressor
  (q.v.), which looks like this:
  
      comentara#ead1#comentar
      9s#32#9
      9#33#9
      9do#0ams
      7áramos#0ead4#comentar
      9eis#35#9
      7aram#36#9
  
  The program undoes the compression, writing to stdout
  the original lines:
  
      comentara#ead1#comentar
      comentaras#ead2#comentar
      comentara#ead3#comentar
      comentarado#ams
      comentáramos#ead4#comentar
      comentáreis#ead5#comentar
      comentaram#ead6#comentar
*/

#include <mlhpack.h>
#include <assert.h>

static char sep;  /* Field separator */      

static char pc[NF][NC];  /* Fields from previous line */ 
static long pn[NF];      /* Their lengths */
static long pf;          /* Their number */

int main(int argc, char *argv[])
  { char c; 
    long kl;  /* Line index */
    long kf;  /* Field index in current line */
    long kc;  /* Byte index in current field */ 
    char *tc; /* points to pc[kf][kc] if it exists */
    bool_t expanding; /* true if next byte must be a prefix length  */
    bool_t done;      /* no more lines to preocess */
    
    setup_digit_tables();

    check(argc <= 2, "usage: mlhpack [separator] < infile > outfile"); 
    if (argc == 2)
      { check (strlen(argv[1]) == 1, "bad separator"); 
        sep = *(argv[1]);
      }
    else
      { sep = '#'; }
    check(digval[((byte)sep)] == -1, "invalid separator");
      
    kl = 0; pf = 0; 
    kf = 0; kc = 0; tc = &(pc[0][0]);
    expanding = false;
    done = false;
    while (!done)
      { /* 
          Invariant: we have already expanded and written out kl lines,
          kf complete fields, and kc characters of the next field.  
          
          Variable pf is the number of fields in the
          previous line (or 0 if we are still on the first line).  
          For any i < NF, element pn[i] is the length of field [i]:
            
            of the current line, if i < kf
            of the previous line, if kf <= i < pf,
          
          and is garbage otherwise.
          
          For any i < NF and j < NC, pc[i][j] is byte [j] of field [i]:
          
            of the current line, if 
              i < kf, or
              i = kf and j < kc;
          
            of the previous line, if 
              i = kf and and i < pf and kc <= j < pn[i], or
              kf < i < pf;
          
          and is garbage otherwise. 
          
          If the variable "expanding" is true, the next character
          to be read (if any) will be the digit code preceding a new field.
          In that case we must have kf < pf and kc = 0.
          If the digit code is not '0', then we must also have
          kf < NF, and the digit's value d must satisfy
          d < MIN(pf[kf], NC).  In that case, the supressed
          prefix is bytes pc[kf][0..d-1].
          
          Finally, whenver kf < NF and kc < NC, we have 
          tc = &(pc[kf,kc]).
        */
        c = getchar(); num_read++;
        if (c == EOF)
          { check((kf == 0) && (kc == 0), "missing final newline"); 
            done = true;
          }
        else if (expanding)
          { int d = digval[(byte)c];
            assert((kf < pf) & (kc == 0));
            if (d != 0)
              { check((kf < NF) && (d <= NC), "mlhpack version mismatch");
                check((d > 0) && (d <= pn[kf]), "invalid prefix length");
                while(kc < d) { putchar(*tc); kc++; tc++; }
              }
            expanding = false;
          }
        else if (c == sep)
          { /* Finish field: */
            if (kf < NF) pn[kf] = kc;
            /* Write field separator */
            putchar(sep); kf++;
            /* Start new field: */
            kc = 0; tc = &(pc[kf][0]);
            expanding = (kf < pf);
          }
        else if (c == '\n')
          { /* Finish field: */
            if (kf < NF) pn[kf] = kc;
            /* Finish line: */
            pf = kf + 1; putchar('\n'); kl++; 
            /* Start new line: */
            kf = 0;
            /* Start new field: */ 
            kc = 0; tc = &(pc[0][0]);
            expanding = true; 
          }
        else
          { /* Copy byte and save if possible: */
            putchar(c); 
            if ((kf < NF) && (kc < NC)) { (*tc) = c; tc++; }
            kc++;
          }
    }
    check (!ferror(stdin), "error while reading stdin");
   
    fprintf(stderr, "%ld lines processed\n", kl);
    return(0);
  } 
  
  
            
  
      
