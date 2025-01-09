/* Last edited on 2024-12-25 09:56:26 by stolfi */
/* 

  Compress vocabulary files, squeezing out common prefixes.
  
  Created by stolfi on 95-06-09. 
  
  Usage: mlhpack [separator] < infile > outfile
  
  This program reads from stdin a sequence of lines,
  preferably sorted, where each line consists of 
  one or more fields delimited by some separator character
  (default '#').  For example,
  
      comentara#ead1#comentar
      comentaras#ead2#comentar
      comentara#ead3#comentar
      comentarado#ams
      comentáramos#ead4#comentar
      comentáreis#ead5#comentar
      comentaram#ead6#comentar
  
  The program writes to stout a compacted version of those
  same lines.  If the kth field of some line begins
  with the same n bytes as the kth field of the preceding line,
  then those bytes are replaced by the single digit n.
  If the previous line has no kth field, this compression is not 
  performed. 

  Note that the value of n is not necessarily maximal.
  Byte counts greater than 9 are encoded by letters:
  a=10, b=11, ..., up to z=35; so, if two fields 
  have a common prefix of more than 35 bytes, only the first 
  35 will be supressed (and replaced by "z"), the 
  rest being copied verbatim.  
  
  Also, if a line has more than NF fields,
  the program will always use n=0 when compressing 
  a field beyond the first NF. (Note that this 0 
  is inserted only if the field exists in the previous line).
  
  Thus, the file shown above would be compressed to
  
      comentara#ead1#comentar
      9s#32#9
      9#33#9
      9do#0ams
      7áramos#0ead4#comentar
      9eis#35#9
      7aram#36#9
  
  The separator must not be a letter or decimal digit.    
  Setting the separator to '\n' is allowed and has the effect
  of treting the entire line as a single field. (How to get
  the '\n' past the shell is YOUR problem.) 
      
AUTHOR
    J. Stolfi, DCC/Unicamp, june 1995.
  
*/

#include <mlhpack.h>

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
    bool_t packing; /* true if shared fiels prefix is being supressed */
    bool_t done;    /* no more lines to process */
    
    setup_digit_tables();
      
    check(argc <= 2, "usage: mlhpack [separator] < infile > outfile"); 
    if (argc == 2)
      { check (strlen(argv[1]) == 1, "separator must be a single character"); 
        sep = *(argv[1]);
       }
    else
      { sep = '#'; }  
    check(digval[((byte)sep)] == -1, "invalid separator");
      
    kl = 0; pf = 0; 
    kf = 0; kc = 0; tc = &(pc[0][0]);
    packing = false;
    done = false;
    while (!done)
      { /* 
          Invariant: we have already read and processed kl lines,
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
            
          The variable "packing" is true only when kf < pf;
          and then if and only if either kc = 0, 
          or kf < NF, kc <= MIN(pn[kf], NC), and 
          bytes [0..kc-1] of the current field [kf] match 
          bytes [0..kc-1] of field [kf] of the previous line. 
        
          Finally, whenver kf < NF and kc < NC, we have 
          tc = &(pc[kf,kc]).
        */
        c = getchar(); num_read++;
        if ((c == EOF) && (kf == 0) && (kc == 0))  
          { done = true; }
        else if ((c == '\n') || (c == EOF))
          { /* Finish field: */
            if (packing) putchar(digit[kc]); 
            if (kf < NF) pn[kf] = kc;
            /* Finish line: */
            pf = kf + 1; putchar('\n'); kl++; 
            /* Start new line: */
            kf = 0;
            /* Start new field: */ 
            kc = 0; tc = &(pc[0][0]);
            packing = true; 
          }
        else if (c == sep)
          { /* Finish field: */
            if (packing) putchar(digit[kc]); 
            if (kf < NF) pn[kf] = kc;
            /* Write field separator */
            putchar(sep); kf++;
            /* Start new field: */
            kc = 0; tc = &(pc[kf][0]);
            packing = (kf < pf);
          }
        else if ((kf >= NF) || (kc >= NC))
          { /* This byte can't be supressed or stored: */
            if (packing) { putchar(digit[kc]); packing = false; }
            putchar(c); kc++;
          }
        else 
          { /* Maybe we can supress this character: */
            if (packing)
              { packing = ((kf < pf) && (kc < pn[kf]) && (c == (*tc)));
                if (!packing) putchar(digit[kc]);
              } 
            /* I we can't, we must write and store it: */
            if (! packing) { putchar(c); (*tc) = c; }
            kc++; tc++;
          }
    }
    check (!ferror(stdin), "error while reading stdin");
   
    fprintf(stderr, "%ld lines processed\n", kl);
    return(0);
  } 
            
  
      
