/* See SOParams.h. */
/* Last edited on 2007-01-04 00:22:40 by stolfi */

/* Based on SOParams.m3 by J.Stolfi, DEC-SRC, 1988.  */

#include <SOParams.h>
#include <SOBasic.h>
#include <vec.h>
#include <r3.h>
#include <r4.h>
#include <affirm.h>
#include <nat.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <values.h>

void SOParams_Halt(SOParams_T *pp);
void SOParams_CheckAllParsed(SOParams_T *pp, unsigned num);

SOParams_T *SOParams_NewT(FILE *wr, int argn, char **argc)
  { SOParams_T *pp = (SOParams_T *)notnull(malloc(sizeof(SOParams_T)), "no mem for SOParams_T");
    pp->wr = wr;
    pp->arg = (string_vec_t){argn, argc};
    pp->parsed = bool_vec_new(argn);
    pp->parsed.e[0] = TRUE;
    { int i; for (i = 1; i < argn; i++) { pp->parsed.e[i] = FALSE; } }
    pp->next = 1;
    pp->nusage = 0;
    pp->usage = string_vec_new(10);
    return pp;
  }

void SOParams_SetUsage(SOParams_T *pp, char *msg)
  { string_vec_expand(&(pp->usage), pp->nusage);
    pp->usage.e[pp->nusage] = msg;
    pp->nusage++;
  }

void SOParams_Error(SOParams_T *pp, char *msg)
  { fprintf(stderr, "%s\n", msg);
    SOParams_Halt(pp);
  }

void SOParams_Halt(SOParams_T *pp)
  { int i;
    fprintf(pp->wr, "usage: ");
    for (i = 0; i < pp->nusage; i++)
      { fprintf(pp->wr, "%s", pp->usage.e[i]); }
    fprintf(pp->wr, "\n");
    exit(1);
  }

bool_t SOParams_KeywordPresent(SOParams_T *pp, char *key)
  { char **a = pp->arg.e;
    bool_t *p = pp->parsed.e;
    int i;
    for (i = 0; i < pp->arg.ne; i++)
      { if ((! p[i]) && (strcmp(key, a[i]) == 0))
          { pp->next = i + 1;
            p[i] = TRUE;
            return TRUE;
          }
      }
    return FALSE;
  }

void SOParams_GetKeyword(SOParams_T *pp, char *key)
  {
    if (! SOParams_KeywordPresent(pp, key))
      { fprintf(pp->wr, "%s: keyword \"%s\" not found.\n", pp->arg.e[0], key);
        SOParams_Halt(pp);
      }
  }

char *SOParams_GetNext(SOParams_T *pp)
  { char **a = pp->arg.e;
    bool_t *p = pp->parsed.e;
    if ((pp->next >= pp->arg.ne) || (p[pp->next]))
      { fprintf(pp->wr, "%s: missing argument after argument %d = \"%s\".\n", 
          a[0], pp->next-1, a[pp->next-1]
        );
        SOParams_Halt(pp); 
      }
    p[pp->next] = TRUE;
    pp->next++;
    return a[pp->next-1];
  }

bool_t SOParams_IsNext(SOParams_T *pp, char *key)
  { char **a = pp->arg.e;
    bool_t *p = pp->parsed.e;
    if (pp->next >= pp->arg.ne) { return FALSE; }
    return (! p[pp->next]) && (strcmp(key, a[pp->next]) == 0);
  }

bool_t SOParams_TestNext(SOParams_T *pp, char *key)
  { bool_t *p = pp->parsed.e;
    if (SOParams_IsNext(pp, key))
      { p[pp->next] = TRUE;
        pp->next++;
        return TRUE;
      }
    else
      { return FALSE; }
  }

void SOParams_MatchNext(SOParams_T *pp, char *key)
  { bool_t *p = pp->parsed.e;
    if (SOParams_IsNext(pp, key))
      { p[pp->next] = TRUE;
        pp->next++;
      }
    else
      { fprintf(stderr, "%s: parameter %d = \"%s\" should be \"%s\".\n",
          pp->arg.e[0], pp->next, pp->arg.e[pp->next], key
        );
        SOParams_Halt(pp);
      }
  }

int SOParams_GetNextInt(SOParams_T *pp, int min, int max) 
  { long int nn;
    char *rest;
    char *txt = SOParams_GetNext(pp);
    nn = strtol(txt, &rest, 10);
    if ((*rest) != '\000')
      { fprintf(stderr, "%s: parameter %d = \"%s\" should be an integer.\n",
          pp->arg.e[0], pp->next-1, txt
        );
        SOParams_Halt(pp); 
      }
    if ((nn < min) || (nn > max))
      { fprintf(stderr, "%s: parameter %d = \"%s\" should be in [%d..%d]\n",
          pp->arg.e[0], pp->next-1, txt, min, max
        );
        SOParams_Halt(pp); 
      }
    return nn;
  }

double SOParams_GetNextDouble(SOParams_T *pp, double min, double max)
  { double x;
    char *txt = SOParams_GetNext(pp);
    char *rest;
    x = strtod(txt, &rest);
    if ((*rest) != '\000')
      { fprintf(stderr, "%s: parameter %d = \"%s\" should be a number.\n",
          pp->arg.e[0], pp->next-1, txt
        );
        SOParams_Halt(pp); 
      }
    if ((x < min) || (x > max))
      { fprintf(stderr, "%s: parameter %d = \"%s\" should be in [%g _ %g]\n",
          pp->arg.e[0], pp->next-1, txt, min, max
        );
        SOParams_Halt(pp); 
      }
    return x;
  }

r2_t SOParams_GetNextR2(SOParams_T *pp, double min, double max)
  { int i;
    r2_t p;
    for (i = 0; i < 2; i++)
      { p.c[i] = SOParams_GetNextDouble(pp, min, max); }
    return p;
  } 

r3_t SOParams_GetNextR3(SOParams_T *pp, double min, double max)
  { int i;
    r3_t p;
    for (i = 0; i < 3; i++)
      { p.c[i] = SOParams_GetNextDouble(pp, min, max); }
    return p;
  } 

r4_t SOParams_GetNextR4(SOParams_T *pp, double min, double max)
  { int i;
    r4_t p;
    for (i = 0; i < 4; i++)
      { p.c[i] = SOParams_GetNextDouble(pp, min, max); }
    return p;
  } 

void SOParams_GetNextRN(SOParams_T *pp, int N, double min, double max, double *x)
  { int i;
    for (i = 0; i < N; i++)
      { x[i] = SOParams_GetNextDouble(pp, min, max); }
  } 
  
r2_t SOParams_GetNextR2Dir(SOParams_T *pp)
  { int i;
    r2_t d;
    for (i = 0; i < 2; i++)
      { d.c[i] = SOParams_GetNextDouble(pp, -INFTY, INFTY); }
    r2_dir(&d, &d);
    return d;
  } 
  
r3_t SOParams_GetNextR3Dir(SOParams_T *pp)
  { int i;
    r3_t d;
    for (i = 0; i < 3; i++)
      { d.c[i] = SOParams_GetNextDouble(pp, -INFTY, INFTY); }
    r3_dir(&d, &d);
    return d;
  } 

int_vec_t SOParams_GetIntList(SOParams_T *pp, char *key, int min, int max)
  { int_vec_t a = int_vec_new(10);
    int nInts = 0;
    while (SOParams_KeywordPresent(pp, key))
      { int_vec_expand(&a, nInts);
        a.e[nInts] = SOParams_GetNextInt(pp, min, max);
        nInts++;
      }
    int_vec_trim(&a, nInts);
    return a;
  }

double_vec_t SOParams_GetDoubleList(SOParams_T *pp, char *key, double min, double max)
  { double_vec_t a = double_vec_new(10);
    int nDoubles = 0;
    while (SOParams_KeywordPresent(pp, key))
      { double_vec_expand(&a, nDoubles);
        a.e[nDoubles] = SOParams_GetNextDouble(pp, min, max);
        nDoubles++;
      }
    double_vec_trim(&a, nDoubles);
    return a;
  }

string_vec_t SOParams_GetList(SOParams_T *pp, char *key)
  { string_vec_t a = string_vec_new(10);
    int nStrings = 0;
    while (SOParams_KeywordPresent(pp, key))
      { string_vec_expand(&a, nStrings);
        a.e[nStrings] = SOParams_GetNext(pp);
        nStrings++;
      }
    string_vec_trim(&a, nStrings);
    return a;
  }

void SOParams_CheckAllParsed(SOParams_T *pp, unsigned num)
  { nat_t MaxBogus = 5;  /* Give up after this many leftovers. */
    unsigned bogus = 0;
    char **a = pp->arg.e;
    bool_t *p = pp->parsed.e;
    int i;
    for (i = 0; i < num; i++)
      { if (! p[i])
          { bogus++;
            if (bogus <= MaxBogus)
              { fprintf(pp->wr, "%s: parameter %d = \"%s\" extraneous or misplaced.\n",
                  a[0], i, a[i]
                );
              }
          }
      }
    if (bogus > MaxBogus) 
      { fprintf(pp->wr, "(and %d more).\n", bogus - MaxBogus); }
    if (bogus > 0) { SOParams_Halt(pp); }
  }

void SOParams_SkipParsed(SOParams_T *pp)
  { bool_t *p = pp->parsed.e;
    pp->next = pp->arg.ne;
    while ((pp->next > 0) && (! p[pp->next-1])) { pp->next--; }
    /* Check for unparsed arguments: */
    SOParams_CheckAllParsed(pp, pp->next);
  }

void SOParams_Finish(SOParams_T *pp)
  { SOParams_CheckAllParsed(pp, pp->arg.ne); }

