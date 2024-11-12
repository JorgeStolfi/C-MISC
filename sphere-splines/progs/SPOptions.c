/* See SPOptions.h. */
/* Last edited on 2008-05-24 12:27:48 by stolfi */

/* Based on SPOptions.m3 by J.Stolfi, DEC-SRC, 1988.  */

#include <SPOptions.h>
#include <SPBasic.h>

#include <vec.h>
#include <r3.h>
#include <r4.h>
#include <bool.h>
#include <nat.h>
#include <affirm.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <values.h>

void SPOptions_Halt(SPOptions_Parser_t *pp);
void SPOptions_CheckAllParsed(SPOptions_Parser_t *pp, nat_t num);

SPOptions_Parser_t *SPOptions_NewParser(FILE *wr, int argn, char **argc)
  { SPOptions_Parser_t *pp = notnull(malloc(sizeof(SPOptions_Parser_t)), "no mem");
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

void SPOptions_SetUsage(SPOptions_Parser_t *pp, char *msg)
  { string_vec_expand(&(pp->usage), pp->nusage);
    pp->usage.e[pp->nusage] = msg;
    pp->nusage++;
  }

void SPOptions_Error(SPOptions_Parser_t *pp, char *msg)
  { fprintf(stderr, "%s\n", msg);
    SPOptions_Halt(pp);
  }

void SPOptions_Halt(SPOptions_Parser_t *pp)
  { int i;
    fprintf(pp->wr, "usage: ");
    for (i = 0; i < pp->nusage; i++)
      { fprintf(pp->wr, "%s", pp->usage.e[i]); }
    fprintf(pp->wr, "\n");
    exit(1);
  }

bool_t SPOptions_IsNext(SPOptions_Parser_t *pp, char *key)
  { char **a = pp->arg.e;
    bool_t *p = pp->parsed.e;
    if (pp->next >= pp->arg.ne) { return FALSE; }
    return (! p[pp->next]) && (strcmp(key, a[pp->next]) == 0);
  }

bool_t SPOptions_TestKeywordGen(SPOptions_Parser_t *pp, char *key, bool_t local)
  { if (local)
      { return SPOptions_TestKeywordNext(pp, key); }
    else
      { return SPOptions_TestKeyword(pp, key); }
  }

bool_t SPOptions_TestKeywordNext(SPOptions_Parser_t *pp, char *key)
  { bool_t *p = pp->parsed.e;
    if (SPOptions_IsNext(pp, key))
      { p[pp->next] = TRUE;
        pp->next++;
        return TRUE;
      }
    else
      { return FALSE; }
  }
  
bool_t SPOptions_TestKeyword(SPOptions_Parser_t *pp, char *key)
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

void SPOptions_GetKeywordGen(SPOptions_Parser_t *pp, char *key, bool_t local)
  { if (local)
      { SPOptions_GetKeywordNext(pp, key); }
    else
      { SPOptions_GetKeyword(pp, key); }
  }

void SPOptions_GetKeywordNext(SPOptions_Parser_t *pp, char *key)
  { bool_t *p = pp->parsed.e;
    if (SPOptions_IsNext(pp, key))
      { p[pp->next] = TRUE;
        pp->next++;
      }
    else
      { fprintf(stderr, "%s: parameter %d = \"%s\" should be \"%s\".\n",
          pp->arg.e[0], pp->next, pp->arg.e[pp->next], key
        );
        SPOptions_Halt(pp);
      }
  }

void SPOptions_GetKeyword(SPOptions_Parser_t *pp, char *key)
  {
    if (! SPOptions_TestKeyword(pp, key))
      { fprintf(pp->wr, "%s: keyword \"%s\" not found.\n", pp->arg.e[0], key);
        SPOptions_Halt(pp);
      }
  }

char *SPOptions_GetNext(SPOptions_Parser_t *pp)
  { char **a = pp->arg.e;
    bool_t *p = pp->parsed.e;
    if ((pp->next >= pp->arg.ne) || (p[pp->next]))
      { fprintf(pp->wr, "%s: missing argument after argument %d = \"%s\".\n", 
          a[0], pp->next-1, a[pp->next-1]
        );
        SPOptions_Halt(pp); 
      }
    p[pp->next] = TRUE;
    pp->next++;
    return a[pp->next-1];
  }

int SPOptions_GetNextInt(SPOptions_Parser_t *pp, int min, int max) 
  { long int nn;
    char *rest;
    char *txt = SPOptions_GetNext(pp);
    nn = strtol(txt, &rest, 10);
    if ((*rest) != '\000')
      { fprintf(stderr, "%s: parameter %d = \"%s\" should be an integer.\n",
          pp->arg.e[0], pp->next-1, txt
        );
        SPOptions_Halt(pp); 
      }
    if ((nn < min) || (nn > max))
      { fprintf(stderr, "%s: parameter %d = \"%s\" should be in [%d..%d]\n",
          pp->arg.e[0], pp->next-1, txt, min, max
        );
        SPOptions_Halt(pp); 
      }
    return nn;
  }

double SPOptions_GetNextDouble(SPOptions_Parser_t *pp, double min, double max)
  { double x;
    char *txt = SPOptions_GetNext(pp);
    char *rest;
    x = strtod(txt, &rest);
    if ((*rest) != '\000')
      { fprintf(stderr, "%s: parameter %d = \"%s\" should be a number.\n",
          pp->arg.e[0], pp->next-1, txt
        );
        SPOptions_Halt(pp); 
      }
    if ((x < min) || (x > max))
      { fprintf(stderr, "%s: parameter %d = \"%s\" should be in [%g _ %g]\n",
          pp->arg.e[0], pp->next-1, txt, min, max
        );
        SPOptions_Halt(pp); 
      }
    return x;
  }

r3_t SPOptions_GetNextR3(SPOptions_Parser_t *pp, double min, double max)
  { int i;
    r3_t p;
    for (i = 0; i < 3; i++)
      { p.c[i] = SPOptions_GetNextDouble(pp, min, max); }
    return p;
  } 

r4_t SPOptions_GetNextR4(SPOptions_Parser_t *pp, double min, double max)
  { int i;
    r4_t p;
    for (i = 0; i < 4; i++)
      { p.c[i] = SPOptions_GetNextDouble(pp, min, max); }
    return p;
  } 
  
r3_t SPOptions_GetNextDir(SPOptions_Parser_t *pp)
  { int i;
    r3_t d;
    for (i = 0; i < 3; i++)
      { d.c[i] = SPOptions_GetNextDouble(pp, -DBL_MAX, DBL_MAX); }
    r3_dir(&d, &d);
    return d;
  } 

int_vec_t SPOptions_GetIntList(SPOptions_Parser_t *pp, char *key, int min, int max)
  { int_vec_t a = int_vec_new(10);
    int nInts = 0;
    while (SPOptions_TestKeyword(pp, key))
      { int_vec_expand(&a, nInts);
        a.e[nInts] = SPOptions_GetNextInt(pp, min, max);
        nInts++;
      }
    int_vec_trim(&a, nInts);
    return a;
  }

void SPOptions_CheckAllParsed(SPOptions_Parser_t *pp, nat_t num)
  { int MaxBogus = 5;  /* Give up after this many leftovers. */
    nat_t bogus = 0;
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
    if (bogus > 0) { SPOptions_Halt(pp); }
  }

void SPOptions_SkipParsed(SPOptions_Parser_t *pp)
  { bool_t *p = pp->parsed.e;
    pp->next = pp->arg.ne;
    while ((pp->next > 0) && (! p[pp->next-1])) { pp->next--; }
    /* Check for unparsed arguments: */
    SPOptions_CheckAllParsed(pp, pp->next);
  }

void SPOptions_Finish(SPOptions_Parser_t *pp)
  { SPOptions_CheckAllParsed(pp, pp->arg.ne);
  }

