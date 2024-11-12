/* See tsp_lib.h */
/* Last edited on 2023-03-31 04:18:23 by stolfi */ 

#define _GNU_SOURCE
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>

#include <fget.h>
#include <nget.h>
#include <bool.h>
#include <affirm.h>

#include <tsp_lib.h>

#define DEBUG FALSE

#define PI M_PI
/* An obscure mathematical constant. */

static char* distr_name[DISTR_MAX+1] = 
  { "GEN_USPHE",
    "GEN_UUNIT",
    "GEN_BIMOD",
    "GEN_GAUSS",
    "EUC_UCUBE",
    "EUC_UBALL",
    "EUC_USPHE",
    "EUC_GAUSS"
  };

char *name_from_distr(Distr_t  distr)
  { return distr_name[distr]; }

Distr_t distr_from_name(char *name)
  { int32_t distr;
    for (distr = 0; distr <= DISTR_MAX; distr++)
      { if (strcmp(distr_name[distr], name) == 0) 
          { return (Distr_t)distr; }
      }
    demand(FALSE, "bad distr name");
    return (Distr_t)0;
  }

void read_params
  ( FILE *rd,        /* Input file. */
    int32_t *seed,       /* Seed for random number generator. */           
    int32_t *nv,         /* Number of vertices. */                         
    bool_t *tour,    /* {TRUE} means tours, {FALSE} means cycloids. */ 
    int32_t *npmax,         /* Max number of perms to generate. */            
    Distr_t *distr,  /* Distribution type. */                          
    int32_t *d           /* A parameter for the cost distribution. */              
  )
  { 
    (*seed)  = nget_int32(rd, "seed"); fget_eol(rd);
    (*nv)    = nget_int32(rd, "vertices"); fget_eol(rd);
    (*tour)  = nget_int32(rd, "tour"); fget_eol(rd);
    (*npmax) = nget_int32(rd, "max-perms"); fget_eol(rd);
    (*distr) = nget_int32(rd, "arc-cost-distr"); fget_eol(rd);
    (*d)     = nget_int32(rd, "arc-cost-param"); fget_eol(rd);
  }

void write_params
  ( FILE *wr,        /* Output file. */                                 
    int32_t seed,        /* Seed for random number generator. */           
    int32_t nv,          /* Number of vertices. */                         
    bool_t tour,     /* {TRUE} means tours, {FALSE} means cycloids. */ 
    int32_t npmax,       /* Max number of perms to generate. */            
    Distr_t distr,   /* Distribution type. */                          
    int32_t d            /* A parameter for the cost distribution. */              
  )
  { 
    fprintf(wr, "seed = %d\n", seed);
    fprintf(wr, "vertices = %d\n", nv);
    fprintf(wr, "tour = %d\n", tour);
    fprintf(wr, "max-perms = %d\n", npmax);
    fprintf(wr, "arc-cost-distr = %d\n", distr);
    fprintf(wr, "arc-cost-param = %d\n", d);
  }

