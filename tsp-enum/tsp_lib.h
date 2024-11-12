/* Typtes and functions for TSP/CycleCover statistics */
/* Last edited on 2023-03-31 04:17:51 by stolfi */
  
#ifndef tsp_lib_H
#define tsp_lib_H

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

#include <bool.h>

typedef unsigned char vtx_t;
/* A {vtx_t} is a vertex number, starting from 0. */

typedef enum { 
    DISTR_GEN_USPHE=0, /* Uniform in unit sphere of {R^{nv × nv}}. */
    DISTR_GEN_UUNIT=1, /* Uniform in {[0,1]}. */
    DISTR_GEN_BIMOD=2, /* Bimodal in {[0,1]} with {1/d} of small values. */
    DISTR_GEN_GAUSS=3, /* Gaussian. */
    DISTR_EUC_UCUBE=4, /* Euclidean with point distr uniform in cube {[-1,+1]^d}. */
    DISTR_EUC_UBALL=5, /* Euclidean with point distr uniform in unit {d}-ball. */
    DISTR_EUC_USPHE=6, /* Euclidean with point distr uniform in unit {d-1}-sphere. */
    DISTR_EUC_GAUSS=7  /* Euclidean with points in {R^d}, gaussian coord distr. */
  } Distr_t;
  /* A {Distr_t} value specifies a probability distribution for the
    arc costs. The distributions {DISTR_EUC_UCUBE, DISTR_EUC_UBALL,
    DISTR_EUC_GAUSS} are defined indirectly: the vertices are mapped to
    random points of {R^d}, for some specified {d}, and the arc costs
    are set to the Euclidean distances between the endpoints. The
    distribution {DISTR_GEN_UUNIT} assigns to each arc a random cost
    uniformly distributed between 0 and 1. The distribution
    {DISTR_GEN_BIMOD} assigns uniform random small costs {1/d} of the
    arcs, and high costs to the remaining arcs. */

#define DISTR_MAX 7
/* Largest numeric value of a {Distr_t}. */

char *name_from_distr(Distr_t  distr);
Distr_t distr_from_name(char *name);
/* Conversion of {Distr_t} values to and from strings. The 
  names are like the value ids of the {Distr_t} type, without the
  {DISTR_} prefix. */

/* 
   PARAMETER I/O */

void read_params
  ( FILE *rd,        /* Input file. */
    int32_t *seed,       /* Seed for random number generator. */           
    int32_t *nv,         /* Number of vertices. */                         
    bool_t *tour,    /* {TRUE} means tours, {FALSE} means cycloids. */ 
    int32_t *npmax,      /* Max number of perms to generate. */            
    Distr_t *distr,  /* Distribution type. */                          
    int32_t *d           /* A parameter for the cost distribution. */              
  );
  /* Reads the parameters of the problem from file {rd},
   one per line, in the '{name} = {value}" format. */

void write_params
  ( FILE *wr,        /* Output file. */                                 
    int32_t seed,        /* Seed for random number generator. */           
    int32_t nv,          /* Number of vertices. */                         
    bool_t tour,     /* {TRUE} means tours, {FALSE} means cycloids. */ 
    int32_t npmax,       /* Max number of perms to generate. */            
    Distr_t distr,   /* Distribution type. */                          
    int32_t d            /* A parameter for the cost distribution. */              
  );
  /* Writes the parameters of the problem to file {wr},
   one per line, in the '{name} = {value}" format. */

#endif
