#define PROG_NAME "unitexture"
#define PROG_DESC "generate a unidimensional texture by iterated local map"
#define PROG_VERS "1.0"

/* Copyright © 2007 by the State University of Campinas (UNICAMP). */
/* See the copyright, authorship, and warranty notice at end of file. */
/* Last edited on 2009-02-09 20:46:33 by stolfi */

/* TO DO:
 + Output the texture (as dynamic plot or static image).
*/

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    -size {N} \\\n" \
  "    -iterations {M} \\\n" \
  "    " argparser_help_info_HELP ""

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  ??.\n" \
  "\n" \
  "OPTIONS\n" \
  "  -size {N}\n" \
  "    This mandatory argument specifies the number of elements" \
  " in the state vector, which must be even.\n" \
  "\n" \
  "  -iterations {M}\n" \
  "    This mandatory argument specifies the number of iterations to perform.\n" \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  ?? (1).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Started around sep/1994 by Jorge Stolfi, UNICAMP, perhaps inspired" \
  " by Quartucci Forster's first undergraduate project on texture synthesis.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  1994-09-28 First draft of {iterate} function.\n" \
  "  2007-08-15 Completed by J.Stolfi into a working program.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  Copyright © 2007 by the State University of Campinas (UNICAMP).\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#define _GNU_SOURCE

#include <argparser.h>
#include <affirm.h>
#include <jsrandom.h>

#include <math.h>
#include <stdio.h>

/* TYPES */

typedef struct utx_options_t 
  { int size;        /* Size of state vector. */
    int iterations;  /* Number of iterations to perform. */
  } utx_options_t;

typedef float op_t(float fm, float fo, float fp, int i);
  /* Type of a local operator as required by {iterate}. The arguments
    are the values {fm,fo,fp} of three consecutive elements of the
    state vector, at positions {i-1,i,i+1}, respectively (with indices
    taken modulo {n}). The function must return the value 
    that should replace {fo} as element {i}. */

utx_options_t *utx_get_options(int argc, char**argv);
  /* Parses the command line options, packs 
    them into a {utx_options_t} record. */

void iterate(op_t *op, float f[], int n);
  /* Applies the local operator {op} toall elements {f[0..n-1]}, in
    staggered fashion. The length {n} must be even, and the array is
    assumed to be cyclic ({f[i+n] = f[i]} for all {i}). */

void compress_range(float f[], int n);
  /* Compresses the range {[-oo _ +oo]} to {[0 _ 1]}. */

void expand_range(float f[], int n);
  /* Expands the range {[0 _ 1]} to {[-oo _ +oo]}
    (actually, to {[-BIG _ +BIG]}). */

float gaussop(float fm, float fo, float fp, int i);
  /* A local operator suitable for {iterate}. */

/* PROTOTYPES */

int main(int argc, char**argv);
utx_options_t *utx_get_options(int argc, char**argv);

void initialize(float f[], int n);
  /* Sets the array {f[0..n-1]} to its initial state. */

/* IMPLEMENTATIONS */

int main(int argc, char**argv)
  {
    utx_options_t *o = utx_get_options(argc, argv);
    
    int n = o->size;
    float f[n];
    
    initialize(f, n);
    int k;
    for (k = 0; k < o->iterations; k++)
      { iterate(&gaussop, f, n);
        if (k % 10 == 0)
          { int i;
            fprintf(stderr, "iteration %5d\n", k);
            for (i = 0; i < n; i++)
              { fprintf(stderr, "  %3d %12.6f\n", i, f[i]); }
            fprintf(stderr, "\n");
          }
      }

    fprintf(stderr, "\n\ndone.\n");
    return 0;
  }

utx_options_t *utx_get_options(int argc, char**argv)
  { 
    utx_options_t *o = (utx_options_t *)notnull(malloc(sizeof(utx_options_t)), "no mem");
    
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    argparser_get_keyword(pp, "-size");
    o->size = argparser_get_next_int(pp, 0, 4096);
    if (o->size % 2 != 0) { argparser_error(pp, "size must be even"); }
    
    argparser_get_keyword(pp, "-iterations");
    o->iterations = argparser_get_next_int(pp, 0, 100000000);
    
    argparser_skip_parsed(pp);
    
    argparser_finish(pp);
    
    return o;
  }

void initialize(float f[], int n)
  {
    int i;
    for (i = 0; i < n; i++) { f[i] = 0; }
  }

void iterate(op_t *op, float f[], int n)
  {
    demand(n % 2 == 0, "{n} must be even");
    int i0, i;
    float fm, fo, fp;
    
    for (i0 = 0; i0 < 2; i0++)
      { for (i = i0; i < n; i += 2)
          {
            if (i == 0) fm = f[n-1]; else fm = f[i-1];
            fo = f[i];
            if (i == n-1) fp = f[0]; else fp = f[i+1];
            f[i] = op(fm, fo, fp, i);
          } 
      }
  }
  
void compress_range(float f[], int n)
  /* Compresses the range {[-oo _ +oo]} to {[0 _ 1]}. */
  {
    int i;
    for (i = 0; i < n; i++)
      { float fi = f[i];
        fi = atan(fi)/M_PI + 0.5;
        if (fi < 0) { fi = 0; }
        if (fi > 1) { fi = 1; }
        f[i] = fi;
      }
  }

#define BIG  1.0e+20
  /* Max absolute value in expanded range. */

#define GOP_SDEV 1.0
  /* Standard deviation of noise for {gaussop}. */

void expand_range(float f[], int n)
  /* Expands the range {[0 _ 1]} to {[-oo _ +oo]} (actually, to {[-BIG _ +BIG]}). */
  {
    int i;
    for (i = 0; i < n; i++)
      { float fi = (f[i] - 0.5)*M_PI;
        fi = (fi <= 0 ? -BIG : (fi >= M_PI/2 ? +BIG : tan(fi)));
        if (fi < -BIG) { fi = -BIG; }
        if (fi > +BIG) { fi = +BIG; }
        f[i] = fi;
      }
  }

float gaussop(float fm, float fo, float fp, int i)
  {
    /* Choose a neighbor {fr} to imitate: */
    int p = random();
    float fr = (p % 16 < 8 ? fm : fp); 
    /* Add Gaussian perturbation: */
    fr += GOP_SDEV*fgaussrand();
    return fr;
  }
