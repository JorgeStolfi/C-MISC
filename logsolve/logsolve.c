#define PROG_NAME "logsolve"
#define PROG_DESC "checks the solver routines and heuristics"
#define PROG_VERS "1.0"

#define logsolve_C_COPYRIGHT \
  "Copyright © 2012 by the State University of Campinas (UNICAMP)"

/* Last edited on 2013-03-17 22:22:00 by stolfilocal */

#define PROG_HELP \
  PROG_NAME " {OPTIONS}..." \

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  "  " PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  This program ..." \
  "\n" \
  argparser_help_info_HELP_INFO "\n" \
  "SEE ALSO\n" \
  "  Hear also. Smell also. Touch also.\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created on or before 2008-02-09 by J. Stolfi (IC-UNICAMP).\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  2008-02-09 (?) created.\n" \
  "  2012-12-15 [J.Stolfi] added equation combination and splitting.\n" \
  "  2012-12-17 [J.Stolfi] added various sorting and initial guessing heuristics.\n" \
  "  2012-12-19 [J.Stolfi] removed pseudorandom calls.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " logsolve_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <bool.h>
#include <sign.h>
#include <argparser.h>

#include <logsys.h>
#include <logsys_arith.h>

#define ONE64 ((uint64_t)1)
#define ALL64 ((uint64_t)(-1LL))

#define logsolve_nx_MAX 32
#define logsolve_ny_MAX 32
#define logsolve_nz_MAX 64

typedef struct logsolve_options_t
  { int nx;               /* Bit width of first multiplier input. */   
    int ny;               /* Bit width of second multiplier input. */
    uint32_t X;           /* Numeric value of first factor. */
    uint32_t Y;           /* Numeric value of second factor. */
    bool_t condense;      /* TRUE to condense the system before solving. */
    int heuristic;        /* Heuristic to use, or -1 if none. */
  } logsolve_options_t;

logsolve_options_t *logsolve_parse_options(int argc, char **argv);

int main(int argc, char **argv)
  {
    logsolve_options_t *o = logsolve_parse_options(argc, argv);
    
    int nx = o->nx;
    int ny = o->ny;
    int nz = nx + ny;
    
    /* Since we use 32/64 bit integers: */
    assert(nx <= logsolve_nx_MAX);
    assert(ny <= logsolve_ny_MAX);
    assert(nz <= logsolve_nz_MAX);
    
    /* Create the basic system: */
    fprintf(stderr, "--- creating a (%d × %d) --> %d multiplier circuit ---\n", nx, ny, nz);
    logsys_va_t *xb[nx]; /* {xb[0..nx-1]} are the bits of {X}. */
    logsys_va_t *yb[ny]; /* {yb[0..ny-1]} are the bits of {Y}. */
    logsys_va_t *zb[nz]; /* {zb[0..nz-1]} are the bits of {Z}. */
    logsys_t *S = logsys_build_mul(nx, ny, xb, yb, zb);
    logsys_print_system(stdout, "--- basic system -----------------", S, "----------------------------------");

    /* Fix the variables: */
    fprintf(stderr, "--- test case ---\n");
    bool_t x[nx]; logsys_get_bits(o->X, nx, x);
    logsys_print_assignments(stderr, "  X = ", nx, x, "\n");
    
    bool_t y[ny]; logsys_get_bits(o->Y, ny, y);
    logsys_print_assignments(stderr, "  Y = ", ny, y, "\n");
    
    uint64_t Z = ((uint64_t)o->X)*((uint64_t)o->Y);
    bool_t z[nz]; logsys_get_bits(Z, nz, z);
    logsys_print_assignments(stderr, "  Z = ", nz, z, "\n");

    /* Fix the output variables in the system itself: */
    fprintf(stderr, "--- fixing the output variables in the system ---\n");
    int k;
    for (k = 0; k < nz; k++)
      { logsys_op3_t op = (z[k] ? logsys_op1_UN0 : logsys_op1_ZR0);
        (void)logsys_add_equation(S, op, 1, &(zb[k]), NULL); 
      }
    logsys_print_system(stdout, "--- forced system ----------------", S, "----------------------------------");

    if (o->condense)
      { fprintf(stderr, "--- condensing and splitting equations ---\n"); 
        logsys_condense_and_split(S);
        logsys_print_system(stdout, "--- condensed system -------------", S, "----------------------------------");
      }
    
    /* Find a solution for the system: */
    fprintf(stderr, "--- trying to find X and Y given Z ---\n");
    { /* Gather all its variables: */
      int nv = 0;
      logsys_va_t **v = NULL;
      logsys_get_variables(S, &nv, &v);
  
      /* Ranges of variables for solving: */
      bool_t lo[nv], hi[nv]; 
      uint64_t r;
      for (r = 0; r < nv; r++) { lo[r] = FALSE; hi[r] = TRUE; }

      if (nv <= 20)
        { /* Enumerate all solutions exhaustively: */
          fprintf(stderr, "--- enumerating solutions exhaustively ---\n");
          
          auto bool_t printsol(int nv, logsys_va_t *va[], bool_t x[]);
          
          bool_t printsol(int nv, logsys_va_t *va[], bool_t x[])
            { fprintf(stderr, "solution:\n");
              logsys_print_var_ranges(stderr, "  X = ", nx, xb, x,x, "\n");
              logsys_print_var_ranges(stderr, "  Y = ", ny, yb, x,x, "\n");
              logsys_print_var_ranges(stderr, "  Z = ", nz, zb, x,x, "\n");
              fprintf(stderr, "\n");
              return TRUE;
            }
          
          logsys_enum_solutions(S, nv, v, lo, hi, &printsol);
        }


      /* Try to sort variables in a suitable order: */
      int *vix = NULL;      /* Variable search order. */
      bool_t *guess= NULL;  /* First guess for each variable. */
      if (o->heuristic >= 0)
        { fprintf(stderr, "sorting vars by heuristic %d\n", o->heuristic); 
          vix = notnull(malloc(nv*sizeof(int)), "no mem");
          guess = notnull(malloc(nv*sizeof(bool_t)), "no mem");
          logsys_sort_variables_for_solver(o->heuristic, S, nv, v, vix, guess);
        }
      else
        { fprintf(stderr, "not using any variable sorting heuristic\n"); }
      
      /* Look for a solution: */
      bool_t ok = logsys_find_solution(S, nv, v, vix, lo, hi, guess);
      if (ok)
        { /* Print and check it: */
          fprintf(stderr, "--- found a presumed solution ---\n");
          logsys_print_var_ranges(stderr, "X = ", nx, xb, lo, hi, "\n");
          logsys_print_var_ranges(stderr, "Y = ", ny, yb, lo, hi, "\n");
          logsys_print_var_ranges(stderr, "Z = ", nz, zb, lo, hi, "\n");
          fprintf(stderr, "--- variables in trial order ---\n");
          int nguess_OK = 0;
          for (r = 0; r <  nv; r++) 
            { int iu = (vix != NULL ? vix[r] : r);
              logsys_va_t *ua = v[iu];
              if (ua != NULL)
                { assert(logsys_variable_id(ua) == iu);
                  logsys_print_variable_id(stderr, "  v", iu, 4, NULL);
                  fprintf(stderr, " val = %c", "01"[lo[iu]]);
                  if (guess != NULL) 
                    { fprintf(stderr, " guess = %c", "01"[guess[iu]]); 
                      if (guess[iu] == lo[iu]) { nguess_OK++; } 
                    }
                  fprintf(stderr, "\n");
                  demand(lo[iu] == hi[iu], "presumed solution is still indeterminate");
                }
            }
          if (guess != NULL)
            { double pct_OK = ((double)nguess_OK)*100/((double)nv);
              fprintf(stderr, "  total %d correct guesses ( %5.1f%% )\n", nguess_OK, pct_OK);
            }
          fprintf(stderr, "checking solution ...");
          if (logsys_is_solution(S, nv, v, lo))
            { fprintf(stderr, " OK!\n"); }
          else
            { fprintf(stderr, " NOT!\n");
              demand(FALSE, "reported solution does not satisfy system");
            }
        }
      else
        { fprintf(stderr, "** no solution found!\n");
          assert(FALSE);
        }
    }

    return 0;
  }
    
logsolve_options_t *logsolve_parse_options(int argc, char **argv)
  {
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);

    logsolve_options_t *o = (logsolve_options_t *)notnull(malloc(sizeof(logsolve_options_t)), "no mem"); 
    
    /* Get keyword arguments: */

    /* Operator kind and parameters: */
    argparser_get_keyword(pp, "-multiplier");
    o->nx = argparser_get_next_int(pp, 1, logsolve_nx_MAX);
    o->ny = argparser_get_next_int(pp, 1, logsolve_ny_MAX);
    if (o->nx + o->ny > logsolve_nz_MAX)
      { argparser_error(pp, "max multiplier output size exceeded"); }

    argparser_get_keyword(pp, "-factors");
    uint64_t XMAX = (ONE64 << (uint64_t)o->nx) + ALL64;
    uint64_t YMAX = (ONE64 << (uint64_t)o->ny) + ALL64;
    o->X = argparser_get_next_uint(pp, 0, XMAX);
    o->Y = argparser_get_next_uint(pp, 0, YMAX);
    
    o->condense = argparser_keyword_present(pp, "-condense");
    
    if (argparser_keyword_present(pp, "-heuristic"))
      { o->heuristic = argparser_get_next_int(pp, -1, 100000); }
    else
      { /* Choose a good scale: */
        o->heuristic = -1;
      }
    
    /* Skip to positional arguments: */
    argparser_skip_parsed(pp);
    
    /* Check for spurious args: */
    argparser_finish(pp);

    return o;
  }
