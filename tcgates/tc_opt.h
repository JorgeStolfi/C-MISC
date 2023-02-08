#ifndef tc_opt_H
#define tc_opt_H

/* Tools for test-and-complement gates and circuits. */
/* Last edited on 2008-02-12 02:06:00 by stolfi */

#define tc_opt_H_COPYRIGHT "Copyright © 2008 by the State University of Campinas (UNICAMP)"
#define tc_opt_H_AUTHORS "Created 2008-fev-11 by Jorge Stolfi, IC-UNICAMP"

#include <bool.h>

#include <stdio.h>
#include <stdint.h>
#include <values.h>

#include <tcgates.h>

typedef struct tc_opt_table_t tc_opt_table_t;
  /* A table of non-ident bifuns with their optimal circuits. */

tc_opt_table_t *tc_opt_table_build(int m, int n, bifun_ct_t max_bifuns, int max_gates);
  /* Tries to compute the optimal TC circuits for all bifuns from {B^m} to {B^n}
    However, stops after generating {max_bifuns} of them. */

void tc_opt_table_print_table(FILE *wr, tc_opt_table_t *tb);
  /* Prints all bifuns in the table {tb} and their optimal TCCs. Also
    prints a summary with the number of optimal TCCs by size: */


#endif
