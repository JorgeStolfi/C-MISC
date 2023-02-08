/* Last edited on 2014-07-22 03:04:28 by stolfilocal */
#ifndef gch_solution_H
#define gch_solution_H

#define _GNU_SOURCE

#include <gem.h>

#include <gch_list.h>
#include <gch_problem_input.h>

typedef struct gch_solution_t 
  { gch_problem_input_t *pi;
    gem_ref_t gem;
    gch_list_t *facets;
  } gch_solution_t;

gem_ref_t gch_solution_duplicate(gch_solution_t *sol, int dim, gch_list_t *facets);

void gch_solution_free(void *sol);

int gch_solution_compare(void *a, void *b);

#endif
