/* Last edited on 2014-07-22 02:29:09 by stolfilocal */
#ifndef gch_computed_solutions_H
#define gch_computed_solutions_H

#define _GNU_SOURCE

#include <gem.h>

#include <gch_tree.h>
#include <gch_list.h>
#include <gch_problem_input.h>

typedef struct gch_computed_solutions_t
  { int dim;
    gch_tree_t **trees;
  } gch_computed_solutions_t;

gch_computed_solutions_t* gch_computed_solutions_new(int dim);

void gch_computed_solutions_free(gch_computed_solutions_t *cs);

void gch_computed_solutions_add(gch_computed_solutions_t *cs, int dim, gch_problem_input_t *pi, gem_ref_t gem, gch_list_t *facets);

gem_ref_t gch_solution_find(gch_computed_solutions_t *cs, int dim, gch_problem_input_t *pi, gch_list_t *facets);

#endif
