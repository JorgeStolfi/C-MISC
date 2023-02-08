/* See {gch_computed_solutions.h} */
/* Last edited on 2014-07-22 02:29:46 by stolfilocal */

#define _GNU_SOURCE
#include <stdlib.h>

#include <gem.h>

#include <gch_list.h>
#include <gch_problem_input.h>
#include <gch_solution.h>

#include <gch_computed_solutions.h>

gch_computed_solutions_t* gch_computed_solutions_new(int dim)
  { gch_computed_solutions_t* cs = (gch_computed_solutions_t*)malloc(sizeof(gch_computed_solutions_t));
    cs->dim = dim;
    cs->trees = (gch_tree_t**)calloc(dim-2,sizeof(gch_tree_t*));
    return cs;
  }

void gch_computed_solutions_free(gch_computed_solutions_t *cs)
  { int i;
    for (i = 0; i < cs->dim-2; i++)
      { gch_tree_gem_traverse(cs->trees[i],gch_solution_free);
        gch_tree_free(cs->trees[i]);
      }
    free(cs->trees);
    free(cs);
  }

gem_ref_t gch_solution_find(gch_computed_solutions_t *cs, int dim, gch_problem_input_t *pi, gch_list_t *facets)
  { gch_solution_t *solution = (gch_solution_t*)gch_tree_find(cs->trees[dim-1],(void*)pi,gch_problem_input_compare);
    if (solution != NULL)
      { return gch_solution_duplicate(solution, dim, facets); }
    return NULL;
  }

void gch_computed_solutions_add(gch_computed_solutions_t *cs, int dim, gch_problem_input_t *pi, gem_ref_t gem, gch_list_t *facets)
  { gch_solution_t *solution = (gch_solution_t*)malloc(sizeof(gch_solution_t));
    solution->pi = pi;
    solution->gem = gem;
    solution->facets = gch_list_new();
    gch_list_copy(solution->facets,facets);
    cs->trees[dim-1] = gch_tree_add(cs->trees[dim-1],solution,gch_solution_compare);
  }


