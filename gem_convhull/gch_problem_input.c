/* See {gch_problem_input.h} */
/* Last edited on 2014-07-22 03:12:39 by stolfilocal */

#define _GNU_SOURCE
#include <stdlib.h>

#include <gem.h>

#include <gch_list.h>

#include <gch_solution.h>
#include <gch_problem_input.h>

gch_problem_input_t* gch_problem_input_new(gch_vector_list_sub_t *vsl)
  {
    int i;
    gch_problem_input_t* pi = (gch_problem_input_t*)malloc(sizeof(gch_problem_input_t));
    pi->count = vsl->count;
    pi->realIndex = (int*)malloc(sizeof(int)*(pi->count));
    for (i = 0; i < pi->count; i++) { pi->realIndex[i] = vsl->realIndex[i]; }
    return pi;
  }

void gch_problem_input_free(gch_problem_input_t *pi)
  { free(pi->realIndex);
    free(pi);
  }

int gch_problem_input_compare(void *sol, void *inp)
  {
    gch_solution_t *solution = (gch_solution_t*)sol;
    gch_problem_input_t *input = (gch_problem_input_t*)inp;
    int i;
    for (i = 0; (i < solution->pi->count)&&(i < input->count); i++)
      { if (solution->pi->realIndex[i] < input->realIndex[i]) { return -1; }
        if (solution->pi->realIndex[i] > input->realIndex[i]) { return +1; }
      }
    return 0;
  }

