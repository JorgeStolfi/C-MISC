/* Last edited on 2014-07-22 02:19:58 by stolfilocal */
#ifndef gch_problem_input_H
#define gch_problem_input_H

#define _GNU_SOURCE

#include <gch_vector_list.h>

typedef struct gch_problem_input_t
  { int count;
    int *realIndex;
  } gch_problem_input_t;

gch_problem_input_t* gch_problem_input_new(gch_vector_list_sub_t *vsl);

void gch_problem_input_free(gch_problem_input_t *pi);

int gch_problem_input_compare(void *sol, void *inp);


#endif
