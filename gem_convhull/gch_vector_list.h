/* Last edited on 2014-07-22 03:11:36 by stolfilocal */
#ifndef gch_vector_list_H
#define gch_vector_list_H

#include <stdio.h>
#include <gch_vector.h>
#include <gch_affine_map.h>

typedef struct gch_vector_list_t
  { int dim;
    int count;
    gch_vector_t **vectors;
  } gch_vector_list_t;

gch_vector_list_t* gch_vector_list_read(char *filename);

void gch_vector_list_free(gch_vector_list_t *vl);

typedef struct gch_vector_list_sub_t 
  { gch_vector_list_t *mainList;
    int count;
    gch_vector_t **vectors;
    int *realIndex;
  } gch_vector_list_sub_t;

gch_vector_list_sub_t* gch_vector_list_sub_complete_new(gch_vector_list_t *vl);
gch_vector_list_sub_t* gch_vector_list_sub_complete_proj_new(gch_vector_list_t *vl, gch_affine_map_t *bt);
gch_vector_list_sub_t* gch_vector_list_sub_proj_new(gch_vector_list_sub_t *vsl, int *listIndexes, int numIndexes, gch_affine_map_t *bt);
void gch_vector_list_sub_free(gch_vector_list_sub_t *vsl);
gch_vector_t* gch_vector_list_real_index(gch_vector_list_sub_t *vsl, int i);

#endif
