/* See {gch_vector_list.h} */
/* Last edited on 2014-07-22 02:28:41 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <gmp.h>

#include <gch_util.h>
#include <gch_affine_map.h>

#include <gch_vector_list.h>

gch_vector_list_t* gch_vector_list_read(char *filename)
  { gch_vector_list_t *vl = (gch_vector_list_t*)malloc(sizeof(gch_vector_list_t));
    int i,j;
    FILE *fp = fopen(filename, "r");
    int nr = fscanf(fp, "%d %d", &(vl->dim), &(vl->count));
    assert(nr == 2);
    vl->vectors = (gch_vector_t**)malloc(sizeof(gch_vector_t*)*(vl->count));
    for (i = 0; i < (vl->count); i++)
      { vl->vectors[i] = gch_vector_new((vl->dim));
        for (j = 0; j < (vl->dim); j++)
          { mpz_inp_str (vl->vectors[i]->coord[j], fp, 10); }
      }
    fclose(fp);
    return vl;
  }

void gch_vector_list_free(gch_vector_list_t *vl)
  { int i;
    for (i = 0; i < vl->count; i++) { gch_vector_free(vl->vectors[i]); }
    free(vl->vectors);
    free(vl);
  }

gch_vector_list_sub_t* gch_vector_list_sub_complete_new(gch_vector_list_t *vl)
  { gch_vector_list_sub_t *vsl = (gch_vector_list_sub_t*)malloc(sizeof(gch_vector_list_sub_t));
    vsl->count = vl->count;
    vsl->mainList = vl;
    vsl->realIndex = (int*)malloc(sizeof(int)*(vsl->count));
    vsl->vectors = (gch_vector_t**)malloc(sizeof(gch_vector_t*)*(vsl->count));
    int i;
    for (i = 0; i < vsl->count; i++)
      { vsl->realIndex[i] = i;
        vsl->vectors[i] = gch_vector_copy(vl->vectors[i]);
      }
    return vsl;
  }

gch_vector_list_sub_t* gch_vector_list_sub_complete_proj_new(gch_vector_list_t *vl, gch_affine_map_t *bt)
  { gch_vector_list_sub_t *vsl = (gch_vector_list_sub_t*)malloc(sizeof(gch_vector_list_sub_t));
    vsl->count = vl->count;
    vsl->mainList = vl;
    vsl->realIndex = (int*)malloc(sizeof(int)*(vsl->count));
    vsl->vectors = (gch_vector_t**)malloc(sizeof(gch_vector_t*)*(vsl->count));
    int i;
    for (i = 0; i < vsl->count; i++)
      { vsl->realIndex[i] = i;
        vsl->vectors[i] = gch_affine_map_map(bt,vl->vectors[i]);
      }
    return vsl;
  }

gch_vector_list_sub_t* gch_vector_list_sub_proj_new(gch_vector_list_sub_t *vsl, int *listIndexes, int numIndexes, gch_affine_map_t *bt)
  { gch_vector_list_sub_t *nvsl = (gch_vector_list_sub_t*)malloc(sizeof(gch_vector_list_sub_t));
    nvsl->count = numIndexes;
    nvsl->mainList = vsl->mainList;
    nvsl->realIndex = (int*)malloc(sizeof(int)*(nvsl->count));
    nvsl->vectors = (gch_vector_t**)malloc(sizeof(gch_vector_t*)*(nvsl->count));
    int *auxIndex = (int*)malloc(sizeof(int)*(nvsl->count));
    int *sortIndex = (int*)malloc(sizeof(int)*(nvsl->count));
    int i;
    for (i = 0; i < numIndexes; i++) { auxIndex[i] = vsl->realIndex[listIndexes[i]]; }
    gch_int_array_sort(auxIndex, numIndexes, sortIndex);
    for (i = 0; i < numIndexes; i++)
      { nvsl->realIndex[i] = auxIndex[sortIndex[i]];
        nvsl->vectors[i] = gch_affine_map_map(bt,vsl->vectors[listIndexes[sortIndex[i]]]);
      }
    free(auxIndex);
    free(sortIndex);
    return nvsl;
  }

void gch_vector_list_sub_free(gch_vector_list_sub_t *vsl)
  {
    int i;
    for (i = 0; i < vsl->count; i++) { gch_vector_free(vsl->vectors[i]); }
    free(vsl->vectors);
    free(vsl->realIndex);
    free(vsl);
  }

gch_vector_t* gch_vector_list_real_index(gch_vector_list_sub_t *vsl,int i)
  {  /* Binary search in array {realIndex}: */
    int left = 0;
    int right = vsl->count-1;
    int index = -1;
    int middle;
    do 
      { middle = (left+right)/2;
        if (vsl->realIndex[middle] == i) 
          { index = middle; }
        else if (vsl->realIndex[middle] > i) 
          { right = middle-1; }
        else 
          { left = middle+1; }
      } 
    while ((left <= right) && (index == -1));
    return vsl->vectors[index];
  }


