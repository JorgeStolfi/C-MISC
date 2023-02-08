/* Last edited on 2014-07-22 12:56:24 by stolfilocal */

#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>

#include <gch_matrix.h>
#include <gch_vector.h>
#include <gch_vector_list.h>

#define DIMENSION_IN 4
#define DIMENSION_OUT 3

gch_matrix_t* gch_point_proj_matrix_make(int distortion);

int main(int argc, char **argv)
  {
    if (argc != 4)
      { fprintf(stderr,"usage: %s <in-file> <out-file> <distortion>\n",argv[0]);
        return -1;
      }
    char *fnameIn = argv[1];        /* Name of input file. */
    char *fnameOut = argv[2];       /* Name of output file. */
    int distortion = atoi(argv[3]); /* Distortion parameter. */

    int i,j;
    gch_matrix_t *projMat = gch_point_proj_matrix_make(distortion);
    gch_vector_list_t *vl = gch_vector_list_read(fnameIn);
    FILE *out = fopen(fnameOut,"w");
    fprintf(out,"%d %d\n", DIMENSION_OUT+1, vl->count);
    for (i = 0; i < vl->count; i++)
      { gch_vector_t *v = vl->vectors[i];
        gch_vector_t *w = gch_vector_new(DIMENSION_IN+1);
        mpz_set_si(w->coord[0],1);
        for (j = 0; j < DIMENSION_IN; j++)
          { mpz_set(w->coord[j+1],v->coord[j]); }
        gch_vector_t *u = gch_matrix_map_col(projMat,w);
        for (j = 0; j < DIMENSION_OUT+1; j++)
          { if (j != 0) fprintf(out," ");
            mpz_out_str (out, 10, u->coord[j]);
          }
        fprintf(out,"\n");
        gch_vector_free(w);
        gch_vector_free(u);
      }
    fclose(out);
    gch_vector_list_free(vl);
    gch_matrix_free(projMat);
    return 0;
  }

gch_matrix_t* gch_point_proj_matrix_make(int distortion)
  {
    int projMat_[4][5] =
      { { 1000,   0,    0,    0, distortion},
        {   0, 1000,    0,    0, 0},
        {   0,    0, 1000,    0, 0},
        {   0,    0,    0, 1000, 0}
      };
    gch_matrix_t* m = gch_matrix_new(DIMENSION_OUT+1, DIMENSION_IN+1);
    int i,j;
    for (i = 0; i < DIMENSION_OUT+1; i++)
      { for (j = 0; j < DIMENSION_IN+1; j++)
          { mpz_set_si(gch_matrix_elem(m,i,j),projMat_[i][j]); }
      }
    return m;
  }
