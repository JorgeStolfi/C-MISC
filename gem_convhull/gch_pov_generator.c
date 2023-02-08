/* Last edited on 2014-07-22 12:52:37 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <gmp.h>

#include <gch_util.h>
#include <gch_vector.h>
#include <gch_vector_list.h>

#define DIMENSION 3

void gch_vector_print_pov(FILE *fp, gch_vector_t *v, double scale, int homo);
  /* Prints to {fp} the vector {v} formatted as a POV-Ray point, with coordinates
    converted to {double} and multiplied by the given {scale}.
    If {homo} is true, assumes that each vector is four homogeneous coordiantes,
    with weight in element 0; otherwise assumes that each vector is three Cartesian
    coordinates.  */

int main(int argc, char **argv)
  { if (argc != 6)
      { fprintf(stderr,"usage: %s <pointFile> <faceFile> <povFile> <scale> <homoflag>\n", argv[0]);
        return -1;
      }
    char *fnamePt = argv[1]; /* Name of file with the point coordinates. */
    char *fnameFc = argv[2]; /* Name of file with the 2-faces. */
    char *fnamePov = argv[3]; /* Name of output POV-Ray file. */
    double scale = atof(argv[4]); /* Scale factor. */
    int homo = atoi(argv[5]); /* True iff input coordinates are homogeneous. */

    FILE *filePov = fopen(fnamePov,"w");
    fprintf(filePov, "#include \"face-gen.inc\"\n");
    fprintf(filePov, "\n");
    fprintf(filePov, " union{\n");

    /* Read the vertex list, output balls for them: */
    gch_vector_list_t *vl = gch_vector_list_read(fnamePt);
    int i,j;
    for (i = 0; i < vl->count; i++)
      { fprintf(filePov,"  object{ Make_Vertex(");
        gch_vector_print_pov(filePov, vl->vectors[i], scale, homo);
        fprintf(filePov,") }\n");
      }

    /* Read the 2-face list, output polygons for them: */
    int nvmax = 100; /* Max vertex count of a face. */
    int ncmax = nvmax * 20; /* Max chars in list of vertex indices. */
    int ixv[nvmax]; /* Indices of face vertices in vertex list. */
    char line[ncmax];
    FILE *fileFc = fopen(fnameFc, "r");
    while (!feof(fileFc))
      { if (fgets(line, 1000, fileFc))
          { int nv = gch_parse_int_vals(line, nvmax, ixv);
            fprintf(filePov,"  object{ Make_Face(%d, array[%d] {", nv, nv);
            for (j = 0; j < nv; j++)
              { gch_vector_print_pov(filePov, vl->vectors[ixv[j]], scale, homo);
                fprintf(filePov, " ");
              }
            fprintf(filePov,"} ) }\n");
          }
      }
    fclose(fileFc);

    gch_vector_list_free(vl);

    fprintf(filePov," }\n"); /* Closes the "union{". */
    fclose(filePov);
    return 0;
  }

void gch_vector_print_pov(FILE *fp, gch_vector_t *v, double scale, int homo)
  { int j;
    double wt = (homo ? mpz_get_d(v->coord[0]) : 1.0); /* Homogeneous weight. */
    fprintf(fp, " < ");
    for (j = 0; j < DIMENSION; j++)
      { if (j != 0) fprintf(fp, ", ");
        double coord;
        if (homo)
          { coord = mpz_get_d(v->coord[j+1])/wt*scale; }
        else
          { coord = mpz_get_d(v->coord[j])*scale; }
        fprintf(fp, "%f", coord);
      }
    fprintf(fp, " > ");
  }
