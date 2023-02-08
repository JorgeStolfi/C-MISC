/* Last edited on 2014-07-22 01:31:57 by stolfilocal */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_DIM 10

int main(int argc, char **argv)
  {
    if (argc != 5)
      { fprintf(stderr,"usage: %s <filename>  <numPoints>  <dimensions>  <radius> ", argv[0]);
        return -1;
      }

    char *fname = argv[1];
    int numPoints = atoi(argv[2]);
    int dim = atoi(argv[3]);
    int radius = atoi(argv[4]);

    /* Writes to file "{fname}" a list of {numPoints} random points of {R^{dim}}
      with integer coordinates, located approximately at distance {radius}
      from the origin. The first line has the two numbers "{dim} {numPoints}". */

    double v[MAX_DIM];
    double length;

    FILE *fp = fopen(fname, "w");

    fprintf(fp, "%d %d\n", dim, numPoints);

    int nt = 0; /* Number of points tried so far. */
    int np = 0; /* Number of points written so far. */
    while (np < numPoints)
      { /* Generate a random point {v} in the cube {[-1 _ +1]^dim}: */
        length = 0.0;
        int j;
        for (j = 0; j < dim; j++)
          { v[j] = ((double)rand()*2.0)/(double)RAND_MAX - 1.0;
            length += v[j]*v[j];
          }
        length = sqrt(length);

        /* If inside the unit ball, project to the surface: */
        if (length <= 1.0)
          { for (j = 0; j < dim; j++)
              { if (j != 0) fprintf(fp, " ");
                int val = (int)((double)radius*v[j]/length);
                fprintf(fp, "%d", val);
              }
            fprintf(fp,"\n");
            np++;
          }
        nt++;
      }

    fclose(fp);
    fprintf(stderr, "generated %d points in %d attempts.\n", np, nt);
    return 0;
  }
