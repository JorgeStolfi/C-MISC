/* Last edited on 2024-12-21 11:26:20 by stolfi */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <jsrandom.h>

#define MAX_DIM 10

int32_t main(int32_t argc, char **argv)
  {
    if (argc != 5)
      { fprintf(stderr,"usage: %s <filename>  <numPoints>  <dimensions>  <radius> ", argv[0]);
        return -1;
      }

    char *fname = argv[1];
    int32_t numPoints = atoi(argv[2]);
    int32_t dim = atoi(argv[3]);
    int32_t radius = atoi(argv[4]);

    /* Writes to file "{fname}" a list of {numPoints} random points of {R^{dim}}
      with integer coordinates, located approximately at distance {radius}
      from the origin. The first line has the two numbers "{dim} {numPoints}". */

    double v[MAX_DIM];
    double length;

    FILE *fp = fopen(fname, "w");

    fprintf(fp, "%d %d\n", dim, numPoints);

    int32_t nt = 0; /* Number of points tried so far. */
    int32_t np = 0; /* Number of points written so far. */
    while (np < numPoints)
      { /* Generate a random point {v} in the cube {[-1 _ +1]^dim}: */
        length = 0.0;
        for (int32_t j = 0; j < dim; j++)
          { v[j] = abrandom(-1.0,+1.0);
            length += v[j]*v[j];
          }
        length = sqrt(length);

        /* If inside the unit ball, project to the surface: */
        if (length <= 1.0)
          { for (int32_t j = 0; j < dim; j++)
              { if (j != 0) fprintf(fp, " ");
                int32_t val = (int32_t)((double)radius*v[j]/length);
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
