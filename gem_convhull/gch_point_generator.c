/* Last edited on 2024-12-21 10:50:26 by stolfi */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include <jsrandom.h>

int32_t main(int32_t argc, char **argv)
  {
    if (argc != 5)
      { fprintf(stderr,"usage: %s <filename> <numPoints> <dim> <maxCoord>", argv[0]);
        return -1;
      }

    char *fname = argv[1];
    int32_t numPoints = atoi(argv[2]);
    int32_t dim = atoi(argv[3]);
    int32_t maxCoord = atoi(argv[4]);

    /* Writes to file "{fname}" a list of {numPoints} points
      of {R^{dim}} with random integer coordinates uniformly
      chosen between {-maxCoord} and {+maxCoord}.
      The first line has the two numbers "{dim} {numPoints}". */

    FILE *fp = fopen(fname, "w");

    fprintf(fp, "%d %d\n", dim, numPoints);

    for (int32_t i = 0; i < numPoints; i++)
      { for (int32_t j = 0; j < dim; j++)
          { if (j != 0) fprintf(fp, " ");
            int32_t val = int32_abrandom(-maxCoord, +maxCoord);
            fprintf(fp, "%d", val);
          }
        fprintf(fp, "\n");
      }
    fclose(fp);
    return 0;
  }
