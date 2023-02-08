/* Last edited on 2014-07-22 12:55:44 by stolfilocal */

#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv)
  {
    if (argc != 5)
      { fprintf(stderr,"usage: %s <filename> <numPoints> <dim> <maxCoord>", argv[0]);
        return -1;
      }

    char *fname = argv[1];
    int numPoints = atoi(argv[2]);
    int dim = atoi(argv[3]);
    int maxCoord = atoi(argv[4]);

    /* Writes to file "{fname}" a list of {numPoints} points
      of {R^{dim}} with random integer coordinates uniformly
      chosen between {-maxCoord} and {+maxCoord}.
      The first line has the two numbers "{dim} {numPoints}". */

    FILE *fp = fopen(fname, "w");

    fprintf(fp, "%d %d\n", dim, numPoints);

    int i, j;
    for (i = 0; i < numPoints; i++)
      { for (j = 0; j < dim; j++)
          { if (j != 0) fprintf(fp, " ");
            int val = (int)((((double)rand()*2.0)/(double)RAND_MAX - 1.0)*maxCoord);
            fprintf(fp, "%d", val);
          }
        fprintf(fp, "\n");
      }
    fclose(fp);
    return 0;
  }
