/* random bipartite graph generation. Based on simulation of a Poisson
random variable, which approximates the binomial random variable in a
real random graph.  This approximation is good if number of vertices
is large and prob. of choosing an edge is small, which is my case,
since I want sparse graphs.  Cf. Ross, A First Course in Probability,
p. 129.

Variation: this program generates graphs with a given number of expected edges,
instead of accepting the probability as input

Output file format specialized for bipartite graphs.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "misc.h"
#include "utility.h"
#include "bgraph.h"
#include "random.h"
#include <sys/time.h>

/* globals */
int Nvertices, Nedges;  /* Nvertices refers to # vertices in ONE partition */
float ene; 		/* expected number of edges for each vertex */
double p;		/* probability that an edge is chosen */
int seed;
#define MAXK  100
int  A[MAXK];	/* stores edges to be chosen by each vertex */
int B[MAX_CARD];

#define MAXRAND (double) ((1<<30)- 1 + (1<<30))

/* prototypes */

void ParseCommandLine(int argc, char **argv);
void GenerateGraph();

/* function definitions */

void main(int argc, char **argv)
{
  ParseCommandLine(argc, argv);
  GenerateGraph();
}

/*	*	*	*	*	*	*/

void ParseCommandLine(int argc, char **argv)
{
  struct timeval tp;

  if ((argc < 3) || (argc > 4))
    {
      printf("Program usage: rbg1 <nvertices/2> <ene> [seed]\n");
      exit(1);
    }

  Nvertices = atoi(argv[1]);
  ene = atof(argv[2]);
  /*  printf("\n%d %f\n",Nvertices,ene); exit(-1);*/
  if (argc == 4)
    seed = atoi(argv[3]);
  else
    {
      gettimeofday(&tp, 0);
      seed = tp.tv_sec + tp.tv_usec;
    }
}

void GenerateGraph()
{

  int v, i, t;
  double u, limit;
  int r, k;
  double lambda;

  /* the simulation of a Poisson random variable is done according to Ross,
     op. cit., p.398. k is the number of edges to be chosen, which will have
     a Poisson distribution.
     */

  InitRandom(seed);
  p = (double) ene/Nvertices;
  lambda = Nvertices * p; /* Poisson parameter */
  limit = 1/exp(lambda);
  /* initialize A */
  for (i = 0; i < MAXK; i++)
    {
      A[i] = i;
    }

  t = Nvertices*2; /* actual number of vertices is t */

  /* generate array B, with vertex numbers randomly distributed */
  for (i = 0; i < Nvertices; i++) B[i] = i;
  RandomPermutation(B, Nvertices);

  printf("c random bipartite graph; U = %d, p = %f, ene = %f, seed = %d\n",
	 Nvertices, p, ene, seed);
  printf("c U: size of one partition (partitions are of equal size)\n");
  printf("c p: probability of choosing an edge\n");
  printf("c ene: expected number of edges for each vertex\n");
 /* other descriptive lines are filled out at the end */

  Nedges = 0;
  for (v = 0; v < Nvertices; v++)
    {
      u = 1;
      k = 0;
      while(1)
	{
	  u *= random()/MAXRAND;
	  if ((u < limit) || k == (Nvertices-1) || k == MAXK) break;
	  k++;
	}
      RandomSubset(0, Nvertices-1, k, A);
      /* now the first k elements of A contain the edges to be chosen
       */
      for (i = 0; i < k; i++)
	{
	  printf("a %d %d\n", B[v]+1, A[i] + 1);
	}
      Nedges += k;
    }

  printf("p bipartite-matching %d %d %d\n", Nvertices, Nvertices, Nedges);
  exit(1);
}
