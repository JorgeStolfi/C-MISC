/* bipartite graph generator.
   Neighborhood by group. Each vertex u in group i has neighbors in groups
   i-1, i, i+1 in V. Wraps around.

   Output file format specialized for bipartite graphs.

   U vertices randomly renumerated in the output
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <misc.h>
#include <bgraph.h>
#include <random.h>
#include <sys/time.h>

#define MAXK 100

/* globals */
int Nvertices, Nedges;  /* Nvertices refers to # vertices in ONE partition */
int N;                  /* number of vertices in a group */
int Ngroups;
double ene; 		/* expected number of edges for each vertex */
double p;		/* probability that an edge is chosen */
int seed;
int  A[MAXK];	/* stores edges to be chosen by each vertex */
int B[MAX_CARD];
double Limit;

#define MAXRAND ((1<<30)-1+(1<<30))

/* prototypes */

void ParseCommandLine(int argc, char **argv);
void GenerateGraph();
void GenerateEdges(int w, int low, int high);
void PrintEdges(int w, int k);
int GetK();

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

  if ((argc < 4) || (argc > 5))
    {
      printf("Program usage: rbg <N> <G> <ene> [seed]\n");
      exit(1);
    }

  N = atoi(argv[1]);
  Ngroups = atoi(argv[2]);
  ene = atof(argv[3]);
  if (argc == 5)
    seed = atoi(argv[4]);
  else
    {
      gettimeofday(&tp, 0);
      seed = tp.tv_sec + tp.tv_usec;
    }
}

void GenerateGraph()
{

  int l, i, u, v;
  double p, lambda;

  Nvertices = N*Ngroups;

  /* simulation of a Poisson random variable */
  InitRandom(seed);
  p = (double) ene/Nvertices;
  lambda = Nvertices * p; /* Poisson parameter */
  Limit = 1/exp(lambda);

  /* initialize A */
  for (i = 0; i < MAXK; i++)
    {
      A[i] = i;
    }

  /* generate array B, with vertex numbers randomly distributed */
  for (i = 0; i < Nvertices; i++) B[i] = i;
  RandomPermutation(B, Nvertices);

  printf("c special random bipartite graph; N=%d, Ngroups=%d, ene=%.1f, seed=%d\n",
	 N, Ngroups, ene, seed);

  printf("c ene: expected number of edges for each vertex\n");
 /* other descriptive lines are filled out at the end */

  Nedges = 0;

  for (i = 0; i < Ngroups; i++)
    {
      for (u = i*N; u < (i+1)*N; u++)
	{
	  GenerateEdges(u, (i-1)*N, (i+2)*N - 1);
	}
    }
  printf("p bipartite-matching %d %d %d\n", Nvertices, Nvertices, Nedges);

}

void GenerateEdges(int w, int low, int high)
{ 
  int k;

  k = GetK();
  RandomSubset(low, high, k, A);
  PrintEdges(w, k);
}

int GetK()
{
  int k;
  double u;

  u = 1;
  k = 0;
  while(1)
    {
      u *= ((double) random()/ MAXRAND);
      if ((u < Limit) || k == (N-1) || k == MAXK) break;
      k++;
    }
  return k;
}

void PrintEdges(int w, int k)
{

  int i;

  for (i = 0; i < k; i++)
    {
      if (A[i] < 0) A[i] += Nvertices;
      else if (A[i] >= Nvertices) A[i] -= Nvertices;
      printf("a %d %d\n", B[w]+1, A[i] + 1);
    }
  Nedges += k;
}



