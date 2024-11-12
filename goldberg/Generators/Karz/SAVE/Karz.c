/* 
   DK
   Generator of hard graphs with unit capasity
   Follows Diniz-Karzanov hard bipartite example
   Written by Andrew V. Goldberg & Boris Cherkassky,
   NEC Research Inst.,  avg@research.nj.nec.com, cher@research.nj.nec.com.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <values.h>

#include "random.h" 

#define DASH '-'
#define MIN(x,y) ((x < y) ? x:y)

/* global variables */
char args[30], errc;
int k=0, d=0, l=0, a=0, f=0, seed; /* input parameters */
int n, m, i, j, seed0, an, i_from, i_to, p, v;
int s, t;          /* source and sink */
int argN;

void error (int errorNo)

{
  switch ( errorNo ) {

  case 1: {
    fprintf ( stderr, "\ndk generates hard flow graph with unit capasities\n");
    fprintf ( stderr, "Usage: dk -kK -lL -aA [-dD] [-fF] [-sSeed]\n");
    fprintf ( stderr, "  k - the number of paths\n");
    fprintf ( stderr, "  l*(i-1)+1 - the length of i-th path\n");
    fprintf ( stderr, "  a - the size of 'wide' piece\n");
    fprintf ( stderr,
              "  d the degree of a node in the wide piece ('a' default)\n");
    fprintf ( stderr,
              "  f the degree of a node in the narrow piece ('a' default)\n");
    break;
  }
  case 2: {
    fprintf ( stderr, "'%c' must be positive\n", errc);
    break;
  }

  case 3: {
    fprintf ( stderr, "'%c' cannot be more than 'a'\n", errc);
    break;
  }
  }

  exit(errorNo);
}

void printArc(int i, int j)

{
  printf("a %d %d 1\n", i, j);
}

void main (int argc, char* argv[])
{
/* initial values */
   seed0 = 3287;

  if ((argc < 4) || (argc > 7)) error(1);

  for ( argN = 1; argN < argc; argN++ ) {
    strcpy ( args, argv[argN] );
    if ( args[0] != DASH ) error (1);
    switch ( args[1] )
      {
      case 'k': {
	k  =  atoi(&args[2]);
        if ( k < 1 ) { errc = 'k'; error(2); }
	break;
      }

      case 'l': {
	l  =  atoi(&args[2]);
        if ( l < 1 )  { errc = 'l'; error(2); }
	break;
      }

      case 'a': {
	a  =  atoi(&args[2]);
        if ( a < 1 )  { errc = 'a'; error(2); }
	break;
      }

      case 'd': {
	d  =  atoi(&args[2]);
        if ( d < 1 )  { errc = 'd'; error(2); }
	break;
      }

      case 'f': {
	f  =  atoi(&args[2]);
        if ( f < 1 )  { errc = 'f'; error(2); }
	break;
      }

      case 's': {
	seed  =  atoi(&args[2]);
        if ( seed <= 0 ) seed = seed0;
	break;
      }

      default:
	error(1);
      }
  }
  /* checking parameters */
  if ( k == 0 || l == 0 || a == 0 ) error(1);
  if ( d == 0 ) d = a;
  if ( f == 0 ) f = a;
  if ( k > a )  { errc = 'k'; error(3); } 
  if ( d > a )  { errc = 'd'; error(3); } 
  if ( f > a )  { errc = 'f'; error(3); } 
  if ( seed == 0 ) seed = seed0;
  SetRandom(seed);

  /* # of source and sink */

  s = 1;
  t = 2;

  /* compute # nodes */

  n = 2 + k + 2*a + (k*(l*(k-1)+2))/2;

  /* compute # arcs */

  m = k + k*a + k*f + a*d + (k*(l*(k-1)+2))/2;

  printf ("c dk k=%d  l=%d  a=%d  d=%d\n", k, l, a, d);
  printf ("p max %d %d\n", n , m);
  printf ("n %d s\n", s);
  printf ("n %d t\n", t);

  /* source arcs */
  i_to = 2;
  for (i = 1; i <= k; i++) {
    printArc(s, i_to + i);
  }

  /* 1st level arcs */
  i_from = 2;
  i_to   = 2 + k;


  if ( f == a )
  for (i = 1; i <= k; i++)
  for (j = 1; j <= a; j++) {
    printArc( i_from + i, i_to + j );
  }
  else 
  for (i = 1; i <= k; i++)
  for (j = 1; j <= f; j++) {
    an = Random ( (long)(i_to+1), (long)(i_to+a) );
    printArc( i_from + i, an );
  }

  /* 2nd level arcs */
  i_from = 2 + k;
  i_to   = 2 + k + a;

  if ( d == a )
  for (i = 1; i <= a; i++)
  for (j = 1; j <= a; j++) {
    printArc( i_from + i, i_to + j );
  }
  else
  for (i = 1; i <= a; i++)
  for (j = 1; j <= d; j++) {
    an = Random ( (long)(i_to+1), (long)(i_to+a) );
    printArc( i_from + i, an );
  }
    
  /* paths */

  i_from = 2 + k + a;
  v = i_from + a + 1;

  for (p = 0; p < k; p++) {

    for (i = 1; i <= a; i ++) {
      printArc ( i_from + i, v );
    }

    for (j = 1; j <= p*l; j++) {
      printArc ( v, v+1 );
      v ++;
    }
    
    printArc ( v, t );
    v ++;
  }

  exit (0);
}

