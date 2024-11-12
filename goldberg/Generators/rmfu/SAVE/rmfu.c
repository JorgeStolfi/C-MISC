/* 
   rmfu
   Generator of networks with unit capasities
   Modification of RMFGEN (Goldfarb and Grigoriadis)
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
int a=0, b=0, c=0, seed; /* input parameters */
int n, m, i, j, l, seed0, a2, an, bn, i_from, i_to;
int s, t;          /* source and sink */
int argN;

void error (int errorNo)

{
  switch ( errorNo ) {

  case 1: {
    fprintf ( stderr, "\nrmfu generates networks with unit capasities\n");
    fprintf ( stderr, "Usage: rmfu -aA -bB -cC  [-sSeed]\n");
    fprintf ( stderr, "  A - the size of layer (a x a)\n");
    fprintf ( stderr, "  B - the number of layers\n");
    fprintf ( stderr, "  C - the number of arcs between two layers\n");
    break;
  }
  case 2: {
    fprintf ( stderr, "'%c' must be positive\n", errc);
    break;
  }

  case 3: {
    fprintf ( stderr, "'c' cannot be more than 'a x a'\n" );
    break;
  }
  }

  exit(errorNo);
}

void printArc(int i, int j)

{
  printf("a %d %d 1\n", i, j);
}

int num ( int i, int j )
{
if (i < 0) i = a-1;
if (j < 0) j = a-1;
if (i >= a) i = 0;
if (j >= a) j = 0;
return i_to + 1 + i*a + j;
}

void main (int argc, char* argv[])
{
/* initial values */
   seed0 = 3287;

  if ((argc < 4) || (argc > 5)) error(1);

  for ( argN = 1; argN < argc; argN++ ) {
    strcpy ( args, argv[argN] );
    if ( args[0] != DASH ) error (1);
    switch ( args[1] )
      {
      case 'a': {
	a  =  atoi(&args[2]);
        if ( a < 1 ) { errc = 'a'; error(2); }
	break;
      }

      case 'b': {
	b  =  atoi(&args[2]);
        if ( b < 1 )  { errc = 'b'; error(2); }
	break;
      }

      case 'c': {
	c  =  atoi(&args[2]);
        if ( c < 1 )  { errc = 'c'; error(2); }
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
  if ( a == 0 || b == 0 || c == 0 ) error(1);
  a2 = a*a;
  if ( c > a2 )  error(3); 
  if ( seed == 0 ) seed = seed0;
  SetRandom(seed);

  /* # of source and sink */

  s = 1;
  t = 2;

  /* compute # nodes */

  n = 2 + a*a*b;

  /* compute # arcs */

  m = 4*a*a*b + c*(b+1);

  printf ("c rmfu a=%d  b=%d  c=%d\n", a, b, c);
  printf ("p max %d %d\n", n , m);
  printf ("n %d s\n", s);
  printf ("n %d t\n", t);

 /* interlayered arcs */
  /* source arcs */
  i_to = 2;
  for (j = 0; j < c; j++) {
    bn = Random ( (long)(i_to+1), (long)(i_to+a2) );
    printArc(s, bn);
  }

  for (l = 1; l < b; l++) {
    i_from = i_to;
    i_to += a2;

    for (j = 0; j < c; j++) {
      an = Random ( (long)(i_from+1), (long)(i_from+a2) );
      bn = Random ( (long)(i_to+1), (long)(i_to+a2) );
      printArc(an, bn);
    }
  }

  /* sink arcs */
  i_from = i_to;
  for (j = 0; j < c; j++) {
    an = Random ( (long)(i_from+1), (long)(i_from+a2) );
    printArc(an, t);
  }

 /* internal arcs */

  i_to   = 2 ;

  for (l = 1; l <= b; l++, i_to += a2) {

    for (i = 0; i < a; i++ )
    for (j = 0; j < a; j++ ) {
      an =  num(i,j);
      printArc (an, num(i-1,j));
      printArc (an, num(i+1,j));
      printArc (an, num(i,j-1));
      printArc (an, num(i,j+1));
    }


  }


  exit (0);
}

