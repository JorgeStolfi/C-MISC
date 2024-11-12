/* 
   H8
   Generator of hard bipartite graphs
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
char args[30];
long int k, d, l, seed, aPerm, sn; /* input parameters */
long int m, seed0, i, _i, j, di, ki, _ki, i0, mi, a, sm, an;
long int s, t;          /* source and sink */
long int argN, YStart;
long int *idPerm;
long int  *iArc, *arcPerm;
long int flowForm;
long int *idPerm2,*loopPerm1,*loopPerm2;
long int XSize, YSize;


void error (int errorNo)

{
  switch ( errorNo ) {

  case 1: {
    fprintf ( stderr, "\nh8 generates bipartite graphs (X+Y, E)\n");
    fprintf ( stderr, "Usage: h8 k [-lL] [-dD] [-p[Seed]] [-n] [-sS]\n");
    fprintf ( stderr, "  |X|=|Y|=k*l\n");
    fprintf ( stderr, "  l - number of pieces (l=1 - default)\n");
    fprintf ( stderr, "  2*d - maximal degree of X-node (d=10 default)\n");
    fprintf ( stderr, "  -p - produces node random permutation\n");
    fprintf ( stderr, "  -n - don't permutate arcs\n");
    fprintf ( stderr, "  s>0 (s<0) - the number of extra-sources (sinks)\n");
    break;
  }
  case 2: {
    fprintf ( stderr, "k must be positive\n");
    break;
  }
  case 3: {
    fprintf ( stderr, "l must be positive\n");
    break;
  }
  case 4: {
    fprintf ( stderr, "d must be more than 1\n");
    break;
  }
  }

  exit(errorNo);
}

void printArc(long int i, long int j)

{
if (flowForm)
  printf("a %ld %ld 1\n",idPerm[i-1], idPerm[j-1]);
else
  printf("a %ld %ld\n",idPerm[i-1], idPerm2[j-1]);
}

void main (int argc, char* argv[])
{
/* defauld values */
  flowForm = 0;
  l = 1;
  d = 10;
  seed = 0;
  seed0 = 3287;
  aPerm = 1;
  sn = 0;
   
  if ((argc < 2) || (argc > 7)) error(1);

  /* first parameter -- x */
  strcpy (args, argv[1]);
  if (( k = atoi(argv[1])) < 1) error (2);
 
  for ( argN = 2; argN < argc; argN++ ) {
    strcpy ( args, argv[argN] );
    if ( args[0] != DASH ) error (1);
    switch ( args[1] )
      {
      case 'p': {
	seed  =  atoi(&args[2]);
        if ( seed <= 0 ) seed = seed0;
	break;
      }

      case 'l': {
	l  =  atoi(&args[2]);
        if ( l < 1 ) error(3);
	break;
      }

      case 'd': {
	d  =  atoi(&args[2]);
        if ( d <= 1 ) error(4);
	break;
      }

      case 's': {
	sn  =  atoi(&args[2]);
	break;
      }

      case 'n': {
	aPerm = 0;
	break;
      }

      case 'f': {
	flowForm = 1;
	break;
      }

      default:
	error(1);
      }
  }

  /* compute # of source and sink */

  if ( d > k ) d = k;

  if (flowForm) {
    s = 2*l*k + 1;
    t = s + 1;
  }
  
  /* compute # arcs */
  if (flowForm)
    m = 2*l*k;
  else
    m = 0;
	
  for (i = 1; i <= k; i++)
    m += (2*l-1) * MIN (i, d);
  
  if (!flowForm) {
    XSize = l*k;
    YSize = l*k;
    if (sn>=0) XSize += sn;
    else YSize += (-sn);
  }

  sm = ( sn >= 0 ) ? sn : -sn;
  m += sm * (d+1);

  iArc    = (long int *) calloc ( k*l+1, sizeof(long int) );
  arcPerm = (long int *) calloc ( k*l+1, sizeof(long int) );
  loopPerm1 = (long int *) calloc (k+1, sizeof(long int));
  loopPerm2 = (long int *) calloc (l+1, sizeof(long int));


  if (flowForm) {
    idPerm = (long int *) calloc(t + sm + 1, sizeof(long int));
    for (i = 0; i < t + sm; i++)
      idPerm[i] = i+1;
    if ( seed != 0 ) {
      SetRandom(seed);
      RandomPermute(idPerm, t + sm);
    }
  }
  else {
    idPerm = (long int *) calloc(XSize, sizeof(long int));
    idPerm2 = (long int *) calloc(YSize, sizeof(long int));
    for (i =0; i<XSize;i++)
      idPerm[i]=i+1;
    for (i=0; i<YSize;i++)
      idPerm2[i]=i+1;
    if (seed != 0) {
      SetRandom(seed);
      RandomPermute(idPerm,XSize);
      RandomPermute(idPerm2,YSize);
    }
  }

  printf("c h8 k=%ld, l=%ld, d=%ld ", k, l, d);
  if ( seed != 0 ) printf ("seed=%ld ", seed);
  if ( aPerm == 0 ) printf ("n ");
  printf ("\n");
  if (flowForm) {
    printf ("p max %ld %ld\n", t+sm , m);
    printf ("n %ld s\n", idPerm[s-1]);
    printf ("n %ld t\n", idPerm[t-1]);
  }
  else {
    printf("p bipartite-matching %ld %ld %ld\n",XSize,YSize,m);
  }


  if (flowForm) {
    /* source arcs */
    for (i = 1; i <= l*k; i++)
      {
	printArc(s, i);
      }
    /* sink arcs */
    for (i = l*k+1; i <= 2*l*k; i++)
      {
	printArc( i, t);
      }
  }
  if (flowForm) YStart = k*l;
  else YStart = 0;

  for (i=0;i<k;i++) loopPerm1[i]=i+1;
  RandomPermute(loopPerm1,k);
  for (i=0;i<l;i++) loopPerm2[i]=i;
  RandomPermute(loopPerm2,l);


  for (_i = 1; _i <= k; _i++) {
    i = loopPerm1[_i-1];
    // i = _i;
    di = MIN ( i, d );

    for (_ki = 0; _ki < l; _ki++) {
      ki = loopPerm2[_ki];
      // ki = _ki;
      i0 = ki*k+i;
      mi = 0;

      if ( ki != l-1 ) 
	for (j = 1; j <= di; j++) {
	  iArc[mi] = YStart+(ki+1)*k+i-di+j;
	  arcPerm [mi] = mi;
	  mi++;
	}

      for (j = 1; j <= di; j++) {
	iArc[mi] = YStart+(ki)*k+i-di+j; 
	arcPerm [mi] = mi;
	mi++;
      }

      if ( aPerm && mi > 1 )
	RandomPermute( arcPerm, mi );

      for ( a = 0; a < mi; a ++ ) {
	printArc ( i0, iArc[ arcPerm [a] ] );
      }
    }

  }

  for ( i = 1; i <= sm; i ++ ) {
    if ( sn > 0 ) {
      if (flowForm) printArc ( s, t+i ); 

      for (j = 1; j <= d; j++) {
	an =  Random ( (long)YStart+1,(long)(YStart+k*l) );
	if (flowForm) printArc ( t+i, an );
	else printArc ( k*l+i, an );
      }
    }
    else {
      if (flowForm) printArc ( t+i, t ); 
      for (j = 1; j <= d; j++)
	{
	  an = Random ( (long)1, (long)(k*l) );
	  if (flowForm) printArc ( an, t+i );
	  else printArc ( an, k*l+i);
	}
    }
  }

  exit (0);
}

