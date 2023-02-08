/* Basic types and procedures. */
/* Last edited on 2008-07-14 22:20:28 by stolfi */

#ifndef basic_H
#define basic_H

#define _GNU_SOURCE
/* To get {asprintf}. */

#include <affirm.h>
#include <bool.h>
#include <sign.h>
#include <vec.h>
#include <r3.h>
#include <r4.h>
#include <interval.h>
#include <quad.h>
#include <frgb.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <values.h>

#define C_SQRT3 (1.73205080756887729352)

#ifndef INF
#define INF INFINITY
  /* Plus infinity as a {double} value. */
#endif

/* DEBUGGING OUTPUT */

#define mumble(...) \
  do { if (verbose) { fprintf(stderr, __VA_ARGS__); } } while(0)

#endif
