/* SPTimeSpaceFuncMap.h -- pointwise operator for sphere-time functions. */
/* Last edited on 2005-08-14 17:30:17 by stolfi */

#ifndef SPTimeSpaceFuncMap_H
#define SPTimeSpaceFuncMap_H

#include <SPBasic.h>
#include <bool.h>

/* 
  An {SPTimeSpaceFuncMap} is an object that describes a pointwise 
  operator {F(u,p,t)} that modifies the value {u} of a sphere-time
  function, possibly depending on the point {p} and time {t}.
  
  This object is used mainly as the right-hand-side of the general
  differential equation {D(f)(p,t) == F(f(p,t),p,t)} where {D} is a
  linear sphere-time differential operator. It is also compatible with
  the {TimeSpaceFuncMap} argument of {SPTimeSpaceFunction_Dot}.
*/

typedef double TimeSpaceFuncMapMap(double u, S2Point *p, double t); 
  /* A procedure that implements a pointwise sphere-time operator. */

typedef struct SPTimeSpaceFuncMap   /* A pointwise sphere-time operator. */
  { TimeSpaceFuncMapMap *map;   /* The map proper. */
    bool_t zeroPres;    /* TRUE implies zero-preserving, i.e. {F(0,p,t) = 0}. */
    char *descr;        /* The formula for {F(u, (x,y,z), t)}. */
  } SPTimeSpaceFuncMap;
  /* A pointwise operator that sphere-time suitable  */

typedef SPTimeSpaceFuncMap TimeSpaceFuncMap; /* A shorter name. */

/* The identity operator: */
#define NoTSFMap (TimeSpaceFuncMap){ NULL, TRUE, "u" }  

SPTimeSpaceFuncMap SPTimeSpaceFuncMap_FromName(char *name);
  /* Returns a {TimeSpaceFuncMap} and its data given its name.
    Bombs out if the name is unrecognized. */

#endif

