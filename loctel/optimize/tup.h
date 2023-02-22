/* Public telephone structures and tools. */
/* Last edited on 2023-02-21 11:40:42 by stolfi */

#ifndef tup_H
#define tup_H

#define _GNU_SOURCE
#include <stdio.h>
#include <values.h>

#include <affirm.h>
#include <stmap.h>

typedef struct phone_t
  { int32_t vertex;     /* Vertex where phone is located. */
    double usage;   /* Predicted usage. */
  } phone_t;
  /* Data for a public telephone. */
    
/* Vectors of {phone_t}: */

vec_typedef(phone_vec_t,phone_vec,phone_t);

void st_write_phones(char *name, Map *m, phone_vec_t *ph);
  /* Writes a list of phones to a file called {name} (or to {stdout} if 
    {name = "-"}).  The phone coordinates are obtained from the map {m}. */

/* Distance matrices: */

typedef struct dpair_t
  { int32_t ui;  /* Map vertex. */
    int32_t pi;  /* Map vertex. */
    double dist; /* Distance on map from {pi} to {ui}. */
  } dpair_t;
  /* A record of walking distance between vertices {ui} and {pi}. */

vec_typedef(dpair_vec_t,dpair_vec,dpair_t);

dpair_vec_t st_vertex_phone_distances(Map *m, double maxDist, phone_vec_t *ph);
  /* Computes a list {H[0..nH-1]} of all user-phone pairs for 
    the street map {m} and phones {ph[0..nph-1]}. By definition,
    the users of a phone are all vertices of {m} which are {maxDist}-reachable
    from the phone's vertex. Allocates the array {nH} and sets {*nH}. */
  
void st_recompute_phone_usage
  ( Map *m, 
    double_vec_t *dem, 
    double K,
    double partDist, 
    phone_vec_t *ph,
    double_vec_t *lost 
  );
  /* Recomputes the estimated usage of each phone {ph[pj].usage},
    given the potential demand {dem[ui]} at each vertex {ui}. Assumes
    that an isolated phone can attract half of the potential demand
    from any vertex whose path-cost is {partDist}. Also stores
    in {lost[ui]} the lost demand at vertex {ui}. */
  
#endif
