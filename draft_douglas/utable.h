#ifndef utable_H
#define utable_H
/* Last edited on 2014-05-15 15:07:53 by stolfilocal */

#include <stdint.h>

#define utable_NULL (~(uint64_t)0)

typedef struct utable_t
  { uint64_t ne;  /* Number of elements allocated for {*e}. */
    uint64_t np;  /* Number of non-null elements in {e[0..ne-1]}. */
    uint64_t *e;  /* Each element is either {utable_NULL} or some other integer. */
  } utable_t;
  /* A table that associates /values/ (64-bit unsiged integers) to a set of
    /identifiers/ (also 64-bit unsiged integers).  */

uint64_t utable_num_entries(utable_t *T);
  /* Returns the number of values stored in {T}. */
  
utable_t *utable_new(uint64_t ne);
  /* Creates a table {T} that can contain up to {ne} 
    non-NULL values, identified by integers in {0..ne-1}.
    Initially there are no values stored. */
    
uint64_t utable_get(utable_t *T, uint64_t id);
  /* Returns the value associated to identifier {id} in table {T}.
    Returns {utable_NULL} if there is no value associated to {id}. */
  
void utable_set(utable_t *T, uint64_t id, uint64_t val);
  /* Associates {val} to {id} in table {T}.  There
    must be no value currently associated to {id}. */   

#endif
