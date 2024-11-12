#include <stdint.h>
#include <stdlib.h>
#include <affirm.h>

#include <utable.h>

uint64_t utable_num_entries(utable_t *T) 
  { return T->np; }
  
utable_t *utable_new(uint64_t ne)
  { utable_t *T = notnull(malloc(sizeof(utable_t)), "no mem");
    T->e = notnull(malloc(ne*sizeof(uint64_t)), "no mem");
    uint64_t i;
    for (i = 0; i < ne; i++) { T->e[i] = utable_NULL; }
    T->ne = ne;
    T->np = 0;
    return T;
  }

uint64_t utable_get(utable_t *T, uint64_t id)
  {
    demand(id < T->ne, "invalid identifier");
    return T->e[id];
  }
  
void utable_set(utable_t *T, uint64_t id, uint64_t val)
  {
    demand(id < T->ne, "invalid identifier");
    demand(T->e[id] == utable_NULL, "identifier already has a value");
    demand(val != utable_NULL, "new value cannot be null");
    T->e[id] = val;
    T->np ++;
  }
