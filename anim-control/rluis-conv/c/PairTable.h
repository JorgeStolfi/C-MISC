
#ifndef PairTable_H
#define PairTable_H



#include <SystemTopology.h>

typedef
  Axis == {X, Y, Z}
  Row == struct ?? {_vec
          v: BITS 1 FOR bool;
          flag: ARRAY Axis OF BITS 1 FOR bool;
        }
  T == struct ?? {_vec
          bool on;
          Row row;
        }

T New(l, nat c);
bool On(T t, nat i);
void SetOn(T t, nat i, bool on);
bool Valid(T t, i, nat j);
void SetValid(T t, i, nat j, bool v);
bool Overlap(T t, i, nat j);
nat BuildVertexFacePairs(T t, SystemTopology.T s);
nat BuildEdgeEdgePairs(T t, SystemTopology.T s);
;
} /* PairTable */.

#endif
