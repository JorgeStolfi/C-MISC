
#ifndef BoxList_H
#define BoxList_H



#include <PairTable.h>
#include <SystemTopology.h>
#include <SystemDynamics.h>

typedef
  Axis == PairTable.Axis;
  
  T == struct ?? {_vec
      min, max: ARRAY Axis OF double;
    }

void BuildVertexBoxes(T bx, SystemDynamics.T s, double dt);
void BuildEdgeBoxes(bx, T vertexBoxes, SystemTopology.T s);
void BuildFaceBoxes(bx, T edgeBoxes, SystemTopology.T s);

#endif
