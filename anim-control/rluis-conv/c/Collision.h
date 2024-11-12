
#ifndef Collision_H
#define Collision_H



#include <VertexData.h>
#include <EdgeData.h>
#include <FaceData.h>

typedef
  VertexToFaceDeviation == struct ?? {
      double H, G1, G2, G3, dH, dG1, dG2, dG3;
    }
  EdgeToEdgeDeviation == struct ?? {
      double H, G1, G2, G3, G4, dH, dG1, dG2, dG3, dG4;
    }
    
PROCEDURE ComputeVFDeviation(READONLY VertexData.Data p,
                             FaceData.Data *f);: VertexToFaceDeviation;

PROCEDURE ComputeEEDeviation(READONLY e1,
                             EdgeData.Data e2): EdgeToEdgeDeviation;

PROCEDURE DetectVFCollision(READONLY pa, VertexData.Data pb,
                            FaceData.Data *fa, fb;
                            double h3): double;

PROCEDURE DetectEECollision(READONLY e1a, e1b, e2a, EdgeData.Data e2b,
                            double h3;
                            bool *?sign): double;

PROCEDURE DetectVFContactBreak(READONLY pa, VertexData.Data pb,
                               FaceData.Data *fa, fb;
                               double h3): double;

PROCEDURE DetectEEContactBreak(READONLY e1a, e1b, e2a, EdgeData.Data e2b,
                               double h3; bool sign): double;
;
} /* Collision */.

#endif
