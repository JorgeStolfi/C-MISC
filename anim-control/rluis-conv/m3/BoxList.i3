INTERFACE BoxList;

IMPORT PairTable, SystemTopology, SystemDynamics;

TYPE
  Axis = PairTable.Axis;
  
  T = REF ARRAY OF RECORD
      min, max: ARRAY Axis OF LONGREAL;
    END;

PROCEDURE BuildVertexBoxes(bx: T; s: SystemDynamics.T; dt: LONGREAL);
PROCEDURE BuildEdgeBoxes(bx, vertexBoxes: T; s: SystemTopology.T);
PROCEDURE BuildFaceBoxes(bx, edgeBoxes: T; s: SystemTopology.T);

END BoxList.
