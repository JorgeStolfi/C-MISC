INTERFACE PairTable;

IMPORT SystemTopology;

TYPE
  Axis = {X, Y, Z};
  Row = REF ARRAY OF RECORD
          v: BITS 1 FOR BOOLEAN;
          flag: ARRAY Axis OF BITS 1 FOR BOOLEAN;
        END;
  T = REF ARRAY OF RECORD
          on: BOOLEAN;
          row: Row;
        END;

PROCEDURE New(l, c: CARDINAL): T;
PROCEDURE On(t: T; i: CARDINAL): BOOLEAN;
PROCEDURE SetOn(t: T; i: CARDINAL; on: BOOLEAN);
PROCEDURE Valid(t: T; i, j: CARDINAL): BOOLEAN;
PROCEDURE SetValid(t: T; i, j: CARDINAL; v: BOOLEAN);
PROCEDURE Overlap(t: T; i, j: CARDINAL): BOOLEAN;
PROCEDURE BuildVertexFacePairs(t: T; s: SystemTopology.T): CARDINAL;
PROCEDURE BuildEdgeEdgePairs(t: T; s: SystemTopology.T): CARDINAL;

END PairTable.
