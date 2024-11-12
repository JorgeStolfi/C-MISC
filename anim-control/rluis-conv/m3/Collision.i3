INTERFACE Collision;

IMPORT VertexData, EdgeData, FaceData;

TYPE
  VertexToFaceDeviation = RECORD
      H, G1, G2, G3, dH, dG1, dG2, dG3: LONGREAL;
    END;
  EdgeToEdgeDeviation = RECORD
      H, G1, G2, G3, G4, dH, dG1, dG2, dG3, dG4: LONGREAL;
    END;
    
PROCEDURE ComputeVFDeviation(READONLY p: VertexData.Data;
                             READONLY f: FaceData.Data): VertexToFaceDeviation;

PROCEDURE ComputeEEDeviation(READONLY e1,
                             e2: EdgeData.Data): EdgeToEdgeDeviation;

PROCEDURE DetectVFCollision(READONLY pa, pb: VertexData.Data;
                            READONLY fa, fb: FaceData.Data;
                            h3: LONGREAL): LONGREAL;

PROCEDURE DetectEECollision(READONLY e1a, e1b, e2a, e2b: EdgeData.Data;
                            h3: LONGREAL;
                            VAR sign: BOOLEAN): LONGREAL;

PROCEDURE DetectVFContactBreak(READONLY pa, pb: VertexData.Data;
                               READONLY fa, fb: FaceData.Data;
                               h3: LONGREAL): LONGREAL;

PROCEDURE DetectEEContactBreak(READONLY e1a, e1b, e2a, e2b: EdgeData.Data;
                               h3: LONGREAL; sign: BOOLEAN): LONGREAL;

END Collision.
