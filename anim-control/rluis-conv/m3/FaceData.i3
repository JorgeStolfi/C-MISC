INTERFACE FaceData;

IMPORT LR3;
FROM Integrator IMPORT Vector;

TYPE
  Data = RECORD
      pa, pb, pc, va, vb, vc,
      n, n1, n2, n3, dn, dn1, dn2, dn3: LR3.T;
    END;
  T = REF ARRAY OF RECORD
      data: Data;
      ready: BOOLEAN;
    END;

PROCEDURE ComputeData(VAR data: Data; a, b, c: CARDINAL; pos, vel: Vector);
PROCEDURE GetData(data: T; f: CARDINAL;
                  a, b, c: CARDINAL; pos, vel: Vector): Data;
PROCEDURE GetReadyData(data: T; f: CARDINAL): Data;
PROCEDURE Discard(data: T);

END FaceData.
