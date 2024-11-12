INTERFACE EdgeData;

IMPORT LR3;
FROM Integrator IMPORT Vector;

TYPE
  Data = RECORD
      pa, pb, va, vb,
      n, dn: LR3.T;
    END;
  T = REF ARRAY OF RECORD
      data: Data;
      ready: BOOLEAN;
    END;

PROCEDURE ComputeData(VAR data: Data; a, b, c, d: CARDINAL; pos, vel: Vector);
PROCEDURE GetData(data: T; l: CARDINAL; a, b, c, d: CARDINAL;
                  pos, vel: Vector): Data;
PROCEDURE GetReadyData(data: T; l: CARDINAL): Data;
PROCEDURE Discard(data: T);

END EdgeData.
