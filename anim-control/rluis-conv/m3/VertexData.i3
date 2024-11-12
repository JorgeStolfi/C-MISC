INTERFACE VertexData;

IMPORT LR3;
FROM Integrator IMPORT Vector;

TYPE
  Data = RECORD
      p, v: LR3.T;
    END;
  T = REF ARRAY OF RECORD
      data: Data;
      ready: BOOLEAN;
    END;

PROCEDURE ComputeData(VAR data: Data; a: CARDINAL; pos, vel: Vector);
PROCEDURE GetData(data: T; p: CARDINAL; a: CARDINAL; 
                  pos, vel: Vector): Data;
PROCEDURE GetReadyData(data: T; p: CARDINAL): Data;
PROCEDURE Discard(data: T);

END VertexData.
