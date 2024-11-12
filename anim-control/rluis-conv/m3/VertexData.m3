MODULE VertexData;

IMPORT LR3;
FROM Integrator IMPORT Vector;

PROCEDURE ComputeData(VAR data: Data; a: CARDINAL; pos, vel: Vector) =
BEGIN
  data.p := LR3.T{pos[a], pos[a+1], pos[a+2]};
  data.v := LR3.T{vel[a], vel[a+1], vel[a+2]};
END ComputeData;

PROCEDURE GetData(data: T; p: CARDINAL; a: CARDINAL; 
                  pos, vel: Vector): Data =
BEGIN
  IF NOT data[p].ready THEN
    ComputeData(data[p].data, a, pos, vel);
    data[p].ready := TRUE;
  END;
  RETURN data[p].data;
END GetData;

PROCEDURE GetReadyData(data: T; p: CARDINAL): Data =
BEGIN RETURN data[p].data END GetReadyData;

PROCEDURE Discard(data: T) =
BEGIN
  FOR i := 0 TO LAST(data^) DO
    data[i].ready := FALSE
  END
END Discard;

BEGIN END VertexData.
