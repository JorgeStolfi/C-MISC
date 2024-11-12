MODULE EdgeData;

IMPORT LR3, LR3Extras;
FROM Integrator IMPORT Vector;

PROCEDURE ComputeData(VAR data: Data; a, b, c, d: CARDINAL; pos, vel: Vector) =
BEGIN
  data.pa := LR3.T{pos[a], pos[a+1], pos[a+2]};
  data.pb := LR3.T{pos[b], pos[b+1], pos[b+2]};
  data.va := LR3.T{vel[a], vel[a+1], vel[a+2]};
  data.vb := LR3.T{vel[b], vel[b+1], vel[b+2]};
  WITH
    pc  = LR3.T{pos[c], pos[c+1], pos[c+2]},
    pd  = LR3.T{pos[d], pos[d+1], pos[d+2]},
    vc  = LR3.T{vel[c], vel[c+1], vel[c+2]},
    vd  = LR3.T{vel[d], vel[d+1], vel[d+2]},
    pm  = LR3.Scale(0.5D0, LR3.Add(pc, pd)),
    vm  = LR3.Scale(0.5D0, LR3.Add(vc, vd)),
    ba  = LR3.Sub(data.pb, data.pa),
    ma  = LR3.Sub(pm, data.pa),
    dba = LR3.Sub(data.vb, data.va),
    dma = LR3.Sub(vm, data.va)
  DO
    data.n  := LR3Extras.Cross(ba, ma);
    data.dn := LR3.Add(LR3Extras.Cross(dba, ma), LR3Extras.Cross(ba, dma));
  END
END ComputeData;

PROCEDURE GetData(data: T; l: CARDINAL; a, b, c, d: CARDINAL; 
                  pos, vel: Vector): Data =
BEGIN
  IF TRUE OR NOT data[l].ready THEN
    ComputeData(data[l].data, a, b, c, d, pos, vel);
    data[l].ready := TRUE;
  END;
  RETURN data[l].data;
END GetData;

PROCEDURE GetReadyData(data: T; l: CARDINAL): Data =
BEGIN RETURN data[l].data END GetReadyData;

PROCEDURE Discard(data: T) =
BEGIN
  FOR i := 0 TO LAST(data^) DO
    data[i].ready := FALSE
  END
END Discard;

BEGIN END EdgeData.
