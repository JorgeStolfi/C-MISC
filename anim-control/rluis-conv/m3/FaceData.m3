MODULE FaceData;

IMPORT LR3, LR3Extras;
FROM Integrator IMPORT Vector;

PROCEDURE ComputeData(VAR data: Data; a, b, c: CARDINAL; pos, vel: Vector) =
BEGIN
  data.pa := LR3.T{pos[a], pos[a+1], pos[a+2]};
  data.pb := LR3.T{pos[b], pos[b+1], pos[b+2]};
  data.pc := LR3.T{pos[c], pos[c+1], pos[c+2]};
  data.va := LR3.T{vel[a], vel[a+1], vel[a+2]};
  data.vb := LR3.T{vel[b], vel[b+1], vel[b+2]};
  data.vc := LR3.T{vel[c], vel[c+1], vel[c+2]};
  WITH
    u  = LR3.Sub(data.pb, data.pa),
    v  = LR3.Sub(data.pc, data.pa),
    w  = LR3.Sub(data.pc, data.pb),
    du = LR3.Sub(data.vb, data.va),
    dv = LR3.Sub(data.vc, data.va),
    dw = LR3.Sub(data.vc, data.vb),
    n  = LR3Extras.Cross(u, v),
    dn = LR3.Add(LR3Extras.Cross(du, v), LR3Extras.Cross(u, dv))
  DO
    data.n   := n;
    data.n1  := LR3Extras.Cross(n, w);
    data.n2  := LR3Extras.Cross(v, n);
    data.n3  := LR3Extras.Cross(n, u);
    data.dn  := dn;
    data.dn1 := LR3.Add(LR3Extras.Cross(dn, w), LR3Extras.Cross(n, dw));
    data.dn2 := LR3.Add(LR3Extras.Cross(dv, n), LR3Extras.Cross(v, dn));
    data.dn3 := LR3.Add(LR3Extras.Cross(dn, u), LR3Extras.Cross(n, du));
  END
END ComputeData;

PROCEDURE GetData(data: T; f: CARDINAL; a, b, c: CARDINAL; 
                  pos, vel: Vector): Data =
BEGIN
  IF NOT data[f].ready THEN
    ComputeData(data[f].data, a, b, c, pos, vel);
    data[f].ready := TRUE;
  END;
  RETURN data[f].data;
END GetData;

PROCEDURE GetReadyData(data: T; f: CARDINAL): Data =
BEGIN RETURN data[f].data END GetReadyData;

PROCEDURE Discard(data: T) =
BEGIN
  FOR i := 0 TO LAST(data^) DO
    data[i].ready := FALSE
  END
END Discard;

BEGIN END FaceData.
