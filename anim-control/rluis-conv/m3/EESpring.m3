MODULE EESpring;

IMPORT LR3;
FROM Force IMPORT Vector;

REVEAL
  T = Public BRANDED OBJECT
      a, b, c, d: CARDINAL;
      a1, a2, a3, a4, k: LONGREAL;
    OVERRIDES
      init    := Init;
      compute := Compute;
    END;

PROCEDURE Init(cf: T; a, b, c, d: CARDINAL; a1, a2, a3, a4, k: LONGREAL) =
BEGIN
  cf.a := a; cf.b := b; cf.c := c; cf.d := d;
  cf.a1 := a1; cf.a2 := a2; cf.a3 := a3; cf.a4 := a4;
  cf.k := k;
END Init;

PROCEDURE Compute(cf: T; pos, vel: Vector; f: Vector; <*UNUSED*> t: LONGREAL) =
VAR 
  a := cf.a; b := cf.b; c := cf.c; d := cf.d;
  a1 := cf.a1; a2 := cf.a2; a3 := cf.a3; a4 := cf.a4;
BEGIN
  WITH
    pa = LR3.T{pos[a], pos[a+1], pos[a+2]},
    pb = LR3.T{pos[b], pos[b+1], pos[b+2]},
    pc = LR3.T{pos[c], pos[c+1], pos[c+2]},
    pd = LR3.T{pos[d], pos[d+1], pos[d+2]},
    q  = LR3.Sub(LR3.Mix(a1, pa, a2, pb), LR3.Mix(a3, pc, a4, pd)),
    fx = cf.k*q[0],
    fy = cf.k*q[1],
    fz = cf.k*q[2]
  DO
    f[a+0] := f[a+0] - a1*fx;
    f[a+1] := f[a+1] - a1*fy;
    f[a+2] := f[a+2] - a1*fz;

    f[b+0] := f[b+0] - a2*fx;
    f[b+1] := f[b+1] - a2*fy;
    f[b+2] := f[b+2] - a2*fz;

    f[c+0] := f[c+0] + a3*fx;
    f[c+1] := f[c+1] + a3*fy;
    f[c+2] := f[c+2] + a3*fz;

    f[d+0] := f[d+0] + a4*fx;
    f[d+1] := f[d+1] + a4*fy;
    f[d+2] := f[d+2] + a4*fz;
  END
END Compute;

BEGIN END EESpring.
