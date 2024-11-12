MODULE Collision;

IMPORT LR3, BezierSearch, VertexData, EdgeData, FaceData;
FROM LR3 IMPORT Sub, Dot;


(*--- Deviations ------------------------------------------------------------*)

PROCEDURE ComputeVFDeviation(READONLY p: VertexData.Data;
                             READONLY f: FaceData.Data): 
                             VertexToFaceDeviation =
VAR dev: VertexToFaceDeviation;
BEGIN
  WITH
    q0  = Sub(p.p, f.pb),
    q1  = Sub(p.p, f.pa),
    dq0 = Sub(p.v, f.vb),
    dq1 = Sub(p.v, f.va)
  DO
    dev.H   := Dot(f.n,  q0);
    dev.G1  := Dot(f.n1, q0);
    dev.G2  := Dot(f.n2, q1);
    dev.G3  := Dot(f.n3, q1);
    dev.dH  := Dot(f.dn,  q0) + Dot(f.n,  dq0);
    dev.dG1 := Dot(f.dn1, q0) + Dot(f.n1, dq0);
    dev.dG2 := Dot(f.dn2, q1) + Dot(f.n2, dq1);
    dev.dG3 := Dot(f.dn3, q1) + Dot(f.n3, dq1);
    RETURN dev
  END
END ComputeVFDeviation;

PROCEDURE ComputeEEDeviation(READONLY e1,
                             e2: EdgeData.Data): EdgeToEdgeDeviation =
VAR dev: EdgeToEdgeDeviation;
BEGIN
  WITH
    q0  = Sub(e1.pa, e2.pa),
    q1  = Sub(e1.pb, e2.pb),
    dq0 = Sub(e1.va, e2.va),
    dq1 = Sub(e1.vb, e2.vb),
    J0 = LR3.T{ e1.pb[1]*e2.pa[2] + e1.pb[2]*e2.pb[1] + e2.pa[1]*e2.pb[2] -
                e1.pb[1]*e2.pb[2] - e1.pb[2]*e2.pa[1] - e2.pa[2]*e2.pb[1],
               -e1.pb[0]*e2.pa[2] - e1.pb[2]*e2.pb[0] - e2.pa[0]*e2.pb[2] +
                e1.pb[0]*e2.pb[2] + e1.pb[2]*e2.pa[0] + e2.pa[2]*e2.pb[0],
                e1.pb[0]*e2.pa[1] + e1.pb[1]*e2.pb[0] + e2.pa[0]*e2.pb[1] -
                e1.pb[0]*e2.pb[1] - e1.pb[1]*e2.pa[0] - e2.pa[1]*e2.pb[0]},
    J1 = LR3.T{-e1.pa[1]*e2.pa[2] - e1.pa[2]*e2.pb[1] - e2.pa[1]*e2.pb[2] +
                e1.pa[1]*e2.pb[2] + e1.pa[2]*e2.pa[1] + e2.pa[2]*e2.pb[1],
                e1.pa[0]*e2.pa[2] + e1.pa[2]*e2.pb[0] + e2.pa[0]*e2.pb[2] -
                e1.pa[0]*e2.pb[2] - e1.pa[2]*e2.pa[0] - e2.pa[2]*e2.pb[0],
               -e1.pa[0]*e2.pa[1] - e1.pa[1]*e2.pb[0] - e2.pa[0]*e2.pb[1] +
                e1.pa[0]*e2.pb[1] + e1.pa[1]*e2.pa[0] + e2.pa[1]*e2.pb[0]},
    J2 = LR3.T{ e1.pa[1]*e1.pb[2] + e1.pa[2]*e2.pb[1] + e1.pb[1]*e2.pb[2] -
                e1.pa[1]*e2.pb[2] - e1.pa[2]*e1.pb[1] - e1.pb[2]*e2.pb[1],
               -e1.pa[0]*e1.pb[2] - e1.pa[2]*e2.pb[0] - e1.pb[0]*e2.pb[2] +
                e1.pa[0]*e2.pb[2] + e1.pa[2]*e1.pb[0] + e1.pb[2]*e2.pb[0],
                e1.pa[0]*e1.pb[1] + e1.pa[1]*e2.pb[0] + e1.pb[0]*e2.pb[1] -
                e1.pa[0]*e2.pb[1] - e1.pa[1]*e1.pb[0] - e1.pb[1]*e2.pb[0]},
    J3 = LR3.T{-e1.pa[1]*e1.pb[2] - e1.pa[2]*e2.pa[1] - e1.pb[1]*e2.pa[2] +
                e1.pa[1]*e2.pa[2] + e1.pa[2]*e1.pb[1] + e1.pb[2]*e2.pa[1],
                e1.pa[0]*e1.pb[2] + e1.pa[2]*e2.pa[0] + e1.pb[0]*e2.pa[2] -
                e1.pa[0]*e2.pa[2] - e1.pa[2]*e1.pb[0] - e1.pb[2]*e2.pa[0],
               -e1.pa[0]*e1.pb[1] - e1.pa[1]*e2.pa[0] - e1.pb[0]*e2.pa[1] +
                e1.pa[0]*e2.pa[1] + e1.pa[1]*e1.pb[0] + e1.pb[1]*e2.pa[0]}
  DO
    dev.H   := -(J0[0]*e1.pa[0] + J1[0]*e1.pb[0] +
                 J2[0]*e2.pa[0] + J3[0]*e2.pb[0]);
    dev.dH  := -(Dot(J0, e1.va) + Dot(J1, e1.vb) +
                 Dot(J2, e2.va) + Dot(J3, e2.vb));
    dev.G1  := -Dot(e1.n, q0);
    dev.G2  := -Dot(e1.n, q1);
    dev.G3  :=  Dot(e2.n, q0);
    dev.G4  :=  Dot(e2.n, q1);
    dev.dG1 := -Dot(e1.dn, q0) - Dot(e1.n, dq0);
    dev.dG2 := -Dot(e1.dn, q1) - Dot(e1.n, dq1);
    dev.dG3 :=  Dot(e2.dn, q0) + Dot(e2.n, dq0);
    dev.dG4 :=  Dot(e2.dn, q1) + Dot(e2.n, dq1);
    RETURN dev
  END
END ComputeEEDeviation;


(*--- Collision Detection ---------------------------------------------------*)

PROCEDURE DetectVFCollision(READONLY pa, pb: VertexData.Data;
                            READONLY fa, fb: FaceData.Data;
                            h3: LONGREAL): LONGREAL =
BEGIN
  WITH
    deva = ComputeVFDeviation(pa, fa),
    devb = ComputeVFDeviation(pb, fb),
    Hc  = deva.H  + h3*deva.dH,
    Hd  = devb.H  - h3*devb.dH,
    G1c = deva.G1 + h3*deva.dG1,
    G1d = devb.G1 - h3*devb.dG1,
    G2c = deva.G2 + h3*deva.dG2,
    G2d = devb.G2 - h3*devb.dG2,
    G3c = deva.G3 + h3*deva.dG3,
    G3d = devb.G3 - h3*devb.dG3,
    t = BezierSearch.FindRootA(deva.H,  Hc,  Hd,  devb.H,
                               deva.G1, G1c, G1d, devb.G1,
                               deva.G2, G2c, G2d, devb.G2,
                               deva.G3, G3c, G3d, devb.G3)
  DO
    RETURN t
  END
END DetectVFCollision;

PROCEDURE DetectEECollision(READONLY e1a, e1b, e2a, e2b: EdgeData.Data;
                            h3: LONGREAL;
                            VAR sign: BOOLEAN): LONGREAL =
BEGIN
  WITH
    deva = ComputeEEDeviation(e1a, e2a),
    devb = ComputeEEDeviation(e1b, e2b),
    Hc  = deva.H  + h3*deva.dH,
    Hd  = devb.H  - h3*devb.dH,
    G1c = deva.G1 + h3*deva.dG1,
    G1d = devb.G1 - h3*devb.dG1,
    G2c = deva.G2 + h3*deva.dG2,
    G2d = devb.G2 - h3*devb.dG2,
    G3c = deva.G3 + h3*deva.dG3,
    G3d = devb.G3 - h3*devb.dG3,
    G4c = deva.G4 + h3*deva.dG4,
    G4d = devb.G4 - h3*devb.dG4,
    t = BezierSearch.FindRootB(deva.H,  Hc,  Hd,  devb.H,
                               deva.G1, G1c, G1d, devb.G1,
                               deva.G2, G2c, G2d, devb.G2,
                               deva.G3, G3c, G3d, devb.G3,
                               deva.G4, G4c, G4d, devb.G4)
  DO
    sign := (deva.H > 0.0d0) OR (devb.H < 0.0d0);
    RETURN t
  END
END DetectEECollision;


(*--- Contact Break Detection -----------------------------------------------*)

PROCEDURE DetectVFContactBreak(READONLY pa, pb: VertexData.Data;
                               READONLY fa, fb: FaceData.Data;
                               h3: LONGREAL): LONGREAL =
BEGIN
  WITH
    deva = ComputeVFDeviation(pa, fa),
    devb = ComputeVFDeviation(pb, fb),
    Hc = deva.H + h3*deva.dH,
    Hd = devb.H - h3*devb.dH,
    t = BezierSearch.FindRoot(deva.H, Hc, Hd, devb.H, TRUE)
  DO
    RETURN t
  END
END DetectVFContactBreak;

PROCEDURE DetectEEContactBreak(READONLY e1a, e1b, e2a, e2b: EdgeData.Data;
                               h3: LONGREAL; sign: BOOLEAN): LONGREAL =
BEGIN
  WITH
    deva = ComputeEEDeviation(e1a, e2a),
    devb = ComputeEEDeviation(e1b, e2b),
    Hc = deva.H + h3*deva.dH,
    Hd = devb.H - h3*devb.dH,
    t = BezierSearch.FindRoot(deva.H, Hc, Hd, devb.H, sign)
  DO
    RETURN t
  END
END DetectEEContactBreak;

BEGIN END Collision.
