MODULE IncompressibleConstraint;

FROM Constraint IMPORT Vector;

REVEAL
  T = Public BRANDED OBJECT
      na, nb, nc, nd: CARDINAL;  (* Vertex indices times 3 *)
      gxa, gya, gza,
      gxb, gyb, gzb,
      gxc, gyc, gzc,
      gxd, gyd, gzd,
      psis: LONGREAL;
    OVERRIDES
      init         := Init;
      computeGrad  := ComputeGrad;
      computePsi   := ComputePsi;
      getGrad       := GetGrad;
      dirDerivative     := DirDerivative;
      getPsi          := GetPsi;
      setPsi       := SetPsi;
      addReactionForce := AddReactionForce;
    END;

PROCEDURE Init(q: T; na, nb, nc, nd: CARDINAL) =
BEGIN
  q.na := 3*na; q.nb := 3*nb; q.nc := 3*nc; q.nd := 3*nd;
END Init;

PROCEDURE ComputeGrad(<*UNUSED*> q: T; <*UNUSED*> pos: Vector;
                      <*UNUSED*> t: LONGREAL) =
BEGIN
  (* It lets all to be computed in ComputePsi. There will be no *
   * problem with this in this latest implementation.           *)
END ComputeGrad;

PROCEDURE ComputePsi(q: T; pos, vel: Vector; <*UNUSED*> t: LONGREAL) =
BEGIN
  (*
    This code could be improved by computing first the 2x2 minors,
    then deriving the 3x3 minors from them.
  *)
  WITH
    na = q.na, nb = q.nb, nc = q.nc, nd = q.nd,
    pxa = pos[na], pya = pos[na+1], pza = pos[na+2],
    pxb = pos[nb], pyb = pos[nb+1], pzb = pos[nb+2],
    pxc = pos[nc], pyc = pos[nc+1], pzc = pos[nc+2],
    pxd = pos[nd], pyd = pos[nd+1], pzd = pos[nd+2],
    vxa = vel[na], vya = vel[na+1], vza = vel[na+2],
    vxb = vel[nb], vyb = vel[nb+1], vzb = vel[nb+2],
    vxc = vel[nc], vyc = vel[nc+1], vzc = vel[nc+2],
    vxd = vel[nd], vyd = vel[nd+1], vzd = vel[nd+2],
    
    gxayb = pzd - pzc, gxazb = pyc - pyd, gxayc = pzb - pzd,
    gxazc = pyd - pyb, gxayd = pzc - pzb, gxazd = pyb - pyc,
    
    gyaxb = pzc - pzd, gyazb = pxd - pxc, gyaxc = pzd - pzb,
    gyazc = pxb - pxd, gyaxd = pzb - pzc,  gyazd = pxc - pxb,
    
    gzaxb = pyd - pyc, gzayb = pxc - pxd, gzaxc = pyb - pyd,
    gzayc = pxd - pxb, gzaxd = pyc - pyb, gzayd = pxb - pxc,
    
    gxbyc = pzd - pza, gxbzc = pya - pyd, gxbyd = pza - pzc,  gxbzd = pyc - pya,
    
    gybxc = pza - pzd, gybzc = pxd - pxa, gybxd = pzc - pza, gybzd = pxa - pxc,
    
    gzbxc = pyd - pya, gzbyc = pxa - pxd, gzbxd = pya - pyc, gzbyd = pxc - pxa,
    
    gxcyd = pzb - pza, gxczd = pya - pyb,

    gycxd = pza - pzb, gyczd = pxb - pxa,
    
    gzcxd = pyb - pya, gzcyd = pxa - pxb,
    
    rxa = gxayb*vyb + gxazb*vzb + gxayc*vyc + gxazc*vzc + gxayd*vyd + gxazd*vzd,
    rya = gyaxb*vxb + gyazb*vzb + gyaxc*vxc + gyazc*vzc + gyaxd*vxd + gyazd*vzd,
    rza = gzaxb*vxb + gzayb*vyb + gzaxc*vxc + gzayc*vyc + gzaxd*vxd + gzayd*vyd,
    rxb = gyaxb*vya + gzaxb*vza + gxbyc*vyc + gxbzc*vzc + gxbyd*vyd + gxbzd*vzd,
    ryb = gxayb*vxa + gzayb*vza + gybxc*vxc + gybzc*vzc + gybxd*vxd + gybzd*vzd,
    rzb = gxazb*vxa + gyazb*vya + gzbxc*vxc + gzbyc*vyc + gzbxd*vxd + gzbyd*vyd,
    rxc = gyaxc*vya + gzaxc*vza + gybxc*vyb + gzbxc*vzb + gxcyd*vyd + gxczd*vzd,
    ryc = gxayc*vxa + gzayc*vza + gxbyc*vxb + gzbyc*vzb + gycxd*vxd + gyczd*vzd,
    rzc = gxazc*vxa + gyazc*vya + gxbzc*vxb + gybzc*vyb + gzcxd*vxd + gzcyd*vyd,
    rxd = gyaxd*vya + gzaxd*vza + gybxd*vyb + gzbxd*vzb + gycxd*vyc + gzcxd*vzc,
    ryd = gxayd*vxa + gzayd*vza + gxbyd*vxb + gzbyd*vzb + gxcyd*vxc + gzcyd*vzc,
    rzd = gxazd*vxa + gyazd*vya + gxbzd*vxb + gybzd*vyb + gxczd*vxc + gyczd*vyc
  DO
    q.gxa :=  pyb*gxayb + pyc*gxayc + pyd*gxayd;
    q.gya :=  pxb*gyaxb + pxc*gyaxc + pxd*gyaxd;
    q.gza := -pxb*gxazb + pxc*gxazc - pxd*gxazd;
    
    q.gxb :=  pya*gyaxb + pyc*gxbyc + pyd*gxbyd;
    q.gyb :=  pxa*gxayb + pxc*gybxc + pxd*gybxd;
    q.gzb :=  pxa*gxazb - pxc*gxbzc - pxd*gxbzd;
    
    q.gxc :=  pya*gyaxc + pyb*gybxc + pyd*gxcyd;
    q.gyc :=  pxa*gxayc + pxb*gxbyc + pxd*gycxd;
    q.gzc := -pxa*gxazc + pxb*gxbzc - pxd*gxczd;
    
    q.gxd := pya*gyaxd + pyb*gzbxd + pyc*gycxd;
    q.gyd := pxa*gxayd + pxb*gxbyd + pxc*gxcyd;
    q.gzd := pxa*gxazd + pxb*gxbzd + pxc*gxczd;
    
    q.psis := vxa*rxa + vya*rya + vza*rza +
              vxb*rxb + vyb*ryb + vzb*rzb +
              vxc*rxc + vyc*ryc + vzc*rzc +
              vxd*rxd + vyd*ryd + vzd*rzd;
  END
END ComputePsi;

PROCEDURE GetGrad(q: T; grad: Vector) =
BEGIN
  grad[q.na] := q.gxa; grad[q.na+1] := q.gya; grad[q.na+2] := q.gza;
  grad[q.nb] := q.gxb; grad[q.nb+1] := q.gyb; grad[q.nb+2] := q.gzb;
  grad[q.nc] := q.gxc; grad[q.nc+1] := q.gyc; grad[q.nc+2] := q.gzc;
  grad[q.nd] := q.gxd; grad[q.nd+1] := q.gyd; grad[q.nd+2] := q.gzd;
END GetGrad;

PROCEDURE DirDerivative(q: T; u: Vector): LONGREAL =
BEGIN
  RETURN q.gxa*u[q.na] + q.gya*u[q.na+1] + q.gza*u[q.na+2] +
         q.gxb*u[q.nb] + q.gyb*u[q.nb+1] + q.gzb*u[q.nb+2] +
         q.gxc*u[q.nc] + q.gyc*u[q.nc+1] + q.gzc*u[q.nc+2] +
         q.gxd*u[q.nd] + q.gyd*u[q.nd+1] + q.gzd*u[q.nd+2]
END DirDerivative;

PROCEDURE GetPsi(q: T; <*UNUSED*> j: CARDINAL): LONGREAL =
BEGIN 
  RETURN q.psis
END GetPsi;

PROCEDURE SetPsi(q: T; <*UNUSED*> j: CARDINAL; psi: LONGREAL) =
BEGIN
  q.psis := psi
END SetPsi;

PROCEDURE AddReactionForce(q: T; lambda: LONGREAL; f: Vector) =
BEGIN
  WITH na = q.na, nb = q.nb, nc = q.nc, nd = q.nd DO
    f[na+0] := f[na+0] - q.gxa*lambda;
    f[na+1] := f[na+1] - q.gya*lambda;
    f[na+2] := f[na+2] - q.gza*lambda;
    f[nb+0] := f[nb+0] - q.gxb*lambda;
    f[nb+1] := f[nb+1] - q.gyb*lambda;
    f[nb+2] := f[nb+2] - q.gzb*lambda;
    f[nc+0] := f[nc+0] - q.gxc*lambda;
    f[nc+1] := f[nc+1] - q.gyc*lambda;
    f[nc+2] := f[nc+2] - q.gzc*lambda;
    f[nd+0] := f[nd+0] - q.gxd*lambda;
    f[nd+1] := f[nd+1] - q.gyd*lambda;
    f[nd+2] := f[nd+2] - q.gzd*lambda;
  END
END AddReactionForce;

BEGIN END IncompressibleConstraint.
