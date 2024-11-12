
#include <Collision.h>



#include <r3.h>
#include <BezierSearch.h>
#include <VertexData.h>
#include <EdgeData.h>
#include <FaceData.h>
#include <r3.h>


/*--- Deviations ------------------------------------------------------------*/

PROCEDURE ComputeVFDeviation(READONLY VertexData.Data p,
                             FaceData.Data *f);: 
                             VertexToFaceDeviation == 
VertexToFaceDeviation *?dev;
{
  { /* with*/ 
    q0 == Sub(p.p, f.pb),
    q1 == Sub(p.p, f.pa),
    dq0 == Sub(p.v, f.vb),
    dq1 == Sub(p.v, f.va)
  ) {
    dev.H = Dot(f.n,  q0);
    dev.G1 = Dot(f.n1, q0);
    dev.G2 = Dot(f.n2, q1);
    dev.G3 = Dot(f.n3, q1);
    dev.dH = Dot(f.dn,  q0) + Dot(f.n,  dq0);
    dev.dG1 = Dot(f.dn1, q0) + Dot(f.n1, dq0);
    dev.dG2 = Dot(f.dn2, q1) + Dot(f.n2, dq1);
    dev.dG3 = Dot(f.dn3, q1) + Dot(f.n3, dq1);
    return dev;
  };
} /* ComputeVFDeviation */;

PROCEDURE ComputeEEDeviation(READONLY e1,
                             EdgeData.Data e2): EdgeToEdgeDeviation == 
EdgeToEdgeDeviation *?dev;
{
  { /* with*/ 
    q0 == Sub(e1.pa, e2.pa),
    q1 == Sub(e1.pb, e2.pb),
    dq0 == Sub(e1.va, e2.va),
    dq1 == Sub(e1.vb, e2.vb),
    J0 == (r3_t){ e1.pb[1]*e2.pa[2] + e1.pb[2]*e2.pb[1] + e2.pa[1]*e2.pb[2] -
                e1.pb[1]*e2.pb[2] - e1.pb[2]*e2.pa[1] - e2.pa[2]*e2.pb[1],
               -e1.pb[0]*e2.pa[2] - e1.pb[2]*e2.pb[0] - e2.pa[0]*e2.pb[2] +
                e1.pb[0]*e2.pb[2] + e1.pb[2]*e2.pa[0] + e2.pa[2]*e2.pb[0],
                e1.pb[0]*e2.pa[1] + e1.pb[1]*e2.pb[0] + e2.pa[0]*e2.pb[1] -
                e1.pb[0]*e2.pb[1] - e1.pb[1]*e2.pa[0] - e2.pa[1]*e2.pb[0]},
    J1 == (r3_t){-e1.pa[1]*e2.pa[2] - e1.pa[2]*e2.pb[1] - e2.pa[1]*e2.pb[2] +
                e1.pa[1]*e2.pb[2] + e1.pa[2]*e2.pa[1] + e2.pa[2]*e2.pb[1],
                e1.pa[0]*e2.pa[2] + e1.pa[2]*e2.pb[0] + e2.pa[0]*e2.pb[2] -
                e1.pa[0]*e2.pb[2] - e1.pa[2]*e2.pa[0] - e2.pa[2]*e2.pb[0],
               -e1.pa[0]*e2.pa[1] - e1.pa[1]*e2.pb[0] - e2.pa[0]*e2.pb[1] +
                e1.pa[0]*e2.pb[1] + e1.pa[1]*e2.pa[0] + e2.pa[1]*e2.pb[0]},
    J2 == (r3_t){ e1.pa[1]*e1.pb[2] + e1.pa[2]*e2.pb[1] + e1.pb[1]*e2.pb[2] -
                e1.pa[1]*e2.pb[2] - e1.pa[2]*e1.pb[1] - e1.pb[2]*e2.pb[1],
               -e1.pa[0]*e1.pb[2] - e1.pa[2]*e2.pb[0] - e1.pb[0]*e2.pb[2] +
                e1.pa[0]*e2.pb[2] + e1.pa[2]*e1.pb[0] + e1.pb[2]*e2.pb[0],
                e1.pa[0]*e1.pb[1] + e1.pa[1]*e2.pb[0] + e1.pb[0]*e2.pb[1] -
                e1.pa[0]*e2.pb[1] - e1.pa[1]*e1.pb[0] - e1.pb[1]*e2.pb[0]},
    J3 == (r3_t){-e1.pa[1]*e1.pb[2] - e1.pa[2]*e2.pa[1] - e1.pb[1]*e2.pa[2] +
                e1.pa[1]*e2.pa[2] + e1.pa[2]*e1.pb[1] + e1.pb[2]*e2.pa[1],
                e1.pa[0]*e1.pb[2] + e1.pa[2]*e2.pa[0] + e1.pb[0]*e2.pa[2] -
                e1.pa[0]*e2.pa[2] - e1.pa[2]*e1.pb[0] - e1.pb[2]*e2.pa[0],
               -e1.pa[0]*e1.pb[1] - e1.pa[1]*e2.pa[0] - e1.pb[0]*e2.pa[1] +
                e1.pa[0]*e2.pa[1] + e1.pa[1]*e1.pb[0] + e1.pb[1]*e2.pa[0]}
  ) {
    dev.H = -(J0[0]*e1.pa[0] + J1[0]*e1.pb[0] +
                 J2[0]*e2.pa[0] + J3[0]*e2.pb[0]);
    dev.dH = -(Dot(J0, e1.va) + Dot(J1, e1.vb) +
                 Dot(J2, e2.va) + Dot(J3, e2.vb));
    dev.G1 = -Dot(e1.n, q0);
    dev.G2 = -Dot(e1.n, q1);
    dev.G3 = Dot(e2.n, q0);
    dev.G4 = Dot(e2.n, q1);
    dev.dG1 = -Dot(e1.dn, q0) - Dot(e1.n, dq0);
    dev.dG2 = -Dot(e1.dn, q1) - Dot(e1.n, dq1);
    dev.dG3 = Dot(e2.dn, q0) + Dot(e2.n, dq0);
    dev.dG4 = Dot(e2.dn, q1) + Dot(e2.n, dq1);
    return dev;
  };
} /* ComputeEEDeviation */;


/*--- Collision Detection ---------------------------------------------------*/

PROCEDURE DetectVFCollision(READONLY pa, VertexData.Data pb,
                            FaceData.Data *fa, fb;
                            double h3): double == 
{
  { /* with*/ 
    deva == ComputeVFDeviation(pa, fa),
    devb == ComputeVFDeviation(pb, fb),
    Hc == deva.H  + h3*deva.dH,
    Hd == devb.H  - h3*devb.dH,
    G1c == deva.G1 + h3*deva.dG1,
    G1d == devb.G1 - h3*devb.dG1,
    G2c == deva.G2 + h3*deva.dG2,
    G2d == devb.G2 - h3*devb.dG2,
    G3c == deva.G3 + h3*deva.dG3,
    G3d == devb.G3 - h3*devb.dG3,
    t == BezierSearch.FindRootA(deva.H,  Hc,  Hd,  devb.H,
                               deva.G1, G1c, G1d, devb.G1,
                               deva.G2, G2c, G2d, devb.G2,
                               deva.G3, G3c, G3d, devb.G3)
  ) {
    return t;
  };
} /* DetectVFCollision */;

PROCEDURE DetectEECollision(READONLY e1a, e1b, e2a, EdgeData.Data e2b,
                            double h3;
                            bool *?sign): double == 
{
  { /* with*/ 
    deva == ComputeEEDeviation(e1a, e2a),
    devb == ComputeEEDeviation(e1b, e2b),
    Hc == deva.H  + h3*deva.dH,
    Hd == devb.H  - h3*devb.dH,
    G1c == deva.G1 + h3*deva.dG1,
    G1d == devb.G1 - h3*devb.dG1,
    G2c == deva.G2 + h3*deva.dG2,
    G2d == devb.G2 - h3*devb.dG2,
    G3c == deva.G3 + h3*deva.dG3,
    G3d == devb.G3 - h3*devb.dG3,
    G4c == deva.G4 + h3*deva.dG4,
    G4d == devb.G4 - h3*devb.dG4,
    t == BezierSearch.FindRootB(deva.H,  Hc,  Hd,  devb.H,
                               deva.G1, G1c, G1d, devb.G1,
                               deva.G2, G2c, G2d, devb.G2,
                               deva.G3, G3c, G3d, devb.G3,
                               deva.G4, G4c, G4d, devb.G4)
  ) {
    sign = (deva.H > 0.0d0)) || ((devb.H < 0.0d0);
    return t;
  };
} /* DetectEECollision */;


/*--- Contact Break Detection -----------------------------------------------*/

PROCEDURE DetectVFContactBreak(READONLY pa, VertexData.Data pb,
                               FaceData.Data *fa, fb;
                               double h3): double == 
{
  { /* with*/ 
    deva == ComputeVFDeviation(pa, fa),
    devb == ComputeVFDeviation(pb, fb),
    Hc == deva.H + h3*deva.dH,
    Hd == devb.H - h3*devb.dH,
    t == BezierSearch.FindRoot(deva.H, Hc, Hd, devb.H, TRUE)
  ) {
    return t;
  };
} /* DetectVFContactBreak */;

PROCEDURE DetectEEContactBreak(READONLY e1a, e1b, e2a, EdgeData.Data e2b,
                               double h3; bool sign): double == 
{
  { /* with*/ 
    deva == ComputeEEDeviation(e1a, e2a),
    devb == ComputeEEDeviation(e1b, e2b),
    Hc == deva.H + h3*deva.dH,
    Hd == devb.H - h3*devb.dH,
    t == BezierSearch.FindRoot(deva.H, Hc, Hd, devb.H, sign)
  ) {
    return t;
  };
} /* DetectEEContactBreak */;

{; } /* Collision */.
