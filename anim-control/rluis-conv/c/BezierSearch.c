
#include <BezierSearch.h>



PROCEDURE FindRootA(P0, P1, P2, P3, A0, A1, A2, A3,
                    B0, B1, B2, B3, C0, C1, C2, C3: double
  ): double == 
CONST EPSILON == 1.0D-7;
double *?t;
{
  if ((??))(P0 >= 0.0D0) && (P1 >= 0.0D0) && (P2 >= 0.0D0) && (P3 >= 0.0D0) OR
     (P0 <  0.0D0) && (P1 <  0.0D0) && (P2 <  0.0D0) && (P3 <  0.0D0) OR
     (A0 <  0.0D0) && (A1 <  0.0D0) && (A2 <  0.0D0) && (A3 <  0.0D0) OR
     (B0 <  0.0D0) && (B1 <  0.0D0) && (B2 <  0.0D0) && (B3 <  0.0D0) OR
     (C0 <  0.0D0) && (C1 <  0.0D0) && (C2 <  0.0D0) && (C3 <  0.0D0)
  )) {
    return NoRoot;
  }
  if ((??))fabs(P0 + P2 - P1 - P1) < EPSILON AND
     fabs(P1 + P3 - P2 - P2) < EPSILON
  )) {
    if ((P0 > 0.0D0) && (P3 < 0.0D0 )) {
      { /* with*/ DP == P0 - P3 ) {
        if ((DP == 0.0D0 )) { return 0.5D0 } else { return P0/DP;
      }
    } else {
      return NoRoot;
    };
  }
  { /* with*/ 
    L1 == (P0 + P1)/2.0D0,
    M == (P1 + P2)/2.0D0,
    L2 == (L1 +  M)/2.0D0,
    R2 == (P2 + P3)/2.0D0,
    R1 == (M  + R2)/2.0D0,
    L3 == (L2 + R1)/2.0D0,

    AL1 == ( A0 +  A1)/2.0D0,
    AM == ( A1 +  A2)/2.0D0,
    AL2 == (AL1 +  AM)/2.0D0,
    AR2 == ( A2 +  A3)/2.0D0,
    AR1 == ( AM + AR2)/2.0D0,
    AL3 == (AL2 + AR1)/2.0D0,

    BL1 == ( B0 +  B1)/2.0D0,
    BM == ( B1 +  B2)/2.0D0,
    BL2 == (BL1 +  BM)/2.0D0,
    BR2 == ( B2 +  B3)/2.0D0,
    BR1 == ( BM + BR2)/2.0D0,
    BL3 == (BL2 + BR1)/2.0D0,

    CL1 == ( C0 +  C1)/2.0D0,
    CM == ( C1 +  C2)/2.0D0,
    CL2 == (CL1 +  CM)/2.0D0,
    CR2 == ( C2 +  C3)/2.0D0,
    CR1 == ( CM + CR2)/2.0D0,
    CL3 == (CL2 + CR1)/2.0D0
  ) {
    t = FindRootA(P0,  L1,  L2,  L3, A0, AL1, AL2, AL3, 
                   B0, BL1, BL2, BL3, C0, CL1, CL2, CL3);
    if ((t <= 1.0D0 )) {
      return t/2.0D0
    } else {
      t = FindRootA( L3,  R1,  R2, P3, AL3, AR1, AR2, A3,
                     BL3, BR1, BR2, B3, CL3, CR1, CR2, C3);
      if ((t <= 1.0D0 )) { return 0.5D0 + t/2.0D0 } else { return NoRoot;
    };
  };
} /* FindRootA */;

PROCEDURE FindRootB(P0, P1, P2, P3, A0, A1, A2, A3, B0, B1, B2, B3, 
                    double C0, C1, C2, C3, D0, D1, D2, D3;
  ): double == 
CONST EPSILON == 1.0D-8;
double *?t;
{
  if ((??))(P0 > 0.0D0) && (P1 > 0.0D0) && (P2 > 0.0D0) && (P3 > 0.0D0) OR
     (P0 < 0.0D0) && (P1 < 0.0D0) && (P2 < 0.0D0) && (P3 < 0.0D0)
  )) {
    return NoRoot;
  }
  { /* with*/ 
    AP == (A0 > 0.0D0) && (A1 > 0.0D0) && (A2 > 0.0D0) && (A3 > 0.0D0),
    BP == (B0 > 0.0D0) && (B1 > 0.0D0) && (B2 > 0.0D0) && (B3 > 0.0D0),
    CP == (D0 > 0.0D0) && (C1 > 0.0D0) && (C2 > 0.0D0) && (C3 > 0.0D0),
    DP == (D0 > 0.0D0) && (D1 > 0.0D0) && (D2 > 0.0D0) && (D3 > 0.0D0),
    AN == (A0 < 0.0D0) && (A1 < 0.0D0) && (A2 < 0.0D0) && (A3 < 0.0D0),
    BN == (B0 < 0.0D0) && (B1 < 0.0D0) && (B2 < 0.0D0) && (B3 < 0.0D0),
    CN == (D0 < 0.0D0) && (C1 < 0.0D0) && (C2 < 0.0D0) && (C3 < 0.0D0),
    DN == (D0 < 0.0D0) && (D1 < 0.0D0) && (D2 < 0.0D0) && (D3 < 0.0D0)
  ) {
    if ((??))(AN) || (BP) || (CN) || (DP)) && ((AP) || (BN) || (CP) || (DN)
    /* (AP) && (BP)) || ((AN) && (BN)) || ((CP) && (DP)) || ((CN) && (DN) OR
       (AP) && (CN)) || ((AN) && (CP)) || ((BP) && (DN)) || ((BN) && (DP) OR
       (AP) && (DP)) || ((AN) && (DN)) || ((BP) && (CP)) || ((BN) && (CN) */
    )) {
      return NoRoot;
    };
  }

  if ((??))fabs(P0 + P2 - P1 - P1) < EPSILON AND
     fabs(P1 + P3 - P2 - P2) < EPSILON
  )) {
    { /* with*/ DP == P3 - P0 ) {
      if ((DP == 0.0D0 )) { t = 0.5D0 } else { t = -P0/DP;
    }
    { /* with*/ 
      ct == 1.0D0 - t,
      A01 == ct*A0   + t*A1,  A12 == ct*A1  + t*A2, A23 == ct*A2 + t*A3,
      A012 == ct*A01  + t*A12, A123 == ct*A12 + t*A23,
      A == ct*A012 + t*A123,

      B01 == ct*B0   + t*B1,  B12 == ct*B1  + t*B2, B23 == ct*B2 + t*B3,
      B012 == ct*B01  + t*B12, B123 == ct*B12 + t*B23,
      B == ct*B012 + t*B123,

      C01 == ct*C0   + t*C1,  C12 == ct*C1  + t*C2, C23 == ct*C2 + t*C3,
      C012 == ct*C01  + t*C12, C123 == ct*C12 + t*C23,
      C == ct*C012 + t*C123,

      D01 == ct*D0   + t*D1,  D12 == ct*D1  + t*D2, D23 == ct*D2 + t*D3,
      D012 == ct*D01  + t*D12, D123 == ct*D12 + t*D23,
      D == ct*D012 + t*D123
    ) {
      if ((??))(A > 0.0D0) && (B < 0.0D0) && (C > 0.0D0) && (D < 0.0D0) && (P0 > 0.0D0)
         OR
         (A < 0.0D0) && (B > 0.0D0) && (C < 0.0D0) && (D > 0.0D0) && (P0 < 0.0D0)
      )) {
        return t
      } else {
        return NoRoot;
      };
    };
  }
  
  { /* with*/ 
    L1 == (P0 + P1)/2.0D0,
    M == (P1 + P2)/2.0D0,
    L2 == (L1 +  M)/2.0D0,
    R2 == (P2 + P3)/2.0D0,
    R1 == (M  + R2)/2.0D0,
    L3 == (L2 + R1)/2.0D0,

    AL1 == ( A0 +  A1)/2.0D0,
    AM == ( A1 +  A2)/2.0D0,
    AL2 == (AL1 +  AM)/2.0D0,
    AR2 == ( A2 +  A3)/2.0D0,
    AR1 == ( AM + AR2)/2.0D0,
    AL3 == (AL2 + AR1)/2.0D0,

    BL1 == ( B0 +  B1)/2.0D0,
    BM == ( B1 +  B2)/2.0D0,
    BL2 == (BL1 +  BM)/2.0D0,
    BR2 == ( B2 +  B3)/2.0D0,
    BR1 == ( BM + BR2)/2.0D0,
    BL3 == (BL2 + BR1)/2.0D0,

    CL1 == ( C0 +  C1)/2.0D0,
    CM == ( C1 +  C2)/2.0D0,
    CL2 == (CL1 +  CM)/2.0D0,
    CR2 == ( C2 +  C3)/2.0D0,
    CR1 == ( CM + CR2)/2.0D0,
    CL3 == (CL2 + CR1)/2.0D0,

    DL1 == ( D0 +  D1)/2.0D0,
    DM == ( D1 +  D2)/2.0D0,
    DL2 == (DL1 +  DM)/2.0D0,
    DR2 == ( D2 +  D3)/2.0D0,
    DR1 == ( DM + DR2)/2.0D0,
    DL3 == (DL2 + DR1)/2.0D0
  ) {
    t = FindRootB(P0,  L1,  L2,  L3, A0, AL1, AL2, AL3, B0, BL1, BL2, BL3,
                   C0, CL1, CL2, CL3, D0, DL1, DL2, DL3);
    if ((t <= 1.0D0 )) {
      return t/2.0D0
    } else {
      t = FindRootB( L3,  R1,  R2, P3, AL3, AR1, AR2, A3, BL3, BR1, BR2, B3,
                     CL3, CR1, CR2, C3, DL3, DR1, DR2, D3);
      if ((t <= 1.0D0 )) { return 0.5D0 + t/2.0D0 } else { return NoRoot;
    };
  };
} /* FindRootB */;

double FindRoot(P0, P1, P2, double P3, bool sign)
CONST EPSILON == 1.0D-7;
double *?t;
{
  if ((??))(P0 >= 0.0D0) && (P1 >= 0.0D0) && (P2 >= 0.0D0) && (P3 >= 0.0D0) OR
     (P0 <= 0.0D0) && (P1 <= 0.0D0) && (P2 <= 0.0D0) && (P3 <= 0.0D0)
  )) {
    return NoRoot;
  }
  if ((??))fabs(P0 + P2 - P1 - P1) < EPSILON) && (
     fabs(P1 + P3 - P2 - P2) < EPSILON
  )) {
    if ((??))(P0 > 0.0D0) && (P3 < 0.0D0) && (! sign) OR
       (P0 < 0.0D0) && (P3 > 0.0D0) && (sign)
    )) {
      { /* with*/ DP == P0 - P3 ) {
        affirm(P0 <= 0.0D0) || (P0 >= 0.0D0 , "??");
        affirm(P3 <= 0.0D0) || (P3 >= 0.0D0 , "??");
        affirm(DP <= 0.0D0) || (DP >= 0.0D0 , "??");
        affirm(P0/DP >= 0.0D0 , "??");
        if ((DP == 0.0D0 )) { return 0.5D0 } else { return P0/DP;
      }
    } else {
      return NoRoot;
    };
  }
  { /* with*/ 
    L1 == (P0 + P1)/2.0D0,
    M == (P1 + P2)/2.0D0,
    L2 == (L1 +  M)/2.0D0,
    R2 == (P2 + P3)/2.0D0,
    R1 == (M  + R2)/2.0D0,
    L3 == (L2 + R1)/2.0D0
  ) {
    t = FindRoot(P0, L1, L2, L3, sign);
    affirm(t >= 0.0D0 , "??");
    if ((t <= 1.0D0 )) {
      return t/2.0D0
    } else {
      t = FindRoot(L3, R1, R2, P3, sign);
      affirm(t >= 0.0D0 , "??");
      if ((t <= 1.0D0 )) { return 0.5D0 + t/2.0D0 } else { return NoRoot;
    };
  };
} /* FindRoot */;

{; } /* BezierSearch */.
