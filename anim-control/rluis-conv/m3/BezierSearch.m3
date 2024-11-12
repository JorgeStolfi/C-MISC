MODULE BezierSearch;

PROCEDURE FindRootA(P0, P1, P2, P3, A0, A1, A2, A3,
                    B0, B1, B2, B3, C0, C1, C2, C3: LONGREAL
  ): LONGREAL =
CONST EPSILON = 1.0D-7;
VAR t: LONGREAL;
BEGIN
  IF (P0 >= 0.0D0 AND P1 >= 0.0D0 AND P2 >= 0.0D0 AND P3 >= 0.0D0) OR
     (P0 <  0.0D0 AND P1 <  0.0D0 AND P2 <  0.0D0 AND P3 <  0.0D0) OR
     (A0 <  0.0D0 AND A1 <  0.0D0 AND A2 <  0.0D0 AND A3 <  0.0D0) OR
     (B0 <  0.0D0 AND B1 <  0.0D0 AND B2 <  0.0D0 AND B3 <  0.0D0) OR
     (C0 <  0.0D0 AND C1 <  0.0D0 AND C2 <  0.0D0 AND C3 <  0.0D0)
  THEN
    RETURN NoRoot
  END;
  IF ABS(P0 + P2 - P1 - P1) < EPSILON AND
     ABS(P1 + P3 - P2 - P2) < EPSILON
  THEN
    IF P0 > 0.0D0 AND P3 < 0.0D0 THEN
      WITH DP = P0 - P3 DO
        IF DP = 0.0D0 THEN RETURN 0.5D0 ELSE RETURN P0/DP END
      END
    ELSE
      RETURN NoRoot
    END
  END;
  WITH
    L1 = (P0 + P1)/2.0D0,
    M  = (P1 + P2)/2.0D0,
    L2 = (L1 +  M)/2.0D0,
    R2 = (P2 + P3)/2.0D0,
    R1 = (M  + R2)/2.0D0,
    L3 = (L2 + R1)/2.0D0,

    AL1 = ( A0 +  A1)/2.0D0,
    AM  = ( A1 +  A2)/2.0D0,
    AL2 = (AL1 +  AM)/2.0D0,
    AR2 = ( A2 +  A3)/2.0D0,
    AR1 = ( AM + AR2)/2.0D0,
    AL3 = (AL2 + AR1)/2.0D0,

    BL1 = ( B0 +  B1)/2.0D0,
    BM  = ( B1 +  B2)/2.0D0,
    BL2 = (BL1 +  BM)/2.0D0,
    BR2 = ( B2 +  B3)/2.0D0,
    BR1 = ( BM + BR2)/2.0D0,
    BL3 = (BL2 + BR1)/2.0D0,

    CL1 = ( C0 +  C1)/2.0D0,
    CM  = ( C1 +  C2)/2.0D0,
    CL2 = (CL1 +  CM)/2.0D0,
    CR2 = ( C2 +  C3)/2.0D0,
    CR1 = ( CM + CR2)/2.0D0,
    CL3 = (CL2 + CR1)/2.0D0
  DO
    t := FindRootA(P0,  L1,  L2,  L3, A0, AL1, AL2, AL3, 
                   B0, BL1, BL2, BL3, C0, CL1, CL2, CL3);
    IF t <= 1.0D0 THEN
      RETURN t/2.0D0
    ELSE
      t := FindRootA( L3,  R1,  R2, P3, AL3, AR1, AR2, A3,
                     BL3, BR1, BR2, B3, CL3, CR1, CR2, C3);
      IF t <= 1.0D0 THEN RETURN 0.5D0 + t/2.0D0 ELSE RETURN NoRoot END
    END
  END
END FindRootA;

PROCEDURE FindRootB(P0, P1, P2, P3, A0, A1, A2, A3, B0, B1, B2, B3, 
                    C0, C1, C2, C3, D0, D1, D2, D3: LONGREAL;
  ): LONGREAL =
CONST EPSILON = 1.0D-8;
VAR t: LONGREAL;
BEGIN
  IF (P0 > 0.0D0 AND P1 > 0.0D0 AND P2 > 0.0D0 AND P3 > 0.0D0) OR
     (P0 < 0.0D0 AND P1 < 0.0D0 AND P2 < 0.0D0 AND P3 < 0.0D0)
  THEN
    RETURN NoRoot
  END;
  WITH
    AP = (A0 > 0.0D0 AND A1 > 0.0D0 AND A2 > 0.0D0 AND A3 > 0.0D0),
    BP = (B0 > 0.0D0 AND B1 > 0.0D0 AND B2 > 0.0D0 AND B3 > 0.0D0),
    CP = (D0 > 0.0D0 AND C1 > 0.0D0 AND C2 > 0.0D0 AND C3 > 0.0D0),
    DP = (D0 > 0.0D0 AND D1 > 0.0D0 AND D2 > 0.0D0 AND D3 > 0.0D0),
    AN = (A0 < 0.0D0 AND A1 < 0.0D0 AND A2 < 0.0D0 AND A3 < 0.0D0),
    BN = (B0 < 0.0D0 AND B1 < 0.0D0 AND B2 < 0.0D0 AND B3 < 0.0D0),
    CN = (D0 < 0.0D0 AND C1 < 0.0D0 AND C2 < 0.0D0 AND C3 < 0.0D0),
    DN = (D0 < 0.0D0 AND D1 < 0.0D0 AND D2 < 0.0D0 AND D3 < 0.0D0)
  DO
    IF (AN OR BP OR CN OR DP) AND (AP OR BN OR CP OR DN)
    (* (AP AND BP) OR (AN AND BN) OR (CP AND DP) OR (CN AND DN) OR
       (AP AND CN) OR (AN AND CP) OR (BP AND DN) OR (BN AND DP) OR
       (AP AND DP) OR (AN AND DN) OR (BP AND CP) OR (BN AND CN) *)
    THEN
      RETURN NoRoot
    END
  END;

  IF ABS(P0 + P2 - P1 - P1) < EPSILON AND
     ABS(P1 + P3 - P2 - P2) < EPSILON
  THEN
    WITH DP = P3 - P0 DO
      IF DP = 0.0D0 THEN t := 0.5D0 ELSE t := -P0/DP END
    END;
    WITH
      ct   = 1.0D0 - t,
      A01  = ct*A0   + t*A1,  A12  = ct*A1  + t*A2, A23 = ct*A2 + t*A3,
      A012 = ct*A01  + t*A12, A123 = ct*A12 + t*A23,
      A    = ct*A012 + t*A123,

      B01  = ct*B0   + t*B1,  B12  = ct*B1  + t*B2, B23 = ct*B2 + t*B3,
      B012 = ct*B01  + t*B12, B123 = ct*B12 + t*B23,
      B    = ct*B012 + t*B123,

      C01  = ct*C0   + t*C1,  C12  = ct*C1  + t*C2, C23 = ct*C2 + t*C3,
      C012 = ct*C01  + t*C12, C123 = ct*C12 + t*C23,
      C    = ct*C012 + t*C123,

      D01  = ct*D0   + t*D1,  D12  = ct*D1  + t*D2, D23 = ct*D2 + t*D3,
      D012 = ct*D01  + t*D12, D123 = ct*D12 + t*D23,
      D    = ct*D012 + t*D123
    DO
      IF (A > 0.0D0 AND B < 0.0D0 AND C > 0.0D0 AND D < 0.0D0 AND P0 > 0.0D0)
         OR
         (A < 0.0D0 AND B > 0.0D0 AND C < 0.0D0 AND D > 0.0D0 AND P0 < 0.0D0)
      THEN
        RETURN t
      ELSE
        RETURN NoRoot
      END
    END
  END;
  
  WITH
    L1 = (P0 + P1)/2.0D0,
    M  = (P1 + P2)/2.0D0,
    L2 = (L1 +  M)/2.0D0,
    R2 = (P2 + P3)/2.0D0,
    R1 = (M  + R2)/2.0D0,
    L3 = (L2 + R1)/2.0D0,

    AL1 = ( A0 +  A1)/2.0D0,
    AM  = ( A1 +  A2)/2.0D0,
    AL2 = (AL1 +  AM)/2.0D0,
    AR2 = ( A2 +  A3)/2.0D0,
    AR1 = ( AM + AR2)/2.0D0,
    AL3 = (AL2 + AR1)/2.0D0,

    BL1 = ( B0 +  B1)/2.0D0,
    BM  = ( B1 +  B2)/2.0D0,
    BL2 = (BL1 +  BM)/2.0D0,
    BR2 = ( B2 +  B3)/2.0D0,
    BR1 = ( BM + BR2)/2.0D0,
    BL3 = (BL2 + BR1)/2.0D0,

    CL1 = ( C0 +  C1)/2.0D0,
    CM  = ( C1 +  C2)/2.0D0,
    CL2 = (CL1 +  CM)/2.0D0,
    CR2 = ( C2 +  C3)/2.0D0,
    CR1 = ( CM + CR2)/2.0D0,
    CL3 = (CL2 + CR1)/2.0D0,

    DL1 = ( D0 +  D1)/2.0D0,
    DM  = ( D1 +  D2)/2.0D0,
    DL2 = (DL1 +  DM)/2.0D0,
    DR2 = ( D2 +  D3)/2.0D0,
    DR1 = ( DM + DR2)/2.0D0,
    DL3 = (DL2 + DR1)/2.0D0
  DO
    t := FindRootB(P0,  L1,  L2,  L3, A0, AL1, AL2, AL3, B0, BL1, BL2, BL3,
                   C0, CL1, CL2, CL3, D0, DL1, DL2, DL3);
    IF t <= 1.0D0 THEN
      RETURN t/2.0D0
    ELSE
      t := FindRootB( L3,  R1,  R2, P3, AL3, AR1, AR2, A3, BL3, BR1, BR2, B3,
                     CL3, CR1, CR2, C3, DL3, DR1, DR2, D3);
      IF t <= 1.0D0 THEN RETURN 0.5D0 + t/2.0D0 ELSE RETURN NoRoot END
    END
  END
END FindRootB;

PROCEDURE FindRoot(P0, P1, P2, P3: LONGREAL; sign: BOOLEAN): LONGREAL =
CONST EPSILON = 1.0D-7;
VAR t: LONGREAL;
BEGIN
  IF (P0 >= 0.0D0 AND P1 >= 0.0D0 AND P2 >= 0.0D0 AND P3 >= 0.0D0) OR
     (P0 <= 0.0D0 AND P1 <= 0.0D0 AND P2 <= 0.0D0 AND P3 <= 0.0D0)
  THEN
    RETURN NoRoot
  END;
  IF ABS(P0 + P2 - P1 - P1) < EPSILON AND 
     ABS(P1 + P3 - P2 - P2) < EPSILON
  THEN
    IF (P0 > 0.0D0 AND P3 < 0.0D0 AND NOT sign) OR
       (P0 < 0.0D0 AND P3 > 0.0D0 AND sign)
    THEN
      WITH DP = P0 - P3 DO
        <* ASSERT P0 <= 0.0D0 OR P0 >= 0.0D0 *>
        <* ASSERT P3 <= 0.0D0 OR P3 >= 0.0D0 *>
        <* ASSERT DP <= 0.0D0 OR DP >= 0.0D0 *>
        <* ASSERT P0/DP >= 0.0D0 *>
        IF DP = 0.0D0 THEN RETURN 0.5D0 ELSE RETURN P0/DP END
      END
    ELSE
      RETURN NoRoot
    END
  END;
  WITH
    L1 = (P0 + P1)/2.0D0,
    M  = (P1 + P2)/2.0D0,
    L2 = (L1 +  M)/2.0D0,
    R2 = (P2 + P3)/2.0D0,
    R1 = (M  + R2)/2.0D0,
    L3 = (L2 + R1)/2.0D0
  DO
    t := FindRoot(P0, L1, L2, L3, sign);
    <* ASSERT t >= 0.0D0 *>
    IF t <= 1.0D0 THEN
      RETURN t/2.0D0
    ELSE
      t := FindRoot(L3, R1, R2, P3, sign);
      <* ASSERT t >= 0.0D0 *>
      IF t <= 1.0D0 THEN RETURN 0.5D0 + t/2.0D0 ELSE RETURN NoRoot END
    END
  END
END FindRoot;

BEGIN END BezierSearch.
