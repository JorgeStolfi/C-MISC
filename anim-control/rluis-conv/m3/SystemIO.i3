INTERFACE SystemIO;

IMPORT Rd, Wr, LR3, LR3x3;

TYPE
  Tetrahedron = RECORD
      p0, p1, p2, p3: CARDINAL;
      A: LR3x3.T;
      density, alpha, beta, eta1, eta2: LONGREAL;
    END;
  Face = RECORD u, v, w, tx: CARDINAL; mu1, mu2, e: LONGREAL END;
  Constraint = RECORD a: CARDINAL END;
  Vectors3D = REF ARRAY OF LR3.T;
  PlainVectors3D = REF ARRAY OF LONGREAL;
  
PROCEDURE ReadVectors(rd: Rd.T; n: CARDINAL; pos, vel: Vectors3D);
PROCEDURE ReadPlainVectors(rd: Rd.T; n: CARDINAL; pos, vel: PlainVectors3D);
PROCEDURE ReadFace(rd: Rd.T; VAR f: Face);
PROCEDURE ReadFixed(rd: Rd.T; VAR i: CARDINAL);
PROCEDURE ReadKinetic(rd: Rd.T; 
  VAR ta, tb: LONGREAL; VAR k: CARDINAL; VAR p, v, xi: LONGREAL;
);
PROCEDURE ReadTexture(rd: Rd.T; VAR p: Texture);
PROCEDURE WriteVectors(wr: Wr.T; n: CARDINAL; pos, vel: Vectors3D; 
                       base: INTEGER);
PROCEDURE WritePlainVectors(wr: Wr.T; n: CARDINAL; pos, vel: PlainVectors3D;
                            base: INTEGER);
PROCEDURE WriteFace(wr: Wr.T; READONLY f: Face; base: CARDINAL);
PROCEDURE WriteFixed(wr: Wr.T; i: CARDINAL; base: CARDINAL);
PROCEDURE WriteKinetic(wr: Wr.T; ta, tb: LONGREAL; k: CARDINAL; p, v, xi: LONGREAL);
PROCEDURE WriteTexture(wr: Wr.T; READONLY p: Texture);

END SystemIO.
