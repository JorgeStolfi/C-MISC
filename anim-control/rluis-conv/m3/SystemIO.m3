MODULE SystemIO;

IMPORT Rd, Wr, Lex, FloatMode, Text, Util;
FROM Thread IMPORT Alerted;
FROM Util IMPORT WriteN, WriteR, WriteLR, WriteSpace, WriteEOL;

CONST (* Original: *)
  TP_VERSION = "95-06-03";
  ST_VERSION = "95-06-03";
  SU_VERSION = "95-09-14";
  CO_VERSION = "95-07-10";
  TX_VERSION = "95-09-14";
  EN_VERSION = "95-07-06";
  IM_VERSION = "95-07-05";

CONST
  ST_VERSION = "95-06-03";
  SU_VERSION = "95-09-14";
  FX_VERSION = "96-12-04";
  KC_VERSION = "96-12-04";
  TX_VERSION = "95-09-14";
  EN_VERSION = "95-07-06";
  IM_VERSION = "95-07-05";
  MM_VERSION = "96-12-04";

(*--- READING AND WRITING OF HEADERS AND FOOTERS ----------------------------*)

PROCEDURE GetVersion(type: TEXT): TEXT =
BEGIN
  IF    Text.Equal(type, "topology")    THEN RETURN TP_VERSION
  ELSIF Text.Equal(type, "state")       THEN RETURN ST_VERSION
  ELSIF Text.Equal(type, "surface")     THEN RETURN SU_VERSION
  ELSIF Text.Equal(type, "fixed")       THEN RETURN FX_VERSION
  ELSIF Text.Equal(type, "kinetics")    THEN RETURN KC_VERSION
  ELSIF Text.Equal(type, "textures")    THEN RETURN TX_VERSION
  ELSIF Text.Equal(type, "energies")    THEN RETURN EN_VERSION
  ELSIF Text.Equal(type, "impulses")    THEN RETURN IM_VERSION
  ELSIF Text.Equal(type, "matrix")      THEN RETURN MM_VERSION
  ELSE <* ASSERT FALSE *>
  END
END GetVersion;

(*--- READING OF TYPES ------------------------------------------------------*)

PROCEDURE ReadVectors(rd: Rd.T; n: CARDINAL; pos, vel: Vectors3D) =
<* FATAL Rd.EndOfFile, Rd.Failure, Alerted, FloatMode.Trap, Lex.Error *>
BEGIN
  FOR i := 0 TO n-1 DO
    EVAL Lex.Int(rd); EVAL Rd.GetChar(rd);
    pos[i][0] := Lex.LongReal(rd); Lex.Skip(rd);
    pos[i][1] := Lex.LongReal(rd); Lex.Skip(rd);
    pos[i][2] := Lex.LongReal(rd); Lex.Skip(rd);
    vel[i][0] := Lex.LongReal(rd); Lex.Skip(rd);
    vel[i][1] := Lex.LongReal(rd); Lex.Skip(rd);
    vel[i][2] := Lex.LongReal(rd); Lex.Skip(rd);
  END
END ReadVectors;

PROCEDURE ReadPlainVectors(rd: Rd.T; n: CARDINAL; pos, vel: PlainVectors3D) =
<* FATAL Rd.EndOfFile, Rd.Failure, Alerted, FloatMode.Trap, Lex.Error *>
BEGIN
  FOR i := 0 TO n-1 BY 3 DO
    EVAL Lex.Int(rd); EVAL Rd.GetChar(rd);
    pos[i+0] := Lex.LongReal(rd); Lex.Skip(rd);
    pos[i+1] := Lex.LongReal(rd); Lex.Skip(rd);
    pos[i+2] := Lex.LongReal(rd); Lex.Skip(rd);
    vel[i+0] := Lex.LongReal(rd); Lex.Skip(rd);
    vel[i+1] := Lex.LongReal(rd); Lex.Skip(rd);
    vel[i+2] := Lex.LongReal(rd); Lex.Skip(rd);
  END
END ReadPlainVectors;

PROCEDURE ReadFace(rd: Rd.T; VAR f: Face) =
<* FATAL Rd.Failure, Alerted, FloatMode.Trap, Lex.Error *>
BEGIN
  f.u   := Lex.Int(rd);      Lex.Skip(rd);
  f.v   := Lex.Int(rd);      Lex.Skip(rd);
  f.w   := Lex.Int(rd);      Lex.Skip(rd);
  f.tx  := Lex.Int(rd);      Lex.Skip(rd);
  f.mu1 := Lex.LongReal(rd); Lex.Skip(rd);
  f.mu1 := Lex.LongReal(rd); Lex.Skip(rd);
  f.e   := Lex.LongReal(rd); Lex.Skip(rd);
  f.e   := (1.0D0 + f.e)/2.0D0;
END ReadFace;

PROCEDURE ReadFixed(rd: Rd.T; VAR i: CARDINAL) =
<* FATAL Rd.Failure, Alerted, FloatMode.Trap, Lex.Error *>
BEGIN
  i := Lex.Int(rd); Lex.Skip(rd);
END ReadFixed;

PROCEDURE ReadKinetic(rd: Rd.T; 
  VAR ta, tb: LONGREAL; VAR k: CARDINAL; VAR p, v, xi: LONGREAL;
) =
VAR i: CARDINAL;
VAR c: CHAR;
<* FATAL FloatMode.Trap, Rd.Failure, Rd.EndOfFile, Alerted, Lex.Error *>
BEGIN
  ta := Lex.LongReal(rd); Lex.Skip(rd);
  tb := Lex.LongReal(rd); Lex.Skip(rd);
  i  := Lex.Int(rd);      Lex.Skip(rd);
  c  := Rd.GetChar(rd);   Lex.Skip(rd);
  p  := Lex.LongReal(rd); Lex.Skip(rd);
  v  := Lex.LongReal(rd); Lex.Skip(rd);
  xi := Lex.LongReal(rd); Lex.Skip(rd);
  IF    c = 'x' THEN k := 3*i+0
  ELSIF c = 'y' THEN k := 3*i+1
  ELSIF c = 'z' THEN k := 3*i+2
  ELSE RAISE Lex.Error
  END
END ReadKinetic;

PROCEDURE ReadTexture(rd: Rd.T; VAR m: Texture) =
<* FATAL Rd.Failure, Alerted, FloatMode.Trap, Lex.Error *>
BEGIN
  m.aR := Lex.Real(rd); Lex.Skip(rd);
  m.aG := Lex.Real(rd); Lex.Skip(rd);
  m.aB := Lex.Real(rd); Lex.Skip(rd);
  m.dR := Lex.Real(rd); Lex.Skip(rd);
  m.dG := Lex.Real(rd); Lex.Skip(rd);
  m.dB := Lex.Real(rd); Lex.Skip(rd);
  m.sR := Lex.Real(rd); Lex.Skip(rd);
  m.sG := Lex.Real(rd); Lex.Skip(rd);
  m.sB := Lex.Real(rd); Lex.Skip(rd);
  m.tR := Lex.Real(rd); Lex.Skip(rd);
  m.tG := Lex.Real(rd); Lex.Skip(rd);
  m.tB := Lex.Real(rd); Lex.Skip(rd);
  m.ir := Lex.Real(rd); Lex.Skip(rd);
  m.n  := Lex.Real(rd); Lex.Skip(rd);
END ReadTexture;


(*--- WRITING OF TYPES ------------------------------------------------------*)

PROCEDURE WriteTetrahedron(wr: Wr.T; READONLY t: Tetrahedron; 
                           nbase: CARDINAL) =
BEGIN
END WriteTetrahedron;
  
PROCEDURE WriteVectors(wr: Wr.T; n: CARDINAL; pos, vel: Vectors3D; 
                       base: INTEGER) =
<* FATAL Wr.Failure, Alerted *>
BEGIN
  WITH d = Util.Digits(FLOAT(n, LONGREAL)) DO
    FOR i := 0 TO n-1 DO
      WriteN(wr, base + i, d); Wr.PutText(wr, ": ");
      WriteLR(wr, pos[i][0]); WriteSpace(wr); 
      WriteLR(wr, pos[i][1]); WriteSpace(wr);
      WriteLR(wr, pos[i][2]); WriteSpace(wr, 2);
      WriteLR(wr, vel[i][0]); WriteSpace(wr);
      WriteLR(wr, vel[i][1]); WriteSpace(wr);
      WriteLR(wr, vel[i][2]); WriteEOL(wr);
    END
  END
END WriteVectors;

PROCEDURE WritePlainVectors(wr: Wr.T; n: CARDINAL; pos, vel: PlainVectors3D;
                            base: INTEGER) =
<* FATAL Wr.Failure, Alerted *>
BEGIN
  WITH d = Util.Digits(FLOAT(n DIV 3, LONGREAL)) DO
    FOR i := 0 TO n-1 BY 3 DO
      WriteN(wr, base + i DIV 3, d); Wr.PutText(wr, ": ");
      WriteLR(wr, pos[i+0]); WriteSpace(wr); 
      WriteLR(wr, pos[i+1]); WriteSpace(wr);
      WriteLR(wr, pos[i+2]); WriteSpace(wr, 2);
      WriteLR(wr, vel[i+0]); WriteSpace(wr);
      WriteLR(wr, vel[i+1]); WriteSpace(wr);
      WriteLR(wr, vel[i+2]); WriteEOL(wr);
    END
  END
END WritePlainVectors;

PROCEDURE WriteFace(wr: Wr.T; READONLY f: Face; base: CARDINAL) =
BEGIN
  WriteN(wr, base + f.u); WriteSpace(wr);
  WriteN(wr, base + f.v); WriteSpace(wr);
  WriteN(wr, base + f.w); WriteSpace(wr);
  WriteN(wr, f.tx);       WriteSpace(wr);
  WriteLR(wr, f.mu1);     WriteSpace(wr);
  WriteLR(wr, f.mu2);     WriteSpace(wr);
  WriteLR(wr, f.e);       WriteEOL(wr);
END WriteFace;

PROCEDURE WriteFixed(wr: Wr.T; i: CARDINAL; base: CARDINAL) =
BEGIN
  WriteN(wr, base + i); WriteEOL(wr);
END WriteFixed;

PROCEDURE WriteKinetic(wr: Wr.T; ta, tb: LONGREAL; k: CARDINAL; p, v, xi: LONGREAL) =
<* FATAL Wr.Failure, Alerted *>
BEGIN
  WriteLR(wr, ta); WriteSpace(wr);
  WriteLR(wr, tb); WriteSpace(wr, 2);
  WriteN (wr, k DIV 3); WriteSpace(wr);
  Wr.PutChar(wr, VAL(ORD('x') + (k MOD 3), CHAR)); WriteSpace(wr, 2);
  WriteLR(wr, p);  WriteSpace(wr);
  WriteLR(wr, v);  WriteSpace(wr);
  WriteLR(wr, xi); WriteSpace(wr);
END WriteKinetic;

PROCEDURE WriteTexture(wr: Wr.T; READONLY m: Texture) =
BEGIN
  WriteR(wr, m.aR); WriteSpace(wr);
  WriteR(wr, m.aG); WriteSpace(wr);
  WriteR(wr, m.aB); WriteSpace(wr, 2);
  WriteR(wr, m.dR); WriteSpace(wr);
  WriteR(wr, m.dG); WriteSpace(wr);
  WriteR(wr, m.dB); WriteSpace(wr, 2);
  WriteR(wr, m.sR); WriteSpace(wr);
  WriteR(wr, m.sG); WriteSpace(wr);
  WriteR(wr, m.sB); WriteSpace(wr, 2);
  WriteR(wr, m.tR); WriteSpace(wr);
  WriteR(wr, m.tG); WriteSpace(wr);
  WriteR(wr, m.tB); WriteSpace(wr, 2);
  WriteR(wr, m.ir); WriteSpace(wr, 2);
  WriteR(wr, m.n);  WriteEOL(wr);
END WriteTexture;

BEGIN END SystemIO.
