MODULE SparseSymetricMatrix;

(* Written by R.L.W.Liesenfeld, 1995 *)

IMPORT Rd, Wr, Fmt, Lex, Util, FileFmt, FPut, FGet, NPut, NGet;
FROM Thread IMPORT Alerted;

PROCEDURE New(n: CARDINAL): T =
  BEGIN
    WITH
      m = NEW(T, n)
    DO
      FOR i := 0 TO n-1 DO
        WITH mi = m[i] DO
          mi.d := 0;
          mi.row := NEW(REF Row, 1);
          mi.row[0] := MatrixElement{j := i, a := 1.0d0}
        END
      END
    END
  END New;

PROCEDURE GetValue(READONLY m: T; i, j: CARDINAL): LONGREAL =
  VAR k, n: CARDINAL;
  BEGIN
    WITH mi = m[i], row = mi.row^ DO
      IF j < i THEN
        k := 0;
        n := mi.d;
      ELSE
        k := mi.d;
        n := NUMBER(row);
      END;
      WHILE k < n AND row[k].j < j DO INC(k) END;
      IF k < n AND row[k].j = j THEN RETURN row[k].a ELSE RETURN 0.0D0 END
    END
  END GetValue;

PROCEDURE GetDiagonalValue(READONLY m: T; i: CARDINAL): LONGREAL =
  BEGIN 
    RETURN m[i].row[m[i].d].a
  END GetDiagonalValue;

PROCEDURE SubMul(READONLY m: T; READONLY v: Vector; VAR x: Vector) =
  BEGIN
    FOR i := 0 TO LAST(m) DO
      WITH r = m[i].row^ DO
        FOR k := 0 TO LAST(r) DO
          x[i] := x[i] - r[k].a * v[r[k].j]
        END
      END
    END
  END SubMul;

PROCEDURE PutRow(VAR m: T; i: CARDINAL; READONLY row: Vector) =
  VAR k := 0;
  BEGIN
    <* ASSERT row[i] # 0.0D0 *>
    FOR j := 0 TO LAST(row) DO IF row[j] # 0.0D0 THEN INC(k) END END;
    WITH crow = NEW(REF Row, k), cr = crow^ DO
      k := 0;
      FOR j := 0 TO LAST(row) DO
        IF row[j] # 0.0D0) OR j = i THEN
          cr[k].a := row[j];
          cr[k].j := j;
          IF j = i THEN m[i].d := k END;
          INC(k);
        END
      END;
      m[i].row := crow;
    END
  END PutRow;

PROCEDURE Factor(READONLY m: T): REF T =
  VAR a, s: LONGREAL;
      p: CARDINAL;
  BEGIN
    WITH
      n   = NUMBER(m),
      rrf = NEW(REF Vector, n), rf = rrf^,
      rmf  = NEW(REF T, n), mf = rmf^
    DO
      FOR i := 0 TO n-1 DO
        WITH 
          ri = m[i].row^,
          jd = m[i].d,
          nz = NUMBER(ri)
        DO
          FOR j := 0 TO i-1 DO rf[j] := GetValue(mf, j, i) END;

          s := 0.0D0;
          FOR j := 0 TO i-1 DO 
            WITH rfj = rf[j] DO s := s + GetDiagonalValue(mf, j)*rfj*rfj END
          END;
          rf[i] := ri[jd].a - s;

          p := jd;
          FOR j := i+1 TO n-1 DO
            WHILE p < nz AND ri[p].j < j DO INC(p) END;
            IF p < nz AND ri[p].j = j THEN a := ri[p].a ELSE a := 0.0D0 END;
            s := 0.0D0;
            FOR k := 0 TO i-1 DO
              s := s + GetDiagonalValue(mf, k)*rf[k]*GetValue(mf, k, j)
            END;
            rf[j] := (a - s)/rf[i];
          END;
          PutRow(mf, i, rf);
        END;
      END;
      RETURN rmf
    END
  END Factor;

PROCEDURE Solve(READONLY mf: T; READONLY v: Vector; VAR x: Vector; i0: INTEGER) =
  VAR s: LONGREAL;
  BEGIN
    WITH 
      n = NUMBER(mf) 
    DO
      FOR i := 0 TO i0-1 DO x[i] := 0.0D0 END;

      (* Solves "LD*y = v", with "y" stored into x *)
      FOR i := i0 TO n-1 DO
        WITH r = mf[i].row^, dj = mf[i].d DO
          s := 0.0D0;
          FOR j := 0 TO dj-1 DO s := s + r[j].a * x[r[j].j] END;
          x[i] := (v[i] - s)/r[dj].a;
        END
      END;

      (* Solves L'x = y *)
      FOR i := n-1 TO 0 BY -1 DO
        WITH r = mf[i].row^, dj = mf[i].d DO
          s := 0.0D0;
          FOR j := dj+1 TO LAST(r) DO s := s + r[j].a*x[r[j].j] END;
          x[i] := x[i] - s;
        END
      END
    END
  END Solve;

PROCEDURE NonZeros(READONLY m: T): CARDINAL =
  VAR nz := 0;
  BEGIN
    FOR i := 0 TO LAST(m) DO 
      WITH r = m[i].row^ DO
        nz := nz + NUMBER(r);
        IF r[m[i].d].a = 0.0d0 THEN DEC(nz) END
      END
    END;
    RETURN nz;
  END NonZeros;

PROCEDURE Print(wr: Wr.T; READONLY m: T) =
  <* FATAL Alerted, Wr.Failure *>
  BEGIN
    WITH
      n  = NUMBER(m),
      nr = FLOAT(n, LONGREAL),
      n2 = 100.0D0/(nr*nr),
      nz = NonZeros(m),
      nzr = FLOAT(nz, LONGREAL),
      d  = Util.Digits(n)
    DO
      Wr.PutText(wr, 
        "Order: " & Fmt.Int(n) & "\n" &
        "Non-zeros: " & Fmt.Int(nz) & " = " & Fmt.LongReal(nzr*n2, prec := 1) & "%\n"
      );
      FOR i := 0 TO n-1 DO
        Wr.PutText(wr, Fmt.Pad(Fmt.Int(i), d) & ":");
        FOR j := 0 TO m[i].d DO
          Wr.PutText(wr, 
            " " & Fmt.Pad(Fmt.Int(m[i].row[j].j), d) & "=" & Fmt.LongReal(m[i].row[j].a)
          )
        END;
        Wr.PutText(wr, "\n");
      END
    END;
    Wr.Flush(wr);
  END Print;

PROCEDURE Remove(VAR m: T; READONLY sel: ARRAY OF BOOLEAN) =
  VAR nk, newnz: CARDINAL;
  BEGIN
    <* ASSERT NUMBER(sel) = NUMBER(m) *>
    WITH
      RM = LAST(CARDINAL)
    DO
      FOR i := 0 TO LAST(m) DO
        WITH
          ri = m[i].row^,
          nz = NUMBER(ri)
        DO
          IF sel[i] THEN
            m[i].row := NEW(REF Row, 1);
            m[i].row[0].j := i;
            m[i].row[0].a := 1.0D0;
            m[i].d := 0;
          ELSE
            newnz := nz;
            FOR k := 0 TO nz-1 DO
              IF ri[k].j = i THEN
                ri[k].a := 1.0d0
              ELSIF sel[ri[k].j] THEN 
                ri[k].j := RM; DEC(newnz)
              END;
            END;
            WITH newrow = NEW(REF Row, newnz), nri = newrow^ DO
              nk := 0;
              FOR k := 0 TO nz-1 DO
                IF ri[k].j # RM THEN
                  nri[nk] := ri[k];
                  IF nri[nk].j = i THEN m[i].d := nk END;
                  INC(nk);
                END
              END;
              m[i].row := newrow;
            END
          END
        END
      END
    END
  END Remove;

CONST MatrixFileVersion = "96-12-04";

PROCEDURE Read(rd: Rd.T; fp: LONGREAL): REF T =
  <* FATAL Alerted, Rd.Failure *>
  BEGIN
    FileFmt.ReadHeader(rd, "matrix", MatrixFileVersion);
    WITH
      n = NGet.Int(rd, "order"),
      fr = NGet.LongReal(rd, "fingerprint"),
      matrix = NEW(REF T, n),
      m = matrix^
    DO
      <* ASSERT fp = 0.0d0 OR fp = fr *>
      FOR i := 0 TO n-1 DO
        Lex.Skip(rd);
        WITH ir = FGet.Int(rd) DO <* ASSERT ir = i *> END;
        FGet.Colon(rd);
        WITH nz = FGet.Int(rd) DO 
          FGet.EOL(rd);
          m[i].row := NEW(REF Row, nz);
          m[i].d := LAST(CARDINAL);
          WITH ri = m[i].row^ DO
            FOR j := 0 TO nz-1 DO
              Lex.Skip(rd);
              ri[j].j := FGet.Int(rd);
              ri[j].a := FGet.LongReal(rd);
              <* ASSERT ri[j].j # i OR ri[j].a # 0.0D0 *>
              IF ri[j].j = i THEN m[i].d := j END;
            END
          END;
          FGet.EOL(rd);
          <* ASSERT m[i].d # LAST(CARDINAL) *>
        END
      END;
      FileFmt.ReadFooter(rd, "matrix");
      RETURN matrix
    END;
  END Read;

PROCEDURE Write(wr: Wr.T; READONLY m: T; fp: LONGREAL) =
  <* FATAL Alerted, Wr.Failure *>
  BEGIN
    FileFmt.WriteHeader(wr, "matrix", MatrixFileVersion);
    NPut.Int(wr, "order", NUMBER(m));
    NPut.LongReal(wr, "fingerprint", fp);
    FOR i := 0 TO LAST(m) DO
      WITH ri = m[i].row^, nz = NUMBER(ri) DO
        FPut.Int(wr, nz); FPut.Colon(wr);
        FOR j := 0 TO nz-1 DO
          IF j MOD 8 = 0 THEN
            FPut.EOL(wr); FPut.Space(wr)
          END;
          FPut.Space(wr);
          FPut.Int(wr, ri[j].j);
          FPut.Space(wr);
          FPut.LongReal(wr, ri[j].a);
        END;
        FPut.EOL(wr);
      END;
    END;
    FileFmt.WriteFooter(wr, "matrix");
    Wr.Flush(wr);
  END Write;

BEGIN END SparseSymetricMatrix.

