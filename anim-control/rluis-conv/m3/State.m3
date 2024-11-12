MODULE State;

IMPORT Rd, Wr, Thread, Lex, FloatMode;

REVEAL
  T = Public BRANDED OBJECT
    OVERRIDES
      alloc := Alloc;
      read := Read;
      write := Write;
    END;

CONST 
  StateFileVersion = "96-12-26";

PROCEDURE Alloc(s: T; nNodes: CARDINAL): T =
  BEGIN
    WITH NC = 3*nNodes DO
      IF s.pos = NIL OR NUMBER(s.pos^) # NC THEN 
        s.pos := NEW(REF Position, NC)
      END;
      IF s.vel = NIL OR NUMBER(s.vel^) # NC THEN
        s.vel := NEW(REF Velocity, NC)
      END;
      RETURN s
    END;
  END Alloc;

PROCEDURE Read(s: T; rd: Rd.T): T =
  <* FATAL Rd.Failure, Rd.EndOfFile, Thread.Alerted, FloatMode.Trap, Lex.Error *>
  BEGIN
    FileFmt.ReadHeader(rd, "state", StateFileVersion);
    s.comment := FileFmt.ReadComment(rd, '|');
    WITH 
      nNodes = NGet.Int(rd, "nodes"),
      NC = 3*nNodes
    DO
      EVAL s.alloc(nNodes);
      s.time := NGet.LongReal(rd, "time");
      FOR i := 0 TO s.nNodes-1 DO
        WITH
          iCheck = Lex.Int(rd),
          kStart := 3*i
        DO
          Lex.Skip(rd);
          WITH iCheck = FGet.Int(rd) DO <* ASSERT iCheck = i *> END; Lex.Skip(rd);
          FGet.Colon(rd);
          FOR k := kStart TO kStart + 2 DO s.pos[k] := FGet.LongReal(rd) END;
          FOR k := kStart TO kStart + 2 DO s.vel[k] := FGet.LongReal(rd) END
        END
      END;
    END;
    FileFmt.ReadFooter(rd, "state");
    RETURN s
  END Read;
  
PROCEDURE Write(s: T; wr: Wr.T) =
  <* FATAL Wr.Failure, Thread.Alerted, FloatMode.Trap, Lex.Error *>
  VAR k: CARDINAL;
  BEGIN
    FileFmt.WriteHeader(wr, "state", StateFileVersion);
    FileFmt.WriteComment(wr, s.comment, '|');
    WITH
      pos = s.pos^,
      vel = s.vel^,
      NC = NUMBER(pos)
      nNodes = NC DIV 3
    DO
      NPut.Int(wr, "nodes", nNodes);
      NPut.LongReal(wr, "time", s.time);
      k := 0;
      FOR i := 0 TO s.nNodes-1 DO
        (* Node number: *)
        FPut.Int(wr, i); FPut.Colon(wr); 
        FPut.Space(wr);
        WITH kStart = 3*i DO
          FOR k := kStart TO kStart + 2 DO 
            FPut.Space(wr);
            FPut.LongReal(wr, pos[k]);
          END;
          FPut.Space(wr);
          FOR k := kStart TO kStart + 2 DO
            FPut.Space(wr);
            FPut.LongReal(wr, vel[k]);
          END
        END;
      END
    END;
    FileFmt.WriteFooter(wr, "state");
  END Write;

BEGIN
END State.
