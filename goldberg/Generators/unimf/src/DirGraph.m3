MODULE DirGraph;

(* Created April 1997 by J. Stolfi  *)
(* See the copyright notice at the end of this file. *)

IMPORT Wr, Fmt, Thread, XRandom, FileFmt;

PROCEDURE Write(
    wr: Wr.T;
    READONLY g: T;
    vBase: INTEGER;
    READONLY omit: SET OF Mark := NoMark;
  ) =
  <* FATAL Thread.Alerted, Wr.Failure *>
  VAR NE: Count;
  BEGIN
    WITH
      NV = g.NV,
      e = g.e^,
      m = g.m^,
      D = Digits(vBase, NV - 1 + vBase)
    DO
      (* Count output edges: *)
      IF omit = SET OF Mark{} THEN 
        NE := g.NE
      ELSE
        (* Count non-omitted edges: *)
        NE := 0;
        FOR i := 0 TO LAST(e) DO
          IF NOT (m[i] IN omit) THEN INC(NE) END
        END
      END;

      FileFmt.WriteComment(wr, g.cmt, 'c');
     
      Wr.PutText(wr, "p ");
      Wr.PutText(wr, g.class);
      Wr.PutText(wr, " " & Fmt.Int(NV));
      Wr.PutText(wr, " " & Fmt.Int(NE));
      Wr.PutText(wr, "\n");
      
      Wr.PutText(wr, "n");
      Wr.PutText(wr, " ");
      Wr.PutText(wr, Fmt.Pad(Fmt.Int(g.s + vBase), D));
      Wr.PutText(wr, " ");
      Wr.PutText(wr, "s");
      Wr.PutText(wr, "\n");
      
      Wr.PutText(wr, "n");
      Wr.PutText(wr, " ");
      Wr.PutText(wr, Fmt.Pad(Fmt.Int(g.t + vBase), D));
      Wr.PutText(wr, " ");
      Wr.PutText(wr, "t");
      Wr.PutText(wr, "\n");
      
      FOR i := 0 TO LAST(e) DO 
        IF NOT m[i] IN omit THEN 
          Wr.PutText(wr, "a");
          Wr.PutText(wr, " ");
          Wr.PutText(wr, Fmt.Pad(Fmt.Int(e[i][0] + vBase), D));
          Wr.PutText(wr, " ");
          Wr.PutText(wr, Fmt.Pad(Fmt.Int(e[i][1] + vBase), D)); 
          IF g.c # NIL THEN
            Wr.PutText(wr, " ");
            Wr.PutText(wr, Fmt.Int(g.c[i]))
          END;
          Wr.PutText(wr, "\n");
        END;
      END;
      Wr.Flush(wr);
    END
  END Write;

PROCEDURE Digits(lo, hi: INTEGER): CARDINAL =
  (* Mas number of bytes needed for an integer in "[lo..hi]" *)
  VAR m, d: CARDINAL := 1;
  BEGIN
    IF lo < 0 THEN m := MAX(ABS(hi), ABS(lo)) ELSE m := ABS(hi) END;
    WHILE m > 10 DO m := m DIV 10; INC(d) END;
    IF lo < 0 THEN RETURN d+1 ELSE RETURN d END;
  END Digits;

PROCEDURE PermuteVertexNums(VAR g: T; rnd: XRandom.T) =
  BEGIN
    WITH
      NV = g.NV,
      NE = g.NE,
      e = g.e^,
      map = NEW(REF ARRAY OF CARDINAL, NV)^
    DO
      FOR v := 0 TO NV-1 DO map[v] := v END; Scramble(rnd, map);
      (* Translate edges: *)
      FOR i := 0 TO NE-1 DO 
        WITH ei = e[i] DO
          ei[0] := map[ei[0]]; ei[1] := map[ei[1]]
        END;
      END;
      (* Translate source and sink: *)
      IF g.s < NV THEN g.s := map[g.s] END;
      IF g.t < NV THEN g.t := map[g.t] END;
    END;
  END PermuteVertexNums;
  
PROCEDURE Scramble(rnd: XRandom.T; VAR v: ARRAY OF CARDINAL) =
  (* Permutes "v" in random order *)
  VAR t: CARDINAL;
  BEGIN
    FOR k := LAST(v) TO 1 BY -1 DO 
      WITH i = rnd.integer(0, k) DO 
        IF i # k THEN
          t := v[i]; v[i] := v[k]; v[k] := t
        END
      END
    END
  END Scramble;
  
PROCEDURE ComplementVertexNums(VAR g: T) =
  BEGIN
    WITH
      NV = g.NV,
      e = g.e^
    DO
      (* Translate edges: *)
      FOR i := 0 TO LAST(e) DO 
        WITH ei = e[i] DO
          WITH x = ei[0] DO x := NV - 1 - x END; 
          WITH x = ei[1] DO x := NV - 1 - x END
        END;
      END;
      (* Translate source and sink: *)
      IF g.s < NV THEN g.s := NV - 1 - g.s END;
      IF g.t < NV THEN g.t := NV - 1 - g.t END;
    END;
  END ComplementVertexNums;
  
PROCEDURE ReverseEdges(VAR e: Edges) =
  BEGIN
    FOR i := 0 TO LAST(e) DO 
      WITH uv = e[i], u = uv[0], v = uv[1] DO
        e[i] := Edge{v, u}
      END
    END
  END ReverseEdges;

BEGIN
END DirGraph.
