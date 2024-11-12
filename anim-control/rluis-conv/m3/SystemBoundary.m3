MODULE SystemBoundary;

IMPORT Word, Util, Fmt;
IMPORT SystemTopology AS ST;
FROM SystemTopology IMPORT
  NodeNum, NoNode,
  CellNum, NoCell;

PROCEDURE ComputeBoundary(t: ST.T): T =
  BEGIN
    WITH b = NEW(T) DO
      FindExposedFaces(t, b);
      FindExposedEdges(t, b);
      FindExposedSites(t, b);
      RETURN b
    END;
  END ComputeBoundary;

PROCEDURE FindExposedFaces(t: ST.T; b: T) =
  (*
    Identifies the exposed faces in "t".
    An exposed face is a triangle that occurs in only one cell.

    Specifically, computes "b.nFaces" and builds the "b.face" table (ecept for
    the "edge" fields).  For the time being, assumes that nodes and sites have
    the same numbering.
  *)
  CONST TwoCells = NoCell-1;
  BEGIN
    b.nFaces := 0;
    WITH
      maxFaces = 4 * t.nCells,
      hashSize = 2 * maxFaces,
      hash = NEW(REF ARRAY OF Face, hashSize)^
    DO
      FOR h := 0 TO LAST(hash) DO hash[h].cell := NoCell END;

      PROCEDURE MarkFace(cn: CellNum; pn, qn, rn: NodeNum) =
      VAR h: CARDINAL;
      BEGIN
        (* Cycle vertices until "pn" is minimum: *)
        WHILE pn > qn OR pn > rn DO 
          VAR tn: NodeNum := pn; BEGIN pn := qn; qn := rn; rn := tn END
        END;
        <* ASSERT pn < qn AND pn < rn *>
        WITH
          (* Hash nodes so that "(pn,qn,rn)" and "(pn,rn,qn)" collide: *)
          ha = Word.Times(417,    pn + 1),
          hb = Word.Times(4615,   qn + rn + 1),
          hc = Word.Times(480317, ABS(rn - qn))
        DO
          h := Word.Mod(Word.Plus(Word.Plus(ha, hb), hc), hashSize);
          (* Look for this triangle, until the first empty slot: *)
          WHILE hash[h].cell # NoCell 
          AND NOT (hash[h].site[0] = pn AND hash[h].site[0] = qn AND hash[h].site[0] = rn)
          AND NOT (hash[h].site[0] = pn AND hash[h].site[0] = rn AND hash[h].site[0] = qn)
          DO
            h := (h + 1) MOD hashSize
          END;
          WITH f = hash[h] DO
            IF f.cell = NoCell THEN
              (* First occurrence of this triangle *)
              f := Face{
                site := ARRAY [0..2] OF SiteNum{pn, qn, rn},
                cell := cn, 
                edge := ARRAY [0..2] OF EdgeNum{NoEdge, NoEdge, NoEdge}
              };
              INC(b.nFaces)
            ELSIF f.cell # TwoCells THEN
              (* This triangle occurred before: *)
              f.cell := TwoCells;
              DEC(b.nFaces);
              (* The previous occurrence must have opposite orientation: *)
              <* ASSERT f.site[1] = rn *>
              <* ASSERT f.site[2] = qn *>
            ELSE 
              (* Oops! three cells sharing the same face *)
              Util.Error("triply-shared triangle " & FmtZ3(pn, qn, rn))
            END
          END
        END
      END MarkFace;

      BEGIN
        FOR cn := 0 TO t.nCells-1 DO
          WITH c = t.cell[cn] DO
            MarkFace(cn, c.node[1], c.node[2], c.node[3]);
            MarkFace(cn, c.node[2], c.node[1], c.node[0]);
            MarkFace(cn, c.node[0], c.node[3], c.node[2]);
            MarkFace(cn, c.node[3], c.node[1], c.node[0]);
          END
        END;
      END;

      b.face := NEW(REF ARRAY OF Face, b.nFaces);
      nFaces := 0;
      VAR fn: FaceNum := 0;
      BEGIN
        FOR h := 0 TO LAST(hash) DO
          WITH f = hash[h] DO
            IF f.cell # NoCell AND f.cell # TwoCell THEN
              b.face[fn] := f; INC(fn); 
            END;
          END;
        END;
        <* ASSERT fn = b.nfaces *>
      END
    END
  END FindExposedFaces;

PROCEDURE FindExposedEdges(t: ST.T; b: T) =
  (* 
    Identifies the exposed edges of "t".
    An edge is exposed if it belongs to one or more exposed faces.

    Specifically, computes "b.nEdges", builds the "b.edge" table, 
    and fills the "edge" field of "b.face[0..b.nFaces]".
    Assumes "b.face" has been built, except for the "edge" fields
    in each face. Assumes also that sites and nodes have the same 
    numbering.
  *)
  BEGIN
    b.nEdges := 0;
    WITH
      maxEdges = 3 * b.nFaces,
      hashSize = 2 * maxEdges,
      hash = NEW(REF ARRAY OF Edge, hashSize)^,
      faces = b.faces^
    DO
      FOR h := 0 TO LAST(hash) DO hash[h].t := NoCell END;

      PROCEDURE MarkEdge(fn: FaceNum; pn, qn: NodeNum) =
      VAR h: CARDINAL;
      BEGIN
        (* Swap vertices until "pn" is minimum: *)
        IF pn > qn THEN 
          VAR tn: NodeNum := pn; BEGIN pn := qn; qn := tn END
        END;
        <* ASSERT pn < qn *>
        WITH
          ph = Word.Times(417,    pn + 1),
          qh = Word.Times(480317, qn + 1)
        DO
          h := Word.Mod(Word.Plus(ph, qh), hashSize);
          (* Look for this edge, until the first empty slot: *)
          WHILE hash[h].face # NoFace 
          AND NOT (hash[h].site[0] = pn AND hash[h].site[1] = qn) DO
            h := (h + 1) MOD hashSize
          END;
          WITH e = hash[h] DO
            IF e.face = NoFace THEN
              (* First occurrence of this edge *)
              WITH cn = b.face[fn].cell DO
                e := Edge{
                  site := ARRAY [0..1] OF Site{pn, qn},
                  face := fn, 
                  cell := cn,
                  dual := GetDualNodes(t.cell[cn], pn, qn)
                };
              END;
              INC(b.nEdges)
            END
          END
        END
      END MarkEdge;

      BEGIN
        FOR fn := 0 TO b.nFaces DO
          WITH f = b.face[fn] DO
            MarkEdge(fn, f.node[0], f.node[1]);
            MarkEdge(fn, f.node[1], f.node[2]);
            MarkEdge(fn, f.node[2], f.node[0]);
          END
        END;
      END;

      b.edge := NEW(REF ARRAY OF Edge, b.nEdges);
      VAR en: EdgeNum := 0;
      BEGIN
        FOR h := 0 TO LAST(hash) DO
          WITH e = hash[h] DO
            IF e.face # NoFace THEN
              b.edge[en] := e; 
              (* Set the edge pointers in the face record: *)
              WITH f = b.face[e.face] DO
                FOR r := 0 TO 2 DO
                  WITH s = (r+1) MOD 3 DO
                    IF (f.site[r] = e.site[0] AND f.site[s] = e.site[1])
                    OR (f.site[r] = e.site[1] AND f.site[s] = e.site[0])
                    THEN
                      f.edge[(r+2) MOD 3] := en
                    END
                  END
                END
              END;
              INC(en)
            END;
          END;
        END;
        <* ASSERT en = b.nEdges *>
      END
    END
  END FindExposedEdges;

PROCEDURE FindExposedSites(t: ST.T; b: T) =
  (*
    Finds the exposed sites of "t", builds the "b.site" table, and
    renumbers all "SiteNum"s in "b.face[fn].site", "b.edge[en].site"
  *)
  BEGIN
    b.nSites := 0;
    WITH
      siteFomNode = NEW(REF ARRAY OF SiteNum, t.nNodes)^
    DO
      (* Find exposed sites, build "siteFomNode" table, renumber sites in "b.face": *)
      FOR nn := 0 TO t.nNodes-1 DO siteFomNode[nn] := NoSite END;
      FOR fn := 0 TO b.nFaces-1 DO
        WITH f = b.face[fn] DO 
          FOR k := 0 TO 2 DO
            WITH nn = f.site[k] DO
              IF siteFomNode[nn] = NoSite THEN
                (* New exposed site *)
                siteFomNode[nn] := b.nSites; INC(b.nSites)
              END;
              f.site[k] := siteFomNode[nn]
            END
          END
        END
      END;
      (* Renumber sites in "b.edge": *)
      FOR en := 0 TO b.nEdges-1 DO 
        WITH e = b.edge[en] DO 
          FOR k := 0 TO 1 DO
            WITH nn = e.site[k] DO
              <* ASSERT siteFomNode[nn] # NoSite *>
              e.site[k] := siteFomNode[nn]
            END
          END
        END
      END;

      (* Build site table: *)
      b.site := NEW(REF ARRAY OF Site, b.nSites);
      VAR sn: SiteNum := 0;
      BEGIN
        FOR nn := 0 TO t.nNodes-1 DO 
          IF siteFomNode[n] # NoSite THEN 
            b.site[sn].node := nn; INC(sn)
          END
        END;
        <* ASSERT sn = b.nSites *>
      END
    END
  END FindExposedSites;

PROCEDURE FmtZ3(p, q, r: INTEGER): TEXT =
  BEGIN
    RETURN "(" & Fmt.Int(p) & ", " & Fmt.Int(q) & ", " & Fmt.Int(r) & ")"
  END FmtZ3;
  
PROCEDURE GetDualNodes(READONLY c: Cell; pn, qn: NodeNum): ARRAY [0..1] OF NodeNum =
  (*
    Given a cell "c" and two of its corners "pn, qn", returns the 
    other two corners "rn, sn", such that "(pn,qn,rn,sn)" is a 
    positive screw (assuming "c.node[0..3]" is already a positive screw).
  *)
  VAR rn, sn: NodeNum;
      xn: ARRAY [0..3] OF NodeNum := c.node;
  BEGIN
     IF pn = xn[1] THEN
       VAR tn := xn[0]; BEGIN xn[0] := xn[1]; xn[1] := xn[2]; xn[2] := tn END;
     END;
     IF pn = xn[2] THEN
       VAR tn := xn[0]; BEGIN xn[0] := xn[2]; xn[2] := xn[1]; xn[1] := tn END;
     END;
     IF pn = xn[3] THEN
       VAR tn := xn[0]; BEGIN xn[0] := xn[3]; xn[3] := xn[1]; xn[1] := tn END;
     END;
     IF qn = xn[2] THEN
       VAR tn := xn[1]; BEGIN xn[1] := xn[2]; xn[2] := xn[3]; xn[3] := tn END;
     END;
     IF qn = xn[3] THEN
       VAR tn := xn[1]; BEGIN xn[1] := xn[3]; xn[3] := xn[2]; xn[2] := tn END;
     END;
     <* ASSERT pn = x[0] AND qn = x[1] *>
     RETURN ARRAY [0..1] OF NodeNum{x[1], x[2]}
    END
  END GetDualNodes;

BEGIN
END SystemBoundary.





