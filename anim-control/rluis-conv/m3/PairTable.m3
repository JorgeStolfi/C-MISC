MODULE PairTable;

IMPORT SystemTopology;

PROCEDURE New(l, c: CARDINAL): T =
BEGIN
  WITH table = NEW(T, l) DO
    FOR i := 0 TO l-1 DO
      table[i].row := NEW(Row, c)
    END;
    RETURN table;
  END
END New;

PROCEDURE On(t: T; i: CARDINAL): BOOLEAN =
BEGIN RETURN t[i].on END On;

PROCEDURE SetOn(t: T; i: CARDINAL; on: BOOLEAN) =
BEGIN t[i].on := on END SetOn;

PROCEDURE Valid(t: T; i, j: CARDINAL): BOOLEAN =
BEGIN RETURN t[i].row[j].v END Valid;

PROCEDURE SetValid(t: T; i, j: CARDINAL; v: BOOLEAN) =
BEGIN t[i].row[j].v := v END SetValid;

PROCEDURE Overlap(t: T; i, j: CARDINAL): BOOLEAN =
BEGIN
  WITH flag = t[i].row[j].flag DO
    RETURN flag[Axis.X] AND flag[Axis.Y] AND flag[Axis.Z]
  END
END Overlap;

PROCEDURE BuildVertexFacePairs(t: T; s: SystemTopology.T): CARDINAL =
VAR n := 0;
BEGIN
  FOR i := 0 TO LAST(s.vertices^) DO
    WITH
      d = s.vertices[i].a DIV 3,
      vertexFixed = s.fixed[d]
    DO
      t[i].on := TRUE;
      FOR j := 0 TO LAST(s.faces^) DO
        WITH
          a = s.faces[j].a DIV 3,
          b = s.faces[j].b DIV 3,
          c = s.faces[j].c DIV 3,
          faceFixed = s.fixed[a] AND s.fixed[b] AND s.fixed[c],
          pair = t[i].row[j]
        DO
          pair.v := NOT (s.sameTetrahedron(a, b, c, d) OR
                         vertexFixed AND faceFixed);
          pair.flag[Axis.X] := FALSE;
          pair.flag[Axis.Y] := FALSE;
          pair.flag[Axis.Z] := FALSE;
          IF pair.v THEN INC(n) END
        END
      END
    END
  END;
  RETURN n
END BuildVertexFacePairs;

PROCEDURE BuildEdgeEdgePairs(t: T; s: SystemTopology.T): CARDINAL =
VAR n := 0;
BEGIN
  FOR i := 0 TO LAST(s.edges^) DO
    t[i].on := TRUE;
    t[i].row[i].v := FALSE;
    WITH
      a = s.edges[i].a DIV 3,
      b = s.edges[i].b DIV 3,
      edge1Fixed = s.fixed[a] AND s.fixed[b]
    DO
      FOR j := i+1 TO LAST(s.edges^) DO
        WITH
          c = s.edges[j].a DIV 3,
          d = s.edges[j].b DIV 3,
          edge2Fixed = s.fixed[c] AND s.fixed[d],
          pair = t[i].row[j]
        DO
          pair.v := NOT (a = c OR a = d OR b = c OR b = d OR
                         s.sameTetrahedron(a, b, c, d) OR
                         edge1Fixed AND edge2Fixed);
          pair.flag[Axis.X] := FALSE;
          pair.flag[Axis.Y] := FALSE;
          pair.flag[Axis.Z] := FALSE;
          IF pair.v THEN INC(n) END
        END
      END
    END
  END;
  RETURN n
END BuildEdgeEdgePairs;

BEGIN END PairTable.
