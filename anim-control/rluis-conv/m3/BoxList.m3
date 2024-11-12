MODULE BoxList;

IMPORT LR3, SystemTopology, SystemDynamics;

PROCEDURE BuildVertexBoxes(bx: T; s: SystemDynamics.T; dt: LONGREAL) =
VAR vertices := s.getVertices();
    min, max: LONGREAL;
BEGIN
  dt := dt/3.0D0;
  FOR i := 0 TO LAST(bx^) DO
    WITH
      a  = vertices[i].a,
      pa = s.getPos0(a),
      va = s.getVel0(a),
      pb = s.getPos(a),
      vb = s.getVel(a),
      qa = LR3.Add(pa, LR3.Scale(dt, va)),
      qb = LR3.Sub(pb, LR3.Scale(dt, vb))
    DO
      IF qb[0] > pb[0] THEN min := pb[0]; max := qb[0]
      ELSE min := qb[0]; max := pb[0] END;
      IF qa[0] < min THEN min := qa[0]
      ELSIF qa[0] > max THEN max := qa[0] END;
      IF pa[0] < min THEN min := pa[0]
      ELSIF pa[0] > max THEN max := pa[0] END;
      bx[i].min[Axis.X] := min;
      bx[i].max[Axis.X] := max;
      
      IF qb[1] > pb[1] THEN min := pb[1]; max := qb[0]
      ELSE min := qb[1]; max := pb[1] END;
      IF qa[1] < min THEN min := qa[1]
      ELSIF qa[1] > max THEN max := qa[1] END;
      IF pa[1] < min THEN min := pa[1]
      ELSIF pa[1] > max THEN max := pa[1] END;
      bx[i].min[Axis.Y] := min;
      bx[i].max[Axis.Y] := max;
      
      IF qb[2] > pb[2] THEN min := pb[2]; max := qb[2]
      ELSE min := qb[2]; max := pb[2] END;
      IF qa[2] < min THEN min := qa[2]
      ELSIF qa[2] > max THEN max := qa[2] END;
      IF pa[2] < min THEN min := pa[2]
      ELSIF pa[2] > max THEN max := pa[2] END;
      bx[i].min[Axis.Z] := min;
      bx[i].max[Axis.Z] := max;
    END
  END
END BuildVertexBoxes;

PROCEDURE BuildEdgeBoxes(bx, vertexBoxes: T; s: SystemTopology.T) =
VAR edges := s.edgeToVertex;
BEGIN
  FOR i := 0 TO LAST(bx^) DO
    WITH
      a = vertexBoxes[edges[i].a],
      b = vertexBoxes[edges[i].b]
    DO
      bx[i].min[Axis.X] := MIN(a.min[Axis.X], b.min[Axis.X]);
      bx[i].max[Axis.X] := MAX(a.max[Axis.X], b.max[Axis.X]);
      bx[i].min[Axis.Y] := MIN(a.min[Axis.Y], b.min[Axis.Y]);
      bx[i].max[Axis.Y] := MAX(a.max[Axis.Y], b.max[Axis.Y]);
      bx[i].min[Axis.Z] := MIN(a.min[Axis.Z], b.min[Axis.Z]);
      bx[i].max[Axis.Z] := MAX(a.max[Axis.Z], b.max[Axis.Z]);
    END
  END
END BuildEdgeBoxes;

PROCEDURE BuildFaceBoxes(bx, edgeBoxes: T; s: SystemTopology.T) =
VAR faces := s.faceToEdge;
BEGIN
  FOR i := 0 TO LAST(bx^) DO
    WITH
      a = edgeBoxes[faces[i].a],
      b = edgeBoxes[faces[i].b]
    DO
      bx[i].min[Axis.X] := MIN(a.min[Axis.X], b.min[Axis.X]);
      bx[i].max[Axis.X] := MAX(a.max[Axis.X], b.max[Axis.X]);
      bx[i].min[Axis.Y] := MIN(a.min[Axis.Y], b.min[Axis.Y]);
      bx[i].max[Axis.Y] := MAX(a.max[Axis.Y], b.max[Axis.Y]);
      bx[i].min[Axis.Z] := MIN(a.min[Axis.Z], b.min[Axis.Z]);
      bx[i].max[Axis.Z] := MAX(a.max[Axis.Z], b.max[Axis.Z]);
    END
  END
END BuildFaceBoxes;

BEGIN END BoxList.
