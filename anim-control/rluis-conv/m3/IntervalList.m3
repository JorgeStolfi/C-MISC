MODULE IntervalList;

IMPORT SystemTopology, BoxList, PairTable, CandidateList, Candidate;
FROM PairTable IMPORT Axis;

PROCEDURE New(s: SystemTopology.T): T =
BEGIN
  WITH
    nv = 2*NUMBER(s.vertices^),
    ne = 2*NUMBER(s.edges^),
    nf = 2*NUMBER(s.faces^),
    n  = nv + ne + nf,
    itv = NEW(T, n)
  DO
    FOR i := 0 TO nv-1 BY 2 DO
      itv[i+0].object := IntervalObject.Vertex;
      itv[i+1].object := IntervalObject.Vertex;
      itv[i+0].begin := TRUE;
      itv[i+1].begin := FALSE;
      itv[i+0].i := i DIV 2;
      itv[i+1].i := itv[i].i;
    END;
    
    FOR i := nv TO nv+ne-1 BY 2 DO
      itv[i+0].object := IntervalObject.Edge;
      itv[i+1].object := IntervalObject.Edge;
      itv[i+0].begin := TRUE;
      itv[i+1].begin := FALSE;
      itv[i+0].i := (i-nv) DIV 2;
      itv[i+1].i := itv[i].i;
    END;
    
    FOR i := nv+ne TO nv+ne+nf-1 BY 2 DO
      itv[i+0].object := IntervalObject.Face;
      itv[i+1].object := IntervalObject.Face;
      itv[i+0].begin := TRUE;
      itv[i+1].begin := FALSE;
      itv[i+0].i := (i-nv-ne) DIV 2;
      itv[i+1].i := itv[i].i;
    END;

    RETURN itv;
  END
END New;

PROCEDURE Process(itv: T; axis: Axis;
                  vertexBox, edgeBox, faceBox: BoxList.T;
                  facePairs, edgePairs: PairTable.T;
                  VAR clist: CandidateList.T) =

  PROCEDURE GetValue(i: CARDINAL): LONGREAL =
  BEGIN
    WITH
      object = itv[i].object,
      begin  = itv[i].begin,
      index  = itv[i].i
    DO
      IF object = IntervalObject.Vertex THEN
        IF begin THEN RETURN vertexBox[index].min[axis]
        ELSE RETURN vertexBox[index].max[axis] END
      ELSIF object = IntervalObject.Edge THEN
        IF begin THEN RETURN edgeBox[index].min[axis]
        ELSE RETURN edgeBox[index].max[axis] END
      ELSE  (* object = IntervalObject.Face *)
        IF begin THEN RETURN faceBox[index].min[axis]
        ELSE RETURN faceBox[index].max[axis] END
      END
    END
  END GetValue;

  PROCEDURE SetStatus(READONLY itvi, itvj: Interval) =
  VAR pairs := facePairs;
      vertexFace := TRUE;
      lmin, lmax, cmin, cmax: LONGREAL;
      obji := itvi.object;
      objj := itvj.object;
      l := itvi.i;
      c := itvj.i;
  BEGIN
    IF obji = IntervalObject.Edge AND objj = IntervalObject.Edge THEN
      pairs := edgePairs;
      vertexFace := FALSE;
      IF l > c THEN
        l := c;
        c := itvi.i;
      END;
      lmin := edgeBox[l].min[axis];
      lmax := edgeBox[l].max[axis];
      cmin := edgeBox[c].min[axis];
      cmax := edgeBox[c].max[axis];
    ELSIF obji = IntervalObject.Vertex AND objj = IntervalObject.Face THEN
      lmin := vertexBox[l].min[axis];
      lmax := vertexBox[l].max[axis];
      cmin := faceBox[c].min[axis];
      cmax := faceBox[c].max[axis];
    ELSIF obji = IntervalObject.Face AND objj = IntervalObject.Vertex THEN
      l := c;
      c := itvi.i;
      lmin := vertexBox[l].min[axis];
      lmax := vertexBox[l].max[axis];
      cmin := faceBox[c].min[axis];
      cmax := faceBox[c].max[axis];
    ELSE
      RETURN
    END;
    WITH overlaped = PairTable.Overlap(pairs, l, c) DO
      pairs[l].row[c].flag[axis] := lmax > cmin AND cmax > lmin;
      IF NOT overlaped AND
         PairTable.On(pairs, l) AND
         PairTable.Valid(pairs, l, c) AND
         PairTable.Overlap(pairs, l, c) 
      THEN
        CandidateList.Insert(clist, Candidate.T{vertexFace, l, c});
      END
    END
  END SetStatus;

BEGIN
  FOR i := 1 TO LAST(itv^) DO
    WITH x = GetValue(i) DO
      VAR interval := itv[i];
          j := i-1;
      BEGIN
        WHILE j >= 0 AND x < GetValue(j) DO
          SetStatus(interval, itv[j]);
          itv[j+1] := itv[j];
          DEC(j);
        END;
        IF j < i-1 THEN itv[j+1] := interval END
      END
    END
  END
END Process;

PROCEDURE ClearCandidates(VAR clist: CandidateList.T;
                          facePairs, edgePairs: PairTable.T) =
VAR ptr := clist;
    prev: CandidateList.T := NIL;
    pairs: PairTable.T;
BEGIN
  WHILE ptr # NIL DO
    WITH cd = CandidateList.GetInfo(ptr) DO
      IF cd.vertexFace THEN pairs := facePairs ELSE pairs := edgePairs END;
      IF
        PairTable.On(pairs, cd.i) AND
        PairTable.Valid(pairs, cd.i, cd.j) AND
        PairTable.Overlap(pairs, cd.i, cd.j)
      THEN
        prev := ptr
      ELSE
        IF prev # NIL THEN
          CandidateList.RemoveNext(prev)
        ELSE
          clist := CandidateList.GetNext(clist)
        END
      END;
      ptr := CandidateList.GetNext(ptr);
    END
  END
END ClearCandidates;
  
BEGIN END IntervalList.
