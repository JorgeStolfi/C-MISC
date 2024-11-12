MODULE System;

IMPORT
  Wr, FileWr, OSError,
  Util, Integrator, SystemIO, SystemDynamics,
  Collision, BoxList, PairTable, IntervalList, CandidateList,
  Candidate, ContactList, Contact, VertexData, EdgeData, FaceData, 
  VFSpring, EESpring, BezierSearch;
FROM Thread IMPORT Alerted;
FROM Stdio IMPORT stderr;
FROM Util IMPORT WriteN, WriteLR, WriteSpace, WriteEOL;
FROM PairTable IMPORT Axis;

TYPE
  Event = {NewFrame, Constraint, Collision, ContactBreak, Stop};
  
REVEAL
  T = Public BRANDED OBJECT
      dataVa, dataVb: VertexData.T;
      dataEa, dataEb: EdgeData.T;
      dataFa, dataFb: FaceData.T;
      vertexBoxes, edgeBoxes, faceBoxes: BoxList.T;
      facePairs, edgePairs: PairTable.T;
      intervalsX, intervalsY, intervalsZ: IntervalList.T;
      cdlist: CandidateList.T;
      contacts: ContactList.T;
    OVERRIDES
      init := Init;
      run  := Run;
    END;

VAR Collide: BOOLEAN;


(*--- INITIALIZATION --------------------------------------------------------*)

PROCEDURE Init(s: T; 
          fps: CARDINAL; t1, dtmin, dtmax, tol, g, spring: LONGREAL;
          READONLY names: SystemTopology.FileNames; 
          print, verbose, collide: BOOLEAN
        ) =
<* FATAL Wr.Failure, Alerted *>
BEGIN
  NARROW(s, SystemDynamics.T).init(fps, t1, dtmin, dtmax, tol, g, Kspring,
                                 names, print, verbose);
  Collide := collide;

  WITH
    np = NUMBER(s.vertices^),
    nl = NUMBER(s.edges^),
    nf = NUMBER(s.faces^)
  DO
    IF verbose THEN
      WriteEOL(stderr);
      Wr.PutText(stderr, "External vertices = "); 
      WriteN(stderr, np); WriteEOL(stderr);
      Wr.PutText(stderr, "External edges    = "); 
      WriteN(stderr, nl); WriteEOL(stderr);
      Wr.PutText(stderr, "External faces    = "); 
      WriteN(stderr, nf); WriteEOL(stderr);
    END;
    
    s.dataVa      := NEW(VertexData.T, np);
    s.dataVb      := NEW(VertexData.T, np);
    s.dataEa      := NEW(EdgeData.T, nl);
    s.dataEb      := NEW(EdgeData.T, nl);
    s.dataFa      := NEW(FaceData.T, nf);
    s.dataFb      := NEW(FaceData.T, nf);
    s.vertexBoxes := NEW(BoxList.T, np);
    s.edgeBoxes   := NEW(BoxList.T, nl);
    s.faceBoxes   := NEW(BoxList.T, nf);
    s.facePairs   := PairTable.New(np, nf);
    s.edgePairs   := PairTable.New(nl, nl);
    s.intervalsX  := IntervalList.New(s);
    s.intervalsY  := IntervalList.New(s);
    s.intervalsZ  := IntervalList.New(s);

    s.cdlist   := NIL;
    s.contacts := NIL;
    
    IF verbose THEN
      WITH
        fp = PairTable.BuildVertexFacePairs(s.facePairs, s),
        lp = PairTable.BuildEdgeEdgePairs(s.edgePairs, s)
      DO
        WriteEOL(stderr);
        Wr.PutText(stderr, "Vertex/face pairs = ");
        WriteN(stderr, fp); WriteEOL(stderr);
        Wr.PutText(stderr, "Edge/edge pairs   = ");
        WriteN(stderr, lp); WriteEOL(stderr);
        Wr.PutText(stderr, "Total pairs       = ");
        WriteN(stderr, fp + lp); WriteEOL(stderr);
        WriteEOL(stderr);
        Wr.Flush(stderr);
      END
    END
  END
END Init;


(*--- RUN -------------------------------------------------------------------*)

PROCEDURE Run(s: T; dest: TEXT; userForce: SystemDynamics.UserForce) =
<* FATAL OSError.E, Wr.Failure, Alerted *>
VAR
  wr: FileWr.T;
  tprev: LONGREAL;
  t0prev := -1.0D100;
  event: Event;
  g    := s.intgr;
  pos0 := g.getPos0();
  vel0 := g.getVel0();
  pos  := g.getPos();
  vel  := g.getVel();
  candidate: Candidate.T;
  sign: BOOLEAN;       (* Spring direction for edge-edge collisions *)
  contact: Contact.T;
  prevContact: ContactList.T;

  PROCEDURE DetectEvent(t0, t1: LONGREAL; VAR te: LONGREAL): Integrator.Event =
  VAR intgrEvent := Integrator.Event.None;
      t := t1 + 1.0D0;
  BEGIN
    IF t1 >= s.stoptime THEN
      t := s.stoptime;
      intgrEvent := Integrator.Event.Stop;
      event := Event.Stop;
    END;
    
    WITH
      tframe = tprev + s.frametime
    DO
      IF tframe > t0 AND tframe <= t1 AND tframe < t THEN
        t := tframe;
        intgrEvent := Integrator.Event.Simple;
        event := Event.NewFrame;
      END
    END;
    
    WITH
      tce = s.cmanager.detectEvent(t0, t1)
    DO
      IF tce > t0 AND tce < t THEN
        t := tce;
        intgrEvent := Integrator.Event.Simple;
        event := Event.Constraint;
      END
    END;
    
    IF Collide THEN
      WITH
        tc = DetectCollision(t0, t1, sign)
      DO
        IF tc >= t0 AND tc < t THEN
          t := tc;
          intgrEvent := Integrator.Event.Special;
          event := Event.Collision;
          <* ASSERT t >= t0 AND t <= t1 *>
        END
      END;
      
      WITH
        tb = DetectContactBreak(t0, t1)
      DO
        IF tb >= t0 AND tb < t THEN
          t := tb;
          intgrEvent := Integrator.Event.Special;
          event := Event.ContactBreak;
          <* ASSERT t >= t0 AND t <= t1 *>
        END
      END
    END;
    
    IF intgrEvent # Integrator.Event.None THEN te := t END;
    RETURN intgrEvent;
  END DetectEvent;
  
  PROCEDURE TreatEvent(t: LONGREAL): BOOLEAN =
  BEGIN
    CASE event OF
      Event.Stop         => RETURN FALSE;
    | Event.NewFrame     => tprev := t; PrintFrame(t); RETURN FALSE;
    | Event.Constraint   => s.cmanager.treatEvent(pos, vel, t); RETURN TRUE;
    | Event.Collision    => TreatCollision(t); RETURN TRUE;
    | Event.ContactBreak => TreatContactBreak(t); RETURN TRUE;
    END
  END TreatEvent;
  
  PROCEDURE PrintFrame(t: LONGREAL) =
  <* FATAL Wr.Failure, Alerted *>
  VAR ii: CARDINAL;
  BEGIN
    SystemIO.WriteHeader(wr, "state", "vertices = ", s.n);
    Wr.PutText(wr, "t = "); WriteLR(wr, t);  WriteEOL(wr);
    WITH
      d = Util.Digits(FLOAT(s.n, LONGREAL))
    DO
      FOR i := 0 TO s.n-1 DO
        WriteN(wr, i, d); Wr.PutText(wr, ": ");
        ii := 3*i;
        WriteLR(wr, g.getx(ii+0)); WriteSpace(wr);
        WriteLR(wr, g.getx(ii+1)); WriteSpace(wr);
        WriteLR(wr, g.getx(ii+2)); WriteSpace(wr, 2);
        WriteLR(wr, g.getv(ii+0)); WriteSpace(wr);
        WriteLR(wr, g.getv(ii+1)); WriteSpace(wr);
        WriteLR(wr, g.getv(ii+2)); WriteEOL(wr);
      END
    END;
    SystemIO.WriteFooter(wr, "state"); WriteEOL(wr);
    Wr.Flush(wr);
  END PrintFrame;

  PROCEDURE DetectCollision(t0, t1: LONGREAL; VAR sign: BOOLEAN): LONGREAL =
  VAR
    dt := 1.0D0;
    h  := t1 - t0;
    h3 := h/3.0D0;
    collisionDetected := FALSE;
    ptr: CandidateList.T;
      
  BEGIN
    IF t0 > t0prev THEN
      VAR
        tempV := s.dataVa;
        tempE := s.dataEa;
        tempF := s.dataFa;
      BEGIN
        s.dataVa := s.dataVb; s.dataVb := tempV;
        s.dataEa := s.dataEb; s.dataEb := tempE;
        s.dataFa := s.dataFb; s.dataFb := tempF;
      END;
      t0prev := t0;
    END;

    VertexData.Discard(s.dataVb);
    EdgeData.Discard(s.dataEb);
    FaceData.Discard(s.dataFb);

    BoxList.BuildVertexBoxes(s.vertexBoxes, s, h);
    BoxList.BuildEdgeBoxes(s.edgeBoxes, s.vertexBoxes, s);
    BoxList.BuildFaceBoxes(s.faceBoxes, s.edgeBoxes,  s);
    
    IntervalList.Process(s.intervalsX, Axis.X,
                         s.vertexBoxes, s.edgeBoxes, s.faceBoxes,
                         s.facePairs,  s.edgePairs, s.cdlist);
    IntervalList.Process(s.intervalsY, Axis.Y,
                         s.vertexBoxes, s.edgeBoxes, s.faceBoxes,
                         s.facePairs,  s.edgePairs, s.cdlist);
    IntervalList.Process(s.intervalsZ, Axis.Z,
                         s.vertexBoxes, s.edgeBoxes, s.faceBoxes,
                         s.facePairs,  s.edgePairs, s.cdlist);
    IntervalList.ClearCandidates(s.cdlist, s.facePairs, s.edgePairs);

    ptr := s.cdlist;
    WHILE ptr # NIL DO
      WITH
        cd = CandidateList.GetInfo(ptr),
        i = cd.i,
        j = cd.j
      DO
        IF cd.vertexFace THEN
          WITH
            a = s.faces[j].a,
            b = s.faces[j].b,
            c = s.faces[j].c,
            d = s.vertices[i].a,
            dataVa = VertexData.GetData(s.dataVa, i, d, pos0, vel0),
            dataVb = VertexData.GetData(s.dataVb, i, d, pos,  vel),
            dataFa = FaceData.GetData(s.dataFa, j, a, b, c, pos0, vel0),
            dataFb = FaceData.GetData(s.dataFb, j, a, b, c, pos,  vel),
            t = Collision.DetectVFCollision(dataVa, dataVb, dataFa, dataFb, h3)
          DO
            IF t <= dt THEN
              collisionDetected := TRUE;
              candidate := cd;
              dt := t;
            END
          END
        ELSE
          WITH
            ai = s.edges[i].a,
            bi = s.edges[i].b,
            ci = s.edges[i].c,
            di = s.edges[i].d,
            aj = s.edges[j].a,
            bj = s.edges[j].b,
            cj = s.edges[j].c,
            dj = s.edges[j].d,
            dataEai = EdgeData.GetData(s.dataEa, i, 
                                       ai, bi, ci, di, pos0, vel0),
            dataEbi = EdgeData.GetData(s.dataEb, i, 
                                       ai, bi, ci, di, pos,  vel),
            dataEaj = EdgeData.GetData(s.dataEa, j,
                                       aj, bj, cj, dj, pos0, vel0),
            dataEbj = EdgeData.GetData(s.dataEb, j, 
                                       aj, bj, cj, dj, pos,  vel),
            t = Collision.DetectEECollision(dataEai, dataEbi,
                                            dataEaj, dataEbj, h3,
                                            sign)
          DO
            IF t <= dt THEN
              collisionDetected := TRUE;
              candidate := cd;
              dt := t;
            END
          END
        END
      END;
      ptr := CandidateList.GetNext(ptr);
    END;

    IF collisionDetected THEN
      RETURN t0 + dt*h
    ELSE
      RETURN t0 - 1.0D0
    END
  END DetectCollision;

  PROCEDURE TreatCollision(<* UNUSED *> t: LONGREAL) =
  VAR
    a, b, c, d: CARDINAL;
    a1, a2, a3, a4: LONGREAL;
  BEGIN
    WITH
      i = candidate.i,
      j = candidate.j
    DO
      IF candidate.vertexFace THEN
        a := s.faces[j].a;
        b := s.faces[j].b;
        c := s.faces[j].c;
        d := s.vertices[i].a;
        WITH
          f = NEW(VFSpring.T)
        DO
          f.init(a, b, c, d, a1, a2, a3, s.Kspring);
          s.fmanager.addForce(f);
          PairTable.SetOn(s.facePairs, i, FALSE);
          PairTable.SetValid(s.facePairs, i, j, FALSE);
          ContactList.Insert(s.contacts, Contact.T{TRUE, i, j, f, TRUE});
        END
      ELSE
        a := s.edges[i].a;
        b := s.edges[i].b;
        c := s.edges[j].a;
        d := s.edges[j].b;
        WITH
          f = NEW(EESpring.T)
        DO
          f.init(a, b, c, d, a1, a2, a3, a4, s.Kspring);
          s.fmanager.addForce(f);
          PairTable.SetValid(s.edgePairs, i, j, FALSE);
          ContactList.Insert(s.contacts, Contact.T{FALSE, i, j, f, sign});
        END
      END
    END
  END TreatCollision;

  PROCEDURE DetectContactBreak(t0, t1: LONGREAL): LONGREAL =
  VAR
    dt := 1.0D0;
    h  := t1 - t0;
    h3 := h/3.0D0;
    contactBreak := FALSE;
    ptr := s.contacts;
    prev: ContactList.T := NIL;
    
  BEGIN
    WHILE ptr # NIL DO
      WITH
        ct = ContactList.GetInfo(ptr),
        i = ct.i,
        j = ct.j
      DO
        IF ct.vertexFace THEN
          WITH
            a = s.faces[j].a,
            b = s.faces[j].b,
            c = s.faces[j].c,
            d = s.vertices[i].a,
            dataVa = VertexData.GetData(s.dataVa, i, d, pos0, vel0),
            dataVb = VertexData.GetData(s.dataVb, i, d, pos, vel),
            dataFa = FaceData.GetData(s.dataFa, j, a, b, c, pos0, vel0),
            dataFb = FaceData.GetData(s.dataFb, j, a, b, c, pos, vel),
            t = Collision.DetectVFContactBreak(dataVa, dataVb, dataFa, dataFb,
                                               h3)
          DO
            IF t <= dt AND t # BezierSearch.NoRoot THEN
              <* ASSERT t >= 0.0D0 *>
              contactBreak := TRUE;
              contact      := ct;
              prevContact  := prev;
              dt := t;
            END
          END
        ELSE
          WITH
            ai = s.edges[i].a,
            bi = s.edges[i].b,
            ci = s.edges[i].c,
            di = s.edges[i].d,
            aj = s.edges[j].a,
            bj = s.edges[j].b,
            cj = s.edges[j].c,
            dj = s.edges[j].d,
            dataEia = EdgeData.GetData(s.dataEa, i, ai, bi, ci, di, 
                                       pos0, vel0),
            dataEib = EdgeData.GetData(s.dataEb, i, ai, bi, ci, di, pos, vel),
            dataEja = EdgeData.GetData(s.dataEa, j, aj, bj, cj, dj, 
                                       pos0, vel0),
            dataEjb = EdgeData.GetData(s.dataEb, j, aj, bj, cj, dj, pos, vel),
            t = Collision.DetectEEContactBreak(dataEia, dataEib, 
                                               dataEja, dataEjb, h3, ct.sign)
          DO
            IF t <= dt AND t # BezierSearch.NoRoot THEN
              <* ASSERT t >= 0.0D0 *>
              contactBreak := TRUE;
              contact      := ct;
              prevContact  := prev;
              dt := t;
            END
          END
        END
      END;
      prev := ptr;
      ptr := ContactList.GetNext(ptr);
    END;
    
    IF contactBreak THEN
      RETURN t0 + dt*h
    ELSE
      RETURN t0 - 1.0D0
    END
  END DetectContactBreak;

  PROCEDURE TreatContactBreak(<* UNUSED *> t: LONGREAL) =
  BEGIN
    s.fmanager.removeForce(contact.force);
    IF contact.vertexFace THEN
      PairTable.SetOn(s.facePairs, contact.i, TRUE);
      PairTable.SetValid(s.facePairs, contact.i, contact.j, TRUE);
    ELSE
      PairTable.SetValid(s.edgePairs, contact.i, contact.j, TRUE);
    END;
    IF prevContact = NIL THEN
      s.contacts := ContactList.GetNext(s.contacts)
    ELSE
      ContactList.RemoveNext(prevContact)
    END
  END TreatContactBreak;
  
BEGIN
  wr := FileWr.Open(dest & ".st");
  tprev := g.getTime();
  PrintFrame(tprev);

  NARROW(s, SystemDynamics.T).run(userForce, DetectEvent, TreatEvent);
  
  Wr.Close(wr);
END Run;

BEGIN END System.
