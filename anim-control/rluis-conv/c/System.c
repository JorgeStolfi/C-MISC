
#include <System.h>



#include <
  Wr, FileWr, OSError,
  Util, Integrator, SystemIO, SystemDynamics,
  Collision, BoxList, PairTable, IntervalList, CandidateList,
  Candidate, ContactList, Contact, VertexData, EdgeData, FaceData, 
  VFSpring, EESpring, BezierSearch;
#include <Thread.h>
#include <stdio.h>
#include <Util.h>
#include <PairTable.h>

typedef
  Event == {NewFrame, Constraint, Collision, ContactBreak, Stop}
  
REVEAL
  T == Public BRANDED OBJECT
      VertexData.T dataVa, dataVb;
      EdgeData.T dataEa, dataEb;
      FaceData.T dataFa, dataFb;
      BoxList.T vertexBoxes, edgeBoxes, faceBoxes;
      PairTable.T facePairs, edgePairs;
      IntervalList.T intervalsX, intervalsY, intervalsZ;
      CandidateList.T cdlist;
      ContactList.T contacts;
    OVERRIDES
      init = Init;
      run = Run;
    }

bool *?Collide;


/*--- INITIALIZATION --------------------------------------------------------*/

PROCEDURE Init(T s, 
          nat fps; t1, dtmin, dtmax, tol, g, double spring,
          SystemTopology.FileNames *names; 
          print, verbose, collide: bool
        ) == 
<* FATAL Wr.Failure, Alerted , "??");
{
  NARROW(s, SystemDynamics.T).init(fps, t1, dtmin, dtmax, tol, g, Kspring,
                                 names, print, verbose);
  Collide = collide;

  { /* with*/ 
    np == NUMBER(s.vertices^),
    nl == NUMBER(s.edges^),
    nf == NUMBER(s.faces^)
  ) {
    if ((verbose )) {
      WriteEOL(stderr);
      fprintf(stderr, "External vertices == "); 
      WriteN(stderr, np); WriteEOL(stderr);
      fprintf(stderr, "External edges == "); 
      WriteN(stderr, nl); WriteEOL(stderr);
      fprintf(stderr, "External faces == "); 
      WriteN(stderr, nf); WriteEOL(stderr);
    }
    
    s.dataVa = NEW(VertexData.T, np);
    s.dataVb = NEW(VertexData.T, np);
    s.dataEa = NEW(EdgeData.T, nl);
    s.dataEb = NEW(EdgeData.T, nl);
    s.dataFa = NEW(FaceData.T, nf);
    s.dataFb = NEW(FaceData.T, nf);
    s.vertexBoxes = NEW(BoxList.T, np);
    s.edgeBoxes = NEW(BoxList.T, nl);
    s.faceBoxes = NEW(BoxList.T, nf);
    s.facePairs = PairTable.New(np, nf);
    s.edgePairs = PairTable.New(nl, nl);
    s.intervalsX = IntervalList.New(s);
    s.intervalsY = IntervalList.New(s);
    s.intervalsZ = IntervalList.New(s);

    s.cdlist = NULL;
    s.contacts = NULL;
    
    if ((verbose )) {
      { /* with*/ 
        fp == PairTable.BuildVertexFacePairs(s.facePairs, s),
        lp == PairTable.BuildEdgeEdgePairs(s.edgePairs, s)
      ) {
        WriteEOL(stderr);
        fprintf(stderr, "Vertex/face pairs == ");
        WriteN(stderr, fp); WriteEOL(stderr);
        fprintf(stderr, "Edge/edge pairs == ");
        WriteN(stderr, lp); WriteEOL(stderr);
        fprintf(stderr, "Total pairs == ");
        WriteN(stderr, fp + lp); WriteEOL(stderr);
        WriteEOL(stderr);
        fflush(stderr);
      };
    };
  };
} /* Init */;


/*--- RUN -------------------------------------------------------------------*/

void Run(T s, char *dest, SystemDynamics.UserForce userForce)
<* FATAL OSError.E, Wr.Failure, Alerted , "??");
VAR
  FileWr.T wr;
  double tprev;
  t0prev = -1.0D100;
  Event event;
  g = s.intgr;
  pos0 = g.getPos0();
  vel0 = g.getVel0();
  pos = g.getPos();
  vel = g.getVel();
  Candidate.T candidate;
  bool sign;       /* Spring direction for edge-edge collisions */
  Contact.T contact;
  ContactList.T prevContact;

  Integrator.Event DetectEvent(t0, double t1, VAR double te)
  VAR intgrEvent = Integrator.Event.None;
      t = t1 + 1.0D0;
  {
    if ((t1 >= s.stoptime )) {
      t = s.stoptime;
      intgrEvent = Integrator.Event.Stop;
      event = Event.Stop;
    }
    
    { /* with*/ 
      tframe == tprev + s.frametime
    ) {
      if ((tframe > t0) && (tframe <= t1) && (tframe < t )) {
        t = tframe;
        intgrEvent = Integrator.Event.Simple;
        event = Event.NewFrame;
      };
    }
    
    { /* with*/ 
      tce == s.cmanager.detectEvent(t0, t1)
    ) {
      if ((tce > t0) && (tce < t )) {
        t = tce;
        intgrEvent = Integrator.Event.Simple;
        event = Event.Constraint;
      };
    }
    
    if ((Collide )) {
      { /* with*/ 
        tc == DetectCollision(t0, t1, sign)
      ) {
        if ((tc >= t0) && (tc < t )) {
          t = tc;
          intgrEvent = Integrator.Event.Special;
          event = Event.Collision;
          affirm(t >= t0) && (t <= t1 , "??");
        };
      }
      
      { /* with*/ 
        tb == DetectContactBreak(t0, t1)
      ) {
        if ((tb >= t0) && (tb < t )) {
          t = tb;
          intgrEvent = Integrator.Event.Special;
          event = Event.ContactBreak;
          affirm(t >= t0) && (t <= t1 , "??");
        };
      };
    }
    
    if ((intgrEvent != Integrator.Event.None )) { te = t; }
    return intgrEvent;
  } /* DetectEvent */;
  
  bool TreatEvent(double t)
  {
    CASE event OF
      Event.Stop == > return FALSE;
    | Event.NewFrame == > tprev = t; PrintFrame(t); return FALSE;
    | Event.Constraint == > s.cmanager.treatEvent(pos, vel, t); return TRUE;
    | Event.Collision == > TreatCollision(t); return TRUE;
    | Event.ContactBreak == > TreatContactBreak(t); return TRUE;
    };
  } /* TreatEvent */;
  
  void PrintFrame(double t)
  <* FATAL Wr.Failure, Alerted , "??");
  nat *?ii;
  {
    SystemIO.WriteHeader(wr, "state", "vertices == ", s.n);
    fprintf(wr, "t == "); WriteLR(wr, t);  WriteEOL(wr);
    { /* with*/ 
      d == Util.Digits(((double)s.n))
    ) {
      for (i = 0;  i < s.n; i++) {
        WriteN(wr, i, d); fprintf(wr, ": ");
        ii = 3*i;
        WriteLR(wr, g.getx(ii+0)); WriteSpace(wr);
        WriteLR(wr, g.getx(ii+1)); WriteSpace(wr);
        WriteLR(wr, g.getx(ii+2)); WriteSpace(wr, 2);
        WriteLR(wr, g.getv(ii+0)); WriteSpace(wr);
        WriteLR(wr, g.getv(ii+1)); WriteSpace(wr);
        WriteLR(wr, g.getv(ii+2)); WriteEOL(wr);
      };
    }
    SystemIO.WriteFooter(wr, "state"); WriteEOL(wr);
    fflush(wr);
  } /* PrintFrame */;

  double DetectCollision(t0, double t1, VAR bool sign)
  VAR
    dt = 1.0D0;
    h = t1 - t0;
    h3 = h/3.0D0;
    collisionDetected = FALSE;
    CandidateList.T ptr;
      
  {
    if ((t0 > t0prev )) {
      VAR
        tempV = s.dataVa;
        tempE = s.dataEa;
        tempF = s.dataFa;
      {
        s.dataVa = s.dataVb; s.dataVb = tempV;
        s.dataEa = s.dataEb; s.dataEb = tempE;
        s.dataFa = s.dataFb; s.dataFb = tempF;
      }
      t0prev = t0;
    }

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

    ptr = s.cdlist;
    while (ptr != NULL ) {
      { /* with*/ 
        cd == CandidateList.GetInfo(ptr),
        i == cd.i,
        j == cd.j
      ) {
        if ((cd.vertexFace )) {
          { /* with*/ 
            a == s.faces[j].a,
            b == s.faces[j].b,
            c == s.faces[j].c,
            d == s.vertices[i].a,
            dataVa == VertexData.GetData(s.dataVa, i, d, pos0, vel0),
            dataVb == VertexData.GetData(s.dataVb, i, d, pos,  vel),
            dataFa == FaceData.GetData(s.dataFa, j, a, b, c, pos0, vel0),
            dataFb == FaceData.GetData(s.dataFb, j, a, b, c, pos,  vel),
            t == Collision.DetectVFCollision(dataVa, dataVb, dataFa, dataFb, h3)
          ) {
            if ((t <= dt )) {
              collisionDetected = TRUE;
              candidate = cd;
              dt = t;
            };
          }
        } else {
          { /* with*/ 
            ai == s.edges[i].a,
            bi == s.edges[i].b,
            ci == s.edges[i].c,
            di == s.edges[i].d,
            aj == s.edges[j].a,
            bj == s.edges[j].b,
            cj == s.edges[j].c,
            dj == s.edges[j].d,
            dataEai == EdgeData.GetData(s.dataEa, i, 
                                       ai, bi, ci, di, pos0, vel0),
            dataEbi == EdgeData.GetData(s.dataEb, i, 
                                       ai, bi, ci, di, pos,  vel),
            dataEaj == EdgeData.GetData(s.dataEa, j,
                                       aj, bj, cj, dj, pos0, vel0),
            dataEbj == EdgeData.GetData(s.dataEb, j, 
                                       aj, bj, cj, dj, pos,  vel),
            t == Collision.DetectEECollision(dataEai, dataEbi,
                                            dataEaj, dataEbj, h3,
                                            sign)
          ) {
            if ((t <= dt )) {
              collisionDetected = TRUE;
              candidate = cd;
              dt = t;
            };
          };
        };
      }
      ptr = CandidateList.GetNext(ptr);
    }

    if ((collisionDetected )) {
      return t0 + dt*h
    } else {
      return t0 - 1.0D0;
    };
  } /* DetectCollision */;

  void TreatCollision(<* UNUSED , "??"); double t)
  VAR
    nat a, b, c, d;
    double a1, a2, a3, a4;
  {
    { /* with*/ 
      i == candidate.i,
      j == candidate.j
    ) {
      if ((candidate.vertexFace )) {
        a = s.faces[j].a;
        b = s.faces[j].b;
        c = s.faces[j].c;
        d = s.vertices[i].a;
        { /* with*/ 
          f == NEW(VFSpring.T)
        ) {
          f.init(a, b, c, d, a1, a2, a3, s.Kspring);
          s.fmanager.addForce(f);
          PairTable.SetOn(s.facePairs, i, FALSE);
          PairTable.SetValid(s.facePairs, i, j, FALSE);
          ContactList.Insert(s.contacts, (Contact.T){TRUE, i, j, f, TRUE});
        }
      } else {
        a = s.edges[i].a;
        b = s.edges[i].b;
        c = s.edges[j].a;
        d = s.edges[j].b;
        { /* with*/ 
          f == NEW(EESpring.T)
        ) {
          f.init(a, b, c, d, a1, a2, a3, a4, s.Kspring);
          s.fmanager.addForce(f);
          PairTable.SetValid(s.edgePairs, i, j, FALSE);
          ContactList.Insert(s.contacts, (Contact.T){FALSE, i, j, f, sign});
        };
      };
    };
  } /* TreatCollision */;

  double DetectContactBreak(t0, double t1)
  VAR
    dt = 1.0D0;
    h = t1 - t0;
    h3 = h/3.0D0;
    contactBreak = FALSE;
    ptr = s.contacts;
    ContactList.T prev = NULL;
    
  {
    while (ptr != NULL ) {
      { /* with*/ 
        ct == ContactList.GetInfo(ptr),
        i == ct.i,
        j == ct.j
      ) {
        if ((ct.vertexFace )) {
          { /* with*/ 
            a == s.faces[j].a,
            b == s.faces[j].b,
            c == s.faces[j].c,
            d == s.vertices[i].a,
            dataVa == VertexData.GetData(s.dataVa, i, d, pos0, vel0),
            dataVb == VertexData.GetData(s.dataVb, i, d, pos, vel),
            dataFa == FaceData.GetData(s.dataFa, j, a, b, c, pos0, vel0),
            dataFb == FaceData.GetData(s.dataFb, j, a, b, c, pos, vel),
            t == Collision.DetectVFContactBreak(dataVa, dataVb, dataFa, dataFb,
                                               h3)
          ) {
            if ((t <= dt) && (t != BezierSearch.NoRoot )) {
              affirm(t >= 0.0D0 , "??");
              contactBreak = TRUE;
              contact = ct;
              prevContact = prev;
              dt = t;
            };
          }
        } else {
          { /* with*/ 
            ai == s.edges[i].a,
            bi == s.edges[i].b,
            ci == s.edges[i].c,
            di == s.edges[i].d,
            aj == s.edges[j].a,
            bj == s.edges[j].b,
            cj == s.edges[j].c,
            dj == s.edges[j].d,
            dataEia == EdgeData.GetData(s.dataEa, i, ai, bi, ci, di, 
                                       pos0, vel0),
            dataEib == EdgeData.GetData(s.dataEb, i, ai, bi, ci, di, pos, vel),
            dataEja == EdgeData.GetData(s.dataEa, j, aj, bj, cj, dj, 
                                       pos0, vel0),
            dataEjb == EdgeData.GetData(s.dataEb, j, aj, bj, cj, dj, pos, vel),
            t == Collision.DetectEEContactBreak(dataEia, dataEib, 
                                               dataEja, dataEjb, h3, ct.sign)
          ) {
            if ((t <= dt) && (t != BezierSearch.NoRoot )) {
              affirm(t >= 0.0D0 , "??");
              contactBreak = TRUE;
              contact = ct;
              prevContact = prev;
              dt = t;
            };
          };
        };
      }
      prev = ptr;
      ptr = ContactList.GetNext(ptr);
    }
    
    if ((contactBreak )) {
      return t0 + dt*h
    } else {
      return t0 - 1.0D0;
    };
  } /* DetectContactBreak */;

  void TreatContactBreak(<* UNUSED , "??"); double t)
  {
    s.fmanager.removeForce(contact.force);
    if ((contact.vertexFace )) {
      PairTable.SetOn(s.facePairs, contact.i, TRUE);
      PairTable.SetValid(s.facePairs, contact.i, contact.j, TRUE);
    } else {
      PairTable.SetValid(s.edgePairs, contact.i, contact.j, TRUE);
    }
    if ((prevContact == NULL )) {
      s.contacts = ContactList.GetNext(s.contacts)
    } else {
      ContactList.RemoveNext(prevContact);
    };
  } /* TreatContactBreak */;
  
{
  wr = FileWr.Open(txtcat(dest, ".st"));
  tprev = g.getTime();
  PrintFrame(tprev);

  NARROW(s, SystemDynamics.T).run(userForce, DetectEvent, TreatEvent);
  
  fclose(wr);
} /* Run */;

{; } /* System */.
