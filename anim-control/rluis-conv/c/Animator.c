
PROCEDURE Simulate(
    T s;
    double t0, t1;       /* Simulation starting and stopping times. */
    EventManager.T em;     /* Detects and handles discrete events. */
    SystemDynamics.T fc;   /* Internal forces and mass matrix */
    FILE *stRd;             /* Initial state. */
    FILE *mmRd;             /* Mass matrix file, or NULL if to be computed. */
    double fps;          /* Frames per second, or 0 if none. */
    double dtMin, dtMax; /* Min and max integration step. */
    double pTol, vTol;   /* Integration error tolerance per coordinate. */ 
    FILE *mmWr;             /* Output factored mass matrix (may be NULL). */
    FILE *stWr;             /* Computed states. */
    bool verbose = FALSE; 
  ) == 
  VAR
    nat frameNumber = 0;   /* Number of next frame */
    double frameTime = t0;    /* Time of next frame. */
  {
    { /* with*/ 
      nCells == s.nCells,
      nNodes == s.nNodes,
      nCoords == 3*s.nNodes,
      state == NEW(State.T).read(stRd)^,
      kc == NEW(KC.T).init(s.nNodes, kcRd, state),
      ef == NEW(EF.T).init(s.nNodes, efRd),
      cf == NEW(CF.T).init(kSpring, mSpring),
      cd == NEW(CD.T).init(s, state.pos),
      cellFixed == ComputeFixedCells(s, kc.fixed^)^,
      f == NEW(REF Force, nCoords)^,
      M == ComputeMassMatrix(s, kc.fixed^)^,
      intgr == NEW(RKF4Integrator.T).init(nCoords)
    ) {    

      PROCEDURE EvalRightHandSide(
          double t; 
          Sim.Position *p;
          Sim.Velocity *v;
          Sim.Acceleration *?a;
        ) == 
        {
          /* Elastic forces: */
          ComputeInternalForces(s, pos, vel, f);
          /* Contact spring forces */
          ct.addSpringForces(pos, vel, f);
          /* External forces */
          ef.addForces(t, pos, vel, f);
          /* Clear out forces on fixed coordinates: */
          for (k = 0;  k <= ((f.nel - 1)|?|MAX_f);  k++) {if ((s.kc.fixed[k] )) { f[k] = 0.0d0; } /*; } */;
          /* Compute reaction forces for kinetic constraints: */
          kc.addForces(M, t, pos, vel, f);
          for (k = 0;  k <= ((f.nel - 1)|?|MAX_f);  k++) {affirm(f[k]*0.0d0 < 0.0d0 , "??"); }
          /* Compute accelerations: */
          SSM.Solve(s.M, f, a);
          for (k = 0;  k <= ((a.nel - 1)|?|MAX_a);  k++) {
            affirm(a[k]*0.0d0 < 0.0d0 , "??");
            if ((s.kc.fixed[k] )) { affirm(a[k] == 0.0d0 , "??");
          };
        } /* EvalRightHandSide */;

      PROCEDURE CheckStep(
          double t0;
          Position *p0;
          Velocity *v0;
          Acceleration *a0;
          double t1;
          Position *?p1;
          Velocity *?v1;
          Acceleration *?a1;
        ): Time == 
        {
          affirm(t0 < t1 , "??");
          affirm(t1 <= frameTime , "??");
          affirm(t1 <= stopTime , "??");
          
          ev = em.next(t0, p0, v0, a0, t1, p1, v1, a1);
          affirm(ev.time > t0 , "??");
          if ((ev.time > t1 )) {
            /* No event foreseen during this step, take it: */
            return LAST(SIM.Time)
          } else if ((t1 - t0 > dtMin )) {
            /*
              A possible event was detected somewhere during this
              step, but since the step is large we don't know for
              sure.  We must try to integrate again until the event.

              To avoid infinite loops here, we stay at least "dtMin" 
              away from the old interval's endpoints.
            */
            { /* with*/ 
              te == max(t0 + dtMin/2.0d0, min(t1 - dtMin/2.0d0, te))
            ) {
              Trace('<');
              return te;
            }
          } else {
            /* 
              One or more events are expected to happen in this step, 
              but the interval it is too small to subdivide. 
              Assume that the events all happen at t1.
            */
            do {
              ev.handle(s);
              ev = em.firstEvent(t0, p0, v0, a0, t1, p1, v1, a1)
            } while (! ((ev.time > t1)))
          
          ;
        } /* CheckStep */;


      Event NextEvent(~)
        {
          { /* with*/ 
            kcEv == kc.firstEvent(t0, t1),
            efEv == ef.firstEvent(t0, t1),
            cdEv == cd.firstEvent(t0, p0, v0, t1, p1, v1),
            cfEv == cf.firstEvent(t0, p0, v0, t1, p1, v1),
            te == min(min(kcEv.time, efEv.time), min(cdEv.time, cfEv.time))
          ) {
            if ((te  )) { ~; }~
          
          ;
        }~

      Run(s, dtMin, dtmax, pTol, vTol, stWr, evWr, verbose);
    };
  } /* Simulate */;

PROCEDURE Init(
    T s; 
    FILE *stRd;   /* Initial state. */
    FILE *mmRd;   /* Mass matrix file, or NULL if to be computed. */
    FILE *kcRd;   /* Kinetic constraints, or NULL if none. */
    FILE *efRd;   /* External forces, or NULL if none. */
    FILE *mmWr;   /* Output factored mass matrix (may be NULL). */
    bool verbose = FALSE; 
  ) == 
  <* FATAL Wr.Failure, Alerted , "??");
  {
    s.verbose = verbose;
    s.nCoords = 3 * s.nNodes;
    { /* with*/ 
      nCoords == s.nCoords
    ) {
      
      s.sim = NEW(Sim.T);
      
      /* Frame times and stopping times:*/
      if ((fps > 0.0d0 )) {
        s.frameInterval = 1.0D0/fps;
        s.frameTime = min(t1, t0 + s.frameInterval)
      } else {
        s.frameTime = t1
        s.frameInterval = t1 - t0;
      }
      s.stopTime = t1;

      /* Initial state: */
      s.state = State.New(s.nNodes);
      if ((verbose )) { Message("Reading initial state..."); }
      State.Read(stRd, s.state);
      affirm(s.state.time == t0 , "??");
      if ((verbose )) { Message("\n"); }
      
      /* Kinematic constraints (including fixed nodes): */
      if ((kcRd != NULL )) {
        if ((verbose )) { Message("Reading kinetic constraints..."); }
        s.kc = KC.Read(kcRd, s.nNodes);
        s.kc.init(s.state.time, s.state.pos, s.state.vel);
        if ((verbose )) { Message("\n"); }
      } else {
        s.kc = NULL;
      }
      
      /* External forces: */
      if ((efRd != NULL )) {
        if ((verbose )) { Message("Reading external forces..."); }
        s.ef = EF.Read(efRd, s.nNodes);
        s.ef.init(s.state.time, ComputeNodeMasses(s)^);
        if ((verbose )) { Message("\n"); }
      } else {
        s.ef = NULL;
      }

      /* Mass matrix: */
      { /* with*/ 
        fp == s.fingerprint()+BoolFingerprint(s.fixed)
      ) {
        if ((mmRd != NULL )) {
          if ((verbose )) { Message("Reading factored mass matrix..."); }
          SSM.Read(mmRd, s.M, fp);
          if ((verbose )) { Message("\n"); }
        } else {
          if ((verbose )) {
            Message("Computing factored mass matrix...");
            ResetCPUTime();
          }
          s.M = ComputeFactoredMassMatrix(s, s.cf.fixed^);
          if ((verbose )) {
            WriteCPUTime(stderr);
            Message("\n");
          };
        }
        if ((mmWr != NULL )) {
          if ((verbose )) { Message("Writing factored mass matrix..."); }
          SSM.Write(mmWr, s.M, fp);
          if ((verbose )) { Message("\n");
        };
      }

      /* Allocate force vector: */
      s.force = NEW(REF Coords, nCoords);
      ;
    };
  } /* Init */;
  
PROCEDURE GetKineticConstraints(
  ): KC.T == 
  {
    { /* with*/ 
      kc == NEW(KC.T)
    ) {
      if ((??));
    };
  } /* GetKineticConstraints */;

PROCEDURE Run(
    T s;
    double dtMin, dtMax;  /* Min and max integration step. */
    double pTol, vTol;    /* Integration error tolerance per coordinate. */ 
    FILE *stWr;              /* Computed states. */
    FILE *evWr;              /* Event trace. */
  ) == 
  <* FATAL Wr.Failure, Alerted , "??");

  PROCEDURE EvalRightHandSide(
      double t; 
      Sim.Position *p;
      Sim.Velocity *v;
      Sim.Acceleration *?a;
    ) == 
    {
      { /* with*/ 
        f == s.force^
      ) {
        ComputeInternalForces(s, pos, vel, f);
        if ((s.ef != NULL )) { s.ef.addForces(t, pos, vel, f); }
        /* Clear out forces on fixed coordinates: */
        for (k = 0;  k <= ((f.nel - 1)|?|MAX_f);  k++) {if ((s.kc.fixed[k] )) { f[k] = 0.0d0; } /*; } */;
        /* Compute reaction forces for other constraints: */
        if ((s.kc != NULL )) {
          s.cmanager.addForces(s.M, t, pos, vel, f);
        }
        for (k = 0;  k <= ((f.nel - 1)|?|MAX_f);  k++) {affirm(f[k]*0.0d0 < 0.0d0 , "??"); }
        /* Compute accelerations: */
        SSM.Solve(s.M, f, a);
        for (k = 0;  k <= ((a.nel - 1)|?|MAX_a);  k++) {
          affirm(a[k]*0.0d0 < 0.0d0 , "??");
          if ((s.kc.fixed[k] )) { affirm(a[k] == 0.0d0 , "??");
        };
      };
    } /* EvalRightHandSide */;
    
  PROCEDURE CheckStep(
      double t0;
      Position *?p0;
      Velocity *?v0;
      Acceleration *?a0;
      double t1;
      Position *?p1;
      Velocity *?v1;
      Acceleration *?a1;
    ): Time == 
    {
      
???
                /*
                  A possible event was detected somewhere during this
                  step, but since the step is large we don't know for
                  sure.  We must try to integrate again until the event.

                  If we integrated blindly up to the predicted time, we
                  might get bogged down in an infinite sequence of steps
                  trying to nail dow some event..

                  To avoid this problem, we set the new integration
                  limit to a trifle before or after the predicted event time
                  "te".  The fudge amount is proportional to a crude
                  estimate of the uncertainty in "te".  For the same
                  reason, we also stay at least "dtMin" away from the
                  old interval's endpoints.
                */
                { /* with*/ 
                  dte == te - t0,
                  dtf == t1 - te,
                  fudge == max(0.5D0 * dtMin, 0.001D0 * min(dte, dtf))
                ) {
                  if ((dte > dtf )) {
                    Trace('<');
                    return max(t0 + dtMin, te - fudge)
                  } else {
                    Trace('>');
                    return min(t1 - dtMin, te + fudge);
                  };
                }
      ;
    } /* CheckStep */;

  {

    if ((s.verbose )) { Message("Running simulation ..."); }



    s.sim.integrate(
      s.state.time,
      s.state.pos,
      s.state.vel,
      dtMin = dtMin, dtMax = dtMax,
      pTol = pTol, vTol = vTol,
      evalRHS = EvalRightHandSide,
      checkState = CheckState,
      checkStep = CheckStep
    );
    if ((s.verbose )) { 
      Util.WriteCPUTime(stderr);
      Message("\n");
    };
  } /* Run */;
