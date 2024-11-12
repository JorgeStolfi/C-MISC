PROCEDURE Simulate(
    s: T;
    t0, t1: LONGREAL;       (* Simulation starting and stopping times. *)
    em: EventManager.T;     (* Detects and handles discrete events. *)
    fc: SystemDynamics.T;   (* Internal forces and mass matrix *)
    stRd: Rd.T;             (* Initial state. *)
    mmRd: Rd.T;             (* Mass matrix file, or NIL if to be computed. *)
    fps: LONGREAL;          (* Frames per second, or 0 if none. *)
    dtMin, dtMax: LONGREAL; (* Min and max integration step. *)
    pTol, vTol: LONGREAL;   (* Integration error tolerance per coordinate. *) 
    mmWr: Wr.T;             (* Output factored mass matrix (may be NIL). *)
    stWr: Wr.T;             (* Computed states. *)
    verbose: BOOL := FALSE; 
  ) =
  VAR
    frameNumber: CARDINAL := 0;   (* Number of next frame *)
    frameTime: LONGREAL := t0;    (* Time of next frame. *)
  BEGIN
    WITH
      nCells = s.nCells,
      nNodes = s.nNodes,
      nCoords = 3*s.nNodes,
      state = NEW(State.T).read(stRd)^,
      kc = NEW(KC.T).init(s.nNodes, kcRd, state),
      ef = NEW(EF.T).init(s.nNodes, efRd),
      cf = NEW(CF.T).init(kSpring, mSpring),
      cd = NEW(CD.T).init(s, state.pos),
      cellFixed = ComputeFixedCells(s, kc.fixed^)^,
      f = NEW(REF Force, nCoords)^,
      M = ComputeMassMatrix(s, kc.fixed^)^,
      intgr = NEW(RKF4Integrator.T).init(nCoords)
    DO    

      PROCEDURE EvalRightHandSide(
          t: LONGREAL; 
          READONLY p: Sim.Position;
          READONLY v: Sim.Velocity;
          VAR a: Sim.Acceleration;
        ) =
        BEGIN
          (* Elastic forces: *)
          ComputeInternalForces(s, pos, vel, f);
          (* Contact spring forces *)
          ct.addSpringForces(pos, vel, f);
          (* External forces *)
          ef.addForces(t, pos, vel, f);
          (* Clear out forces on fixed coordinates: *)
          FOR k := 0 TO LAST(f) DO IF s.kc.fixed[k] THEN f[k] := 0.0d0 END END;
          (* Compute reaction forces for kinetic constraints: *)
          kc.addForces(M, t, pos, vel, f);
          FOR k := 0 TO LAST(f) DO <* ASSERT f[k]*0.0d0 < 0.0d0 *> END;
          (* Compute accelerations: *)
          SSM.Solve(s.M, f, a);
          FOR k := 0 TO LAST(a) DO 
            <* ASSERT a[k]*0.0d0 < 0.0d0 *>
            IF s.kc.fixed[k] THEN <* ASSERT a[k] = 0.0d0 *> END
          END;
        END EvalRightHandSide;

      PROCEDURE CheckStep(
          t0: LONGREAL;
          READONLY p0: Position;
          READONLY v0: Velocity;
          READONLY a0: Acceleration;
          t1: LONGREAL;
          VAR p1: Position;
          VAR v1: Velocity;
          VAR a1: Acceleration;
        ): Time =
        BEGIN
          <* ASSERT t0 < t1 *>
          <* ASSERT t1 <= frameTime *>
          <* ASSERT t1 <= stopTime *>
          
          ev := em.next(t0, p0, v0, a0, t1, p1, v1, a1);
          <* ASSERT ev.time > t0 *>
          IF ev.time > t1 THEN
            (* No event foreseen during this step, take it: *)
            RETURN LAST(SIM.Time)
          ELSIF t1 - t0 > dtMin THEN
            (*
              A possible event was detected somewhere during this
              step, but since the step is large we don't know for
              sure.  We must try to integrate again until the event.

              To avoid infinite loops here, we stay at least "dtMin" 
              away from the old interval's endpoints.
            *)
            WITH
              te = MAX(t0 + dtMin/2.0d0, MIN(t1 - dtMin/2.0d0, te))
            DO
              Trace('<');
              RETURN te
            END
          ELSE
            (* 
              One or more events are expected to happen in this step, 
              but the interval it is too small to subdivide. 
              Assume that the events all happen at t1.
            *)
            REPEAT
              ev.handle(s);
              ev := em.firstEvent(t0, p0, v0, a0, t1, p1, v1, a1)
            UNTIL ev.time > t1
          
          
        END CheckStep;


      PROCEDURE NextEvent(~): Event =
        BEGIN
          WITH
            kcEv = kc.firstEvent(t0, t1),
            efEv = ef.firstEvent(t0, t1),
            cdEv = cd.firstEvent(t0, p0, v0, t1, p1, v1),
            cfEv = cf.firstEvent(t0, p0, v0, t1, p1, v1),
            te = MIN(MIN(kcEv.time, efEv.time), MIN(cdEv.time, cfEv.time))
          DO
            IF te  THEN ~ END~
          
          
        END~

      Run(s, dtMin, dtmax, pTol, vTol, stWr, evWr, verbose)
    END
  END Simulate;

PROCEDURE Init(
    s: T; 
    stRd: Rd.T;   (* Initial state. *)
    mmRd: Rd.T;   (* Mass matrix file, or NIL if to be computed. *)
    kcRd: Rd.T;   (* Kinetic constraints, or NIL if none. *)
    efRd: Rd.T;   (* External forces, or NIL if none. *)
    mmWr: Wr.T;   (* Output factored mass matrix (may be NIL). *)
    verbose: BOOL := FALSE; 
  ) =
  <* FATAL Wr.Failure, Alerted *>
  BEGIN
    s.verbose := verbose;
    s.nCoords := 3 * s.nNodes;
    WITH
      nCoords = s.nCoords
    DO
      
      s.sim := NEW(Sim.T);
      
      (* Frame times and stopping times:*)
      IF fps > 0.0d0 THEN
        s.frameInterval := 1.0D0/fps;
        s.frameTime := MIN(t1, t0 + s.frameInterval)
      ELSE
        s.frameTime := t1
        s.frameInterval := t1 - t0;
      END;
      s.stopTime := t1;

      (* Initial state: *)
      s.state := State.New(s.nNodes);
      IF verbose THEN Message("Reading initial state...") END;
      State.Read(stRd, s.state);
      <* ASSERT s.state.time = t0 *>
      IF verbose THEN Message("\n") END;
      
      (* Kinematic constraints (including fixed nodes): *)
      IF kcRd # NIL THEN
        IF verbose THEN Message("Reading kinetic constraints...") END;
        s.kc := KC.Read(kcRd, s.nNodes);
        s.kc.init(s.state.time, s.state.pos, s.state.vel);
        IF verbose THEN Message("\n") END;
      ELSE
        s.kc := NIL;
      END;
      
      (* External forces: *)
      IF efRd # NIL THEN
        IF verbose THEN Message("Reading external forces...") END;
        s.ef := EF.Read(efRd, s.nNodes);
        s.ef.init(s.state.time, ComputeNodeMasses(s)^);
        IF verbose THEN Message("\n") END;
      ELSE
        s.ef := NIL;
      END;

      (* Mass matrix: *)
      WITH
        fp = s.fingerprint()+BoolFingerprint(s.fixed)
      DO
        IF mmRd # NIL THEN
          IF verbose THEN Message("Reading factored mass matrix...") END;
          SSM.Read(mmRd, s.M, fp);
          IF verbose THEN Message("\n") END;
        ELSE
          IF verbose THEN
            Message("Computing factored mass matrix...");
            ResetCPUTime();
          END;
          s.M := ComputeFactoredMassMatrix(s, s.cf.fixed^);
          IF verbose THEN
            WriteCPUTime(stderr);
            Message("\n")
          END;
        END;
        IF mmWr # NIL THEN
          IF verbose THEN Message("Writing factored mass matrix...") END;
          SSM.Write(mmWr, s.M, fp);
          IF verbose THEN Message("\n") END;
        END
      END;

      (* Allocate force vector: *)
      s.force := NEW(REF Coords, nCoords);
      
    END;
  END Init;
  
PROCEDURE GetKineticConstraints(
  ): KC.T =
  BEGIN
    WITH
      kc = NEW(KC.T)
    DO
      IF 
    END;
  END GetKineticConstraints;

PROCEDURE Run(
    s: T;
    dtMin, dtMax: LONGREAL;  (* Min and max integration step. *)
    pTol, vTol: LONGREAL;    (* Integration error tolerance per coordinate. *) 
    stWr: Wr.T;              (* Computed states. *)
    evWr: Wr.T;              (* Event trace. *)
  ) =
  <* FATAL Wr.Failure, Alerted *>

  PROCEDURE EvalRightHandSide(
      t: LONGREAL; 
      READONLY p: Sim.Position;
      READONLY v: Sim.Velocity;
      VAR a: Sim.Acceleration;
    ) =
    BEGIN
      WITH
        f = s.force^
      DO
        ComputeInternalForces(s, pos, vel, f);
        IF s.ef # NIL THEN s.ef.addForces(t, pos, vel, f) END;
        (* Clear out forces on fixed coordinates: *)
        FOR k := 0 TO LAST(f) DO IF s.kc.fixed[k] THEN f[k] := 0.0d0 END END;
        (* Compute reaction forces for other constraints: *)
        IF s.kc # NIL THEN
          s.cmanager.addForces(s.M, t, pos, vel, f);
        END;
        FOR k := 0 TO LAST(f) DO <* ASSERT f[k]*0.0d0 < 0.0d0 *> END;
        (* Compute accelerations: *)
        SSM.Solve(s.M, f, a);
        FOR k := 0 TO LAST(a) DO 
          <* ASSERT a[k]*0.0d0 < 0.0d0 *>
          IF s.kc.fixed[k] THEN <* ASSERT a[k] = 0.0d0 *> END
        END;
      END;
    END EvalRightHandSide;
    
  PROCEDURE CheckStep(
      t0: LONGREAL;
      VAR p0: Position;
      VAR v0: Velocity;
      VAR a0: Acceleration;
      t1: LONGREAL;
      VAR p1: Position;
      VAR v1: Velocity;
      VAR a1: Acceleration;
    ): Time =
    BEGIN
      
???
                (*
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
                *)
                WITH
                  dte = te - t0,
                  dtf = t1 - te,
                  fudge = MAX(0.5D0 * dtMin, 0.001D0 * MIN(dte, dtf))
                DO
                  IF dte > dtf THEN
                    Trace('<');
                    RETURN MAX(t0 + dtMin, te - fudge)
                  ELSE
                    Trace('>');
                    RETURN MIN(t1 - dtMin, te + fudge)
                  END;
                END
      
    END CheckStep;

  BEGIN

    IF s.verbose THEN Message("Running simulation ...") END;



    s.sim.integrate(
      s.state.time,
      s.state.pos,
      s.state.vel,
      dtMin := dtMin, dtMax := dtMax,
      pTol := pTol, vTol := vTol,
      evalRHS := EvalRightHandSide,
      checkState := CheckState,
      checkStep := CheckStep
    );
    IF s.verbose THEN 
      Util.WriteCPUTime(stderr);
      Message("\n")
    END
  END Run;

