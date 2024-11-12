MODULE Simulator;

IMPORT Integrator AS INTG, RKF4Integrator, LRN;

IMPORT Wr, Thread;        (*DEBUG*)
FROM Stdio IMPORT stderr; (*DEBUG*)

REVEAL
  T = Public BRANDED OBJECT
      s0, s1: REF INTG.State;
      v0, v1: REF INTG.State;
      e1: REF INTG.Error;
      intg: INTG.T;
    OVERRIDES
      integrate := Integrate;
    END;

PROCEDURE Integrate(
    g: T; 
    VAR t: LONGREAL;    (* IN: starting time; OUT: final time *)
    VAR p: Position;    (* IN: initial position; OUT: final position *)
    VAR v: Velocity;    (* IN: initial velocity; OUT: final velocity *)
    dtMin, dtMax: Time; (* Min and max time steps *)
    pTol: Coord;        (* Max integration error in "p[i]" allowed per step *)
    vTol: Speed;        (* Max integration error in "v[i]" allowed per step *)
    evalRHS: EvalRHSProc;            (* Computes acceleration from time+state *)
    checkStep: CheckStepProc;        (* Checks and handles discrete events. *)
  ) =
  BEGIN
    (* Allocate work areas if necessary: *)
    IF g.intg = NIL THEN g.intg := NEW(RKF4Integrator.T) END;
    WITH n = NUMBER(p) DO

      IF g.s0 = NIL OR NUMBER(g.s0^) # 2*n THEN
        g.s0 := NEW(REF INTG.State, 2*n);
        g.s1 := NEW(REF INTG.State, 2*n);
        g.v0 := NEW(REF INTG.State, 2*n);
        g.v1 := NEW(REF INTG.State, 2*n)
      END;
      <* ASSERT g.s1 # NIL AND NUMBER(g.s1^) = 2*n *>
      <* ASSERT g.v0 # NIL AND NUMBER(g.v0^) = 2*n *>
      <* ASSERT g.v1 # NIL AND NUMBER(g.v1^) = 2*n *>

      PROCEDURE Diff(t: Time; READONLY s: INTG.State; VAR v: INTG.Velocity) =
        BEGIN
          WITH
            sp = SUBARRAY(s, 0, n),
            sv = SUBARRAY(s, n, n),
            vv = SUBARRAY(v, 0, n),
            va = SUBARRAY(v, n, n)
          DO
            vv := sv;
            evalRHS(t, sp, vv, va)
          END;
        END Diff;

      VAR t0, dt, t1, te, rError: LONGREAL;
      BEGIN
        WITH
          s0 = g.s0^, s0p = SUBARRAY(s0, 0, n), s0v = SUBARRAY(s0, n, n),
          v0 = g.v0^, v0v = SUBARRAY(v0, 0, n), v0a = SUBARRAY(v0, n, n),
          s1 = g.s1^, s1p = SUBARRAY(s1, 0, n), s1v = SUBARRAY(s1, n, n),
          v1 = g.v1^, v1v = SUBARRAY(v1, 0, n), v1a = SUBARRAY(v1, n, n),
          e1 = g.e1^, e1p = SUBARRAY(e1, 0, n), e1v = SUBARRAY(e1, n, n)
        DO
          t0 := t;
          s0p := p; s0v := v;
          v0v := s0v; evalRHS(t0, s0p, v0v, v0a);
          t1 := t0 + dtMin;
          LOOP
            dt := t1 - t0;
            <* ASSERT dt > 0.0d0 *>
            g.intg.step(Diff, t0,s0,v0, t1,s1, e1);
            (* Check accuracy: *)
            WITH
              pError = LRN.LInfNorm(e1p),
              vError = LRN.LInfNorm(e1v)
            DO
              rError := MAX(pError/pTol, vError/vTol)
            END;
            (* Compute next ideal step size, from accuracy standpoint: *)
            dt := g.intg.adjustStepSize(t1 - t0, rError, 1.0D0, dtMin, dtMax);
            IF rError > 1.0D0 AND t1 - t0 > dtMin THEN
              (* Too much error, retry with smaller step: *)
              Trace('*'); (*DEBUG*)
              <* ASSERT dt = dtMin OR t0 + dt < t1 *>
              t1 := MAX(t0 + dtMin, MIN(t0 + dt, t1 - dtMin))
            ELSE
              (* Look for discrete events: *)
              evalRHS(t1,  s1p, s1v, (*OUT*) v1a);
              te := checkStep(t0, s0p, s0v, v0a, t1, s1p, s1v, v1a);
              IF te <= t0 THEN
                (* Stop at "t0": *)
                Trace('!');
                EXIT
              ELSIF te < t1 THEN
                (* Redo step partially, up to "te": *)
                Trace('<');
                t1 := te
              ELSE
                (* Step accepted, advance: *)
                Trace('-');
                t0 := t1;
                s0 := s1;
                v0v := s1v;
                v0a := v1a;
                IF te = t1 THEN
                  (* Stop at "t1" *)
                  Trace('!');
                  EXIT
                ELSE
                  t1 := MIN(t0 + dt, te)
                END
              END
            END
          END
        END
      END
    END
  END Integrate;

PROCEDURE Trace(c: CHAR) =
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    Wr.PutChar(stderr, c)
  END Trace;

BEGIN
END Simulator.
