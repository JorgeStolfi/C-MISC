
#include <Simulator.h>



#include <Integrator AS INTG.h>
#include <RKF4Integrator.h>
#include <LRN.h>

#include <Wr.h>
#include <Thread.h>        /*DEBUG*/
#include <stdio.h>

REVEAL
  T == Public BRANDED OBJECT
      s0, s1: REF INTG.State;
      v0, v1: REF INTG.State;
      e1: REF INTG.Error;
      INTG.T intg;
    OVERRIDES
      integrate = Integrate;
    }

PROCEDURE Integrate(
    T g; 
    double *?t;    /* IN: starting time; OUT: final time */
    Position *?p;    /* IN: initial position; OUT: final position */
    Velocity *?v;    /* IN: initial velocity; OUT: final velocity */
    Time dtMin, dtMax; /* Min and max time steps */
    Coord pTol;        /* Max integration error in "p[i]" allowed per step */
    Speed vTol;        /* Max integration error in "v[i]" allowed per step */
    EvalRHSProc evalRHS;            /* Computes acceleration from time+state */
    CheckStepProc checkStep;        /* Checks and handles discrete events. */
  ) == 
  {
    /* Allocate work areas if necessary: */
    if ((g.intg == NULL )) { g.intg = NEW(RKF4Integrator.T); }
    { /* with*/ n == p.nel ) {

      if ((g.s0 == NULL) || (NUMBER(g.s0^) != 2*n )) {
        g.s0 = NEW(REF INTG.State, 2*n);
        g.s1 = NEW(REF INTG.State, 2*n);
        g.v0 = NEW(REF INTG.State, 2*n);
        g.v1 = NEW(REF INTG.State, 2*n);
      }
      affirm(g.s1 != NULL) && (NUMBER(g.s1^) == 2*n , "??");
      affirm(g.v0 != NULL) && (NUMBER(g.v0^) == 2*n , "??");
      affirm(g.v1 != NULL) && (NUMBER(g.v1^) == 2*n , "??");

      void Diff(Time t, READONLY INTG.State s, VAR INTG.Velocity v)
        {
          { /* with*/ 
            sp == SUBARRAY(s, 0, n),
            sv == SUBARRAY(s, n, n),
            vv == SUBARRAY(v, 0, n),
            va == SUBARRAY(v, n, n)
          ) {
            vv = sv;
            evalRHS(t, sp, vv, va);
          };
        } /* Diff */;

      double *?t0, dt, t1, te, rError;
      {
        { /* with*/ 
          s0 == g.s0^, s0p == SUBARRAY(s0, 0, n), s0v == SUBARRAY(s0, n, n),
          v0 == g.v0^, v0v == SUBARRAY(v0, 0, n), v0a == SUBARRAY(v0, n, n),
          s1 == g.s1^, s1p == SUBARRAY(s1, 0, n), s1v == SUBARRAY(s1, n, n),
          v1 == g.v1^, v1v == SUBARRAY(v1, 0, n), v1a == SUBARRAY(v1, n, n),
          e1 == g.e1^, e1p == SUBARRAY(e1, 0, n), e1v == SUBARRAY(e1, n, n)
        ) {
          t0 = t;
          s0p = p; s0v = v;
          v0v = s0v; evalRHS(t0, s0p, v0v, v0a);
          t1 = t0 + dtMin;
          while (1) {
            dt = t1 - t0;
            affirm(dt > 0.0d0 , "??");
            g.intg.step(Diff, t0,s0,v0, t1,s1, e1);
            /* Check accuracy: */
            { /* with*/ 
              pError == LRN.LInfNorm(e1p),
              vError == LRN.LInfNorm(e1v)
            ) {
              rError = max(pError/pTol, vError/vTol);
            }
            /* Compute next ideal step size, from accuracy standpoint: */
            dt = g.intg.adjustStepSize(t1 - t0, rError, 1.0D0, dtMin, dtMax);
            if ((rError > 1.0D0) && (t1 - t0 > dtMin )) {
              /* Too much error, retry with smaller step: */
              Trace('*'); /*DEBUG*/
              affirm(dt == dtMin) || (t0 + dt < t1 , "??");
              t1 = max(t0 + dtMin, min(t0 + dt, t1 - dtMin))
            } else {
              /* Look for discrete events: */
              evalRHS(t1,  s1p, s1v, /*OUT*/ v1a);
              te = checkStep(t0, s0p, s0v, v0a, t1, s1p, s1v, v1a);
              if ((te <= t0 )) {
                /* Stop at "t0": */
                Trace('!');
                EXIT
              } else if ((te < t1 )) {
                /* Redo step partially, up to "te": */
                Trace('<');
                t1 = te
              } else {
                /* Step accepted, advance: */
                Trace('-');
                t0 = t1;
                s0 = s1;
                v0v = s1v;
                v0a = v1a;
                if ((te == t1 )) {
                  /* Stop at "t1" */
                  Trace('!');
                  EXIT
                } else {
                  t1 = min(t0 + dt, te);
                };
              };
            };
          };
        };
      };
    };
  } /* Integrate */;

void Trace(char c)
  <* FATAL Wr.Failure, Thread.Alerted , "??");
  {
    Wr.PutChar(stderr, c);
  } /* Trace */;

{;
} /* Simulator */.
