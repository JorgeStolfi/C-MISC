
#ifndef EventManager_H
#define EventManager_H



typedef
  T <: Public;
  Public == OBJECT
    METHODS
      detect(
          double t0;
          Position *p0;
          Velocity *v0;
          double t1;
          Position *p1;
          Velocity *v1;
        ): Time;
        /*
          Returns the time of the next event to happen after time "t0", 
          given the states "(p0,v0)" and "(p1,v1)" at times "t0" and "t1>t0".
          
          In general, event detection is only approximate, and its
          precision decreases as the interval "t1-t0" increases.
          Moreover, when the method returns a time "te" greater than
          "t1", that means only that (probably) there are no events in
          the interval "(t0_ t1]", and integration may (hopefully) continue
          from state "(p1,v1)" for some positive time.
          
          In that case, "te" is only a suggested endpoint for the
          next integration step, subject to revision.
        */
        
      advance(
          double t0;
          Position *p0;
          double t1;
          Position *p1;
          Velocity *v;
        ): Event.T;
        /*
          Handles all events that happen as the system moves from
          position "p0" at time "t0" to position "p1" at time "t1".
          Assumes the velocity is constant ("v") in the interval,
          so that the position changes affinely: "p(t) == p0 + (t-t0)*v". 
          
          Handling an event usually means modifying the system 
          of differential equations, and/or some other auxiliary data
          structures.
          
          This method should be called only when "t1-t0" is small
          enough for the final state "(p1,v)" to be considered valid,
          even if it lies a little outside the set of valid states.
       */;
    }
    
T Join(T_vec em);
  /*
    Joins several event managers into a single one.
    The "detect" method will merely call "em[i].detect"
    for all "i" and, return the earliest of the resulting events.
    The "advance" method will call "em[i].advance" for all "i".
    If "em[i]" is NULL. it is assumed to generate no events,
    and require no "advance". */
;
} /* EventManager */.

#endif
