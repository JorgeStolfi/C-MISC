    kcRd: Rd.T;             (* Kinetic constraints, or NIL if none. *)
    efRd: Rd.T;             (* External forces, or NIL if none. *)
    evWr: Wr.T;             (* Event trace. *)
      
IF s.kc # NIL THEN
        s.cmanager.addForces(s.M, t, pos, vel, f);
      END;
      
