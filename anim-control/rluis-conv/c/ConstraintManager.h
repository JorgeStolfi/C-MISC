
    FILE *kcRd;             /* Kinetic constraints, or NULL if none. */
    FILE *efRd;             /* External forces, or NULL if none. */
    FILE *evWr;             /* Event trace. */
      
if ((s.kc != NULL )) {
        s.cmanager.addForces(s.M, t, pos, vel, f);
      }
      
