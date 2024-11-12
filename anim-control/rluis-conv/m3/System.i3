INTERFACE System;

IMPORT SystemDynamics;

TYPE
  T <: Public;
  Public = SystemDynamics.T OBJECT METHODS
      init(
          fps: CARDINAL; t1, dtmin, dtmax, tol, g, spring: LONGREAL;
          READONLY names: SystemTopology.FileNames; 
          print, verbose, collide: BOOLEAN
        );
      run(dest: TEXT; userForce: SystemDynamics.UserForce);
    END;

END System.
