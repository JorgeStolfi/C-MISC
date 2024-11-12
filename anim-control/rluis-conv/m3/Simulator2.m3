MODULE Simulator;

IMPORT LR3, Params, Wr, Text, Scan, FloatMode, Lex, Process, 
       Constraint, Force, SystemDynamics, System;
FROM Stdio IMPORT stderr;
FROM Thread IMPORT Alerted;

TYPE
  FileNames = RECORD
      tp: TEXT := NIL;  (* Topology and matrix *)
      su: TEXT := NIL;  (* Surface *)
      st: TEXT := NIL;  (* Initial state *)
      fx: TEXT := NIL;  (* Fixed vertex file *)
      kc: TEXT := NIL;  (* Kinetic constraint file *)
    END;
  
VAR
  s: System.T;
  names: FileNames;
  sfname: TEXT := NIL;  (* Output state file *)
  xi: LONGREAL := 0.0D0;
    
PROCEDURE NullForce(<*UNUSED*> t: LONGREAL): Vectors3D =
BEGIN RETURN force END NullForce;

PROCEDURE GetNumber(): CARDINAL = BEGIN RETURN s.getNumber() END GetNumber;
PROCEDURE GetPos(i: CARDINAL): LR3.T = BEGIN RETURN s.getPos(3*i) END GetPos;
PROCEDURE GetVel(i: CARDINAL): LR3.T = BEGIN RETURN s.getVel(3*i) END GetVel;

PROCEDURE GetVertices(): Vertices =
BEGIN RETURN s.getVertices() END GetVertices;

PROCEDURE GetInitialPos(): Vector =
BEGIN RETURN s.intgr.getPos0() END GetInitialPos;

PROCEDURE GetInitialVel(): Vector =
BEGIN RETURN s.intgr.getVel0() END GetInitialVel;

PROCEDURE GetInitialTime(): LONGREAL =
BEGIN RETURN s.intgr.getTime() END GetInitialTime;

PROCEDURE Getxi(): LONGREAL = BEGIN RETURN xi END Getxi;

PROCEDURE Init() =
VAR t1      := 5.0D0;
    fps     := 30;
    tol     := 0.1D0;
    dt_min  := 0.0001D0;
    dt_max  := 0.01D0;
    g       := 980.0D0;
    spring  := 1000000.0D0;
    print   := FALSE;
    collide := TRUE;
    quiet   := FALSE;
  
  PROCEDURE ReadParams() =
  <* FATAL FloatMode.Trap, Lex.Error, Alerted, Wr.Failure *>
  VAR i := 1;
      op: TEXT;
  BEGIN
    op := Params.Get(i);
    WHILE Text.GetChar(op, 0) = '-' DO
      CASE Text.GetChar(op, 1) OF
        'b' => INC(i); t1     := Scan.LongReal(Params.Get(i));
      | 'c' => collide := FALSE;
      | 'e' => INC(i); tol    := Scan.LongReal(Params.Get(i));
      | 'f' => INC(i); fps    := Scan.Int(Params.Get(i));
      | 'g' => INC(i); g      := Scan.LongReal(Params.Get(i));
      | 'k' => INC(i); spring := Scan.LongReal(Params.Get(i));
      | 'm' => INC(i); dt_min := Scan.LongReal(Params.Get(i));
      | 'n' => INC(i); dt_max := Scan.LongReal(Params.Get(i));
      | 'p' => print := TRUE;
      | 'q' => quiet := TRUE;
      | 'x' => INC(i); xi     := Scan.LongReal(Params.Get(i));
      | 'F' => 
          CASE Text.GetChar(op, 2) OF
           'X' => INC(i); names.fx := Params.Get(i)
          ELSE
            Wr.PutText(stderr, "Invalid option: " & op & "\n");
          END;
      | 'K' => 
          CASE Text.GetChar(op, 2) OF
           'C' => INC(i); names.kc := Params.Get(i)
          ELSE
            Wr.PutText(stderr, "Invalid option: " & op & "\n");
          END;
      | 'S' => 
          CASE Text.GetChar(op, 2) OF
          | 'T' => INC(i); names.st := Params.Get(i);
          | 'F' => INC(i); sfname   := Params.Get(i);
          | 'U' => INC(i); names.su := Params.Get(i);
          ELSE
            Wr.PutText(stderr, "Invalid option: " & op & "\n");
          END;
      ELSE
        Wr.PutText(stderr, "Invalid option: " & op & "\n");
      END;
      INC(i);
      IF i < Params.Count THEN op := Params.Get(i) ELSE RETURN END
    END;
    names.tp := op;
  END ReadParams;
  
  PROCEDURE ShowHelp() =
  <* FATAL Alerted, Wr.Failure *>
  BEGIN
    Wr.PutText(stderr, "Syntax: sim [options] <filename>\n\n" &
                       "Option       Description        Default value\n" &
                       "-b r   final time                   5 s\n" &
                       "-f n   frames per second            30\n" &
                       "-e r   estimated error tolerance    0.1\n" &
                       "-m r   minimum time step            0.0001 s\n" &
                       "-n r   maximum time step            0.01 s\n" &
                       "-g r   gravity aceleration          980 cm/s^2\n" &
                       "-k r   contact spring constant      1000000 g/s^2\n" &
                       "-x r   constraint spring constant   0 s^-1\n" &
                       "-c     ignore collisions\n" &
                       "-p     print mass matrix\n" &
                       "-q     quiet\n\n" &
                       "-FX f  fixed vertex file\n\n" &
                       "-KC f  kinematic constraint file\n\n" &
                       "-ST f  initial state file\n\n" &
                       "-SF f  final state file\n\n" &
                       "-SU f  surface file\n\n" &
                       "Obs.: n is a positive integer\n" &
                       "      r is a positive real\n" &
                       "      f is a file name (without extension)\n");
  END ShowHelp;
  
BEGIN
  IF Params.Count > 1 THEN ReadParams() END;
  IF name = NIL THEN
    ShowHelp();
    Process.Exit(1);    
  ELSE
    s := NEW(System.T);
    s.init(fps, t1, dt_min, dt_max, tol, g, spring,
           names, print, NOT quiet, collide);
    force := NEW(Vectors3D, s.n);
    FOR i := 0 TO LAST(force^) DO force[i] := LR3.Zero() END;
  END
END Init;

PROCEDURE Go(force: SystemDynamics.UserForce) =
BEGIN s.run(name, force) END Go;

PROCEDURE AddForce(f: Force.T) =
BEGIN s.fmanager.addForce(f) END AddForce;

BEGIN END Simulator.
