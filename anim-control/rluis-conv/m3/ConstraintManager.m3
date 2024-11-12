      coordFixed: REF BOOLS; (* "fixed[k] = TRUE" iff system coordinate "k" is fixed *)
      s.kc.fixed := NEW(REF BOOLS, nCoords);
      FOR k := 0 TO nCoords-1 DO s.coordFixed[k] := FALSE END;

PROCEDURE ReadDynamics(s: T; fxName, kcName: TEXT) =
  <* FATAL OSError.E, Rd.Failure, Alerted *>
  BEGIN
    WITH f = FileRd.Open(stname) DO ReadState(s, f);    Rd.Close(f) END;
    WITH f = FileRd.Open(suname) DO ReadSurface(s, f);  Rd.Close(f) END;
    BuildGraph(s);
    BuildVertices(s);
    BuildEdges(s);
    BuildInversions(s);
    s.fixed := NEW(REF ARRAY OF BOOLEAN, s.n);
    IF FileExists(fxname) THEN
      WITH f = FileRd.Open(fxname) DO ReadFixedVertices(s, f);   Rd.Close(f) END
    ELSE
      s.nfixed := 0;
      FOR i := 0 TO s.n-1 DO s.fixed[i] := FALSE END;
    END;

    s.cmanager := NEW(ConstraintManager.T);
    IF FileExists(kcname) THEN
      WITH f = FileRd.Open(kcname) DO ReadKineticConstraints(s, f); Rd.Close(f) END
    ELSE
      s.cmanager.init(3*s.n, 0);
    END;
    s.cmanager.start(s.starttime, s.pos,s.vel);

    s.fmanager := NEW(ForceManager.T);
    s.fmanager.init(10);
  END ReadDynamics;

PROCEDURE ReadFixedVertices(s: T; rd: Rd.T) =
  <* FATAL Alerted, Rd.Failure, Rd.EndOfFile, Trap, Lex.Error *>
  BEGIN
    s.nfixed := ReadHeader(rd, "fixed", "fixed = ");
    FOR i := 0 TO s.n-1 DO s.fixed[i] := FALSE END;
    FOR f := 1 TO s.nfixed DO
      WITH i = Lex.Int(rd) DO s.fixed[i] := TRUE END;
      Lex.Skip(rd);
    END;
    ReadFooter(rd, "fixed");
  END ReadFixedVertices;

PROCEDURE ReadKineticConstraints(s: T; rd: Rd.T) =
  <* FATAL Alerted, Rd.Failure, Rd.EndOfFile, Trap, Lex.Error *>

  PROCEDURE ReadKinetic(): KinematicConstraint.T =
    VAR ta, tb: LONGREAL;
        k: CARDINAL;
        p, v: LONGREAL;
        xi: LONGREAL;
    BEGIN
      SystemIO.ReadKinetic(rd, ta, tb, k, p, v, xi);
      WITH c = NEW(KinematicConstraint.T) DO
        c.init(k, ta, tb, p, v, xi);
        RETURN c
      END;
    END ReadKinetic;

  BEGIN
    WITH h = ReadHeader(rd, "kinetic", "kinetic = ") DO
      s.cmanager.init(3*s.n, h);
      FOR f := 1 TO h DO
        WITH c = ReadKinetic() DO
          s.cmanager.addConstraint(c)
        END
      END;
    ReadFooter(rd, "kinetic");
  END ReadKineticConstraints;


