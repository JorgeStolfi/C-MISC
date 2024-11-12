
      coordFixed: REF BOOLS; /* "fixed[k] == TRUE" iff system coordinate "k" is fixed */
      s.kc.fixed = NEW(REF BOOLS, nCoords);
      for (k = 0;  k < nCoords; k++) {s.coordFixed[k] = FALSE; }

void ReadDynamics(T s, fxName, char *kcName)
  <* FATAL OSError.E, Rd.Failure, Alerted , "??");
  {
    { /* with*/ f == FileRd.Open(stname) ) { ReadState(s, f);    Rd.Close(f); }
    { /* with*/ f == FileRd.Open(suname) ) { ReadSurface(s, f);  Rd.Close(f); }
    BuildGraph(s);
    BuildVertices(s);
    BuildEdges(s);
    BuildInversions(s);
    s.fixed = bool_vec_new(s.n);
    if ((FileExists(fxname) )) {
      { /* with*/ f == FileRd.Open(fxname) ) { ReadFixedVertices(s, f);   Rd.Close(f); }
    } else {
      s.nfixed = 0;
      for (i = 0;  i < s.n; i++) {s.fixed[i] = FALSE;
    }

    s.cmanager = NEW(ConstraintManager.T);
    if ((FileExists(kcname) )) {
      { /* with*/ f == FileRd.Open(kcname) ) { ReadKineticConstraints(s, f); Rd.Close(f); }
    } else {
      s.cmanager.init(3*s.n, 0);
    }
    s.cmanager.start(s.starttime, s.pos,s.vel);

    s.fmanager = NEW(ForceManager.T);
    s.fmanager.init(10);
  } /* ReadDynamics */;

void ReadFixedVertices(T s, FILE *rd)
  <* FATAL Alerted, Rd.Failure, Rd.EndOfFile, Trap, Lex.Error , "??");
  {
    s.nfixed = ReadHeader(rd, "fixed", "fixed == ");
    for (i = 0;  i < s.n; i++) {s.fixed[i] = FALSE; }
    for (f = 1;  f <= s.nfixed;  f++) {
      { /* with*/ i == Lex.Int(rd) ) { s.fixed[i] = TRUE; }
      Lex.Skip(rd);
    }
    ReadFooter(rd, "fixed");
  } /* ReadFixedVertices */;

void ReadKineticConstraints(T s, FILE *rd)
  <* FATAL Alerted, Rd.Failure, Rd.EndOfFile, Trap, Lex.Error , "??");

  KinematicConstraint.T ReadKinetic()
    double *?ta, tb;
        nat k;
        double p, v;
        double xi;
    {
      SystemIO.ReadKinetic(rd, ta, tb, k, p, v, xi);
      { /* with*/ c == NEW(KinematicConstraint.T) ) {
        c.init(k, ta, tb, p, v, xi);
        return c;
      };
    } /* ReadKinetic */;

  {
    { /* with*/ h == ReadHeader(rd, "kinetic", "kinetic == ") ) {
      s.cmanager.init(3*s.n, h);
      for (f = 1;  f <= h;  f++) {
        { /* with*/ c == ReadKinetic() ) {
          s.cmanager.addConstraint(c);
        };
      }
    ReadFooter(rd, "kinetic");
  } /* ReadKineticConstraints */;

