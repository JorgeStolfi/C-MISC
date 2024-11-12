PROCEDURE ReadSurface(s: T; rd: Rd.T) =
<* FATAL Rd.EndOfFile *>
VAR n: CARDINAL;

  PROCEDURE ReadFace(VAR f: Face) =
  VAR fIO: SystemIO.Face;
  BEGIN
    SystemIO.ReadFace(rd, fIO);
    f.a   := 3*fIO.u;
    f.b   := 3*fIO.v;
    f.c   := 3*fIO.w;
    f.mu1 := fIO.mu1;
    f.mu2 := fIO.mu2;
    f.e   := fIO.e;
  END ReadFace;
  
BEGIN
  n := ReadHeader(rd, "surface", "faces = ");
  s.faces := NEW(REF ARRAY OF Face, n);
  FOR i := 0 TO n-1 DO ReadFace(s.faces[i]) END;
  ReadFooter(rd, "surface");
END ReadSurface;

