
void ReadSurface(T s, FILE *rd)
<* FATAL Rd.EndOfFile , "??");
nat *?n;

  void ReadFace(VAR Face f)
  SystemIO.Face *?fIO;
  {
    SystemIO.ReadFace(rd, fIO);
    f.a = 3*fIO.u;
    f.b = 3*fIO.v;
    f.c = 3*fIO.w;
    f.mu1 = fIO.mu1;
    f.mu2 = fIO.mu2;
    f.e = fIO.e;
  } /* ReadFace */;
  
{
  n = ReadHeader(rd, "surface", "faces == ");
  s.faces = Face_vec_new(n);
  for (i = 0;  i < n; i++) {ReadFace(s.faces[i]); }
  ReadFooter(rd, "surface");
} /* ReadSurface */;
