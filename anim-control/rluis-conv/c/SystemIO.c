
#include <SystemIO.h>



#include <Rd.h>
#include <Wr.h>
#include <Lex.h>
#include <FloatMode.h>
#include <Text.h>
#include <Util.h>
#include <Thread.h>
#include <Util.h>

CONST /* Original: */
  TP_VERSION == "95-06-03";
  ST_VERSION == "95-06-03";
  SU_VERSION == "95-09-14";
  CO_VERSION == "95-07-10";
  TX_VERSION == "95-09-14";
  EN_VERSION == "95-07-06";
  IM_VERSION == "95-07-05";

CONST
  ST_VERSION == "95-06-03";
  SU_VERSION == "95-09-14";
  FX_VERSION == "96-12-04";
  KC_VERSION == "96-12-04";
  TX_VERSION == "95-09-14";
  EN_VERSION == "95-07-06";
  IM_VERSION == "95-07-05";
  MM_VERSION == "96-12-04";

/*--- READING) && (WRITING OF HEADERS) && (FOOTERS ----------------------------*/

char *GetVersion(char *type)
{
  if ((Text.Equal(type, "topology")    )) { return TP_VERSION
  } else if ((Text.Equal(type, "state")       )) { return ST_VERSION
  } else if ((Text.Equal(type, "surface")     )) { return SU_VERSION
  } else if ((Text.Equal(type, "fixed")       )) { return FX_VERSION
  } else if ((Text.Equal(type, "kinetics")    )) { return KC_VERSION
  } else if ((Text.Equal(type, "textures")    )) { return TX_VERSION
  } else if ((Text.Equal(type, "energies")    )) { return EN_VERSION
  } else if ((Text.Equal(type, "impulses")    )) { return IM_VERSION
  } else if ((Text.Equal(type, "matrix")      )) { return MM_VERSION
  } else { affirm(FALSE , "??");
  };
} /* GetVersion */;

/*--- READING OF TYPES ------------------------------------------------------*/

void ReadVectors(FILE *rd, nat n, pos, Vectors3D vel)
<* FATAL Rd.EndOfFile, Rd.Failure, Alerted, FloatMode.Trap, Lex.Error , "??");
{
  for (i = 0;  i < n; i++) {
    EVAL Lex.Int(rd); EVAL Rd.GetChar(rd);
    pos[i][0] = Lex.Double(rd); Lex.Skip(rd);
    pos[i][1] = Lex.Double(rd); Lex.Skip(rd);
    pos[i][2] = Lex.Double(rd); Lex.Skip(rd);
    vel[i][0] = Lex.Double(rd); Lex.Skip(rd);
    vel[i][1] = Lex.Double(rd); Lex.Skip(rd);
    vel[i][2] = Lex.Double(rd); Lex.Skip(rd);
  };
} /* ReadVectors */;

void ReadPlainVectors(FILE *rd, nat n, pos, PlainVectors3D vel)
<* FATAL Rd.EndOfFile, Rd.Failure, Alerted, FloatMode.Trap, Lex.Error , "??");
{
  for (i = 0;  i <= n-1 BY 3;  i++) {
    EVAL Lex.Int(rd); EVAL Rd.GetChar(rd);
    pos[i+0] = Lex.Double(rd); Lex.Skip(rd);
    pos[i+1] = Lex.Double(rd); Lex.Skip(rd);
    pos[i+2] = Lex.Double(rd); Lex.Skip(rd);
    vel[i+0] = Lex.Double(rd); Lex.Skip(rd);
    vel[i+1] = Lex.Double(rd); Lex.Skip(rd);
    vel[i+2] = Lex.Double(rd); Lex.Skip(rd);
  };
} /* ReadPlainVectors */;

void ReadFace(FILE *rd, VAR Face f)
<* FATAL Rd.Failure, Alerted, FloatMode.Trap, Lex.Error , "??");
{
  f.u = Lex.Int(rd);      Lex.Skip(rd);
  f.v = Lex.Int(rd);      Lex.Skip(rd);
  f.w = Lex.Int(rd);      Lex.Skip(rd);
  f.tx = Lex.Int(rd);      Lex.Skip(rd);
  f.mu1 = Lex.Double(rd); Lex.Skip(rd);
  f.mu1 = Lex.Double(rd); Lex.Skip(rd);
  f.e = Lex.Double(rd); Lex.Skip(rd);
  f.e = (1.0D0 + f.e)/2.0D0;
} /* ReadFace */;

void ReadFixed(FILE *rd, VAR nat i)
<* FATAL Rd.Failure, Alerted, FloatMode.Trap, Lex.Error , "??");
{
  i = Lex.Int(rd); Lex.Skip(rd);
} /* ReadFixed */;

PROCEDURE ReadKinetic(FILE *rd, 
  double *?ta, tb; VAR nat k, VAR p, v, double xi,
) == 
nat *?i;
char *?c;
<* FATAL FloatMode.Trap, Rd.Failure, Rd.EndOfFile, Alerted, Lex.Error , "??");
{
  ta = Lex.Double(rd); Lex.Skip(rd);
  tb = Lex.Double(rd); Lex.Skip(rd);
  i = Lex.Int(rd);      Lex.Skip(rd);
  c = Rd.GetChar(rd);   Lex.Skip(rd);
  p = Lex.Double(rd); Lex.Skip(rd);
  v = Lex.Double(rd); Lex.Skip(rd);
  xi = Lex.Double(rd); Lex.Skip(rd);
  if ((c == 'x' )) { k = 3*i+0
  } else if ((c == 'y' )) { k = 3*i+1
  } else if ((c == 'z' )) { k = 3*i+2
  } else { RAISE Lex.Error;
  };
} /* ReadKinetic */;

void ReadTexture(FILE *rd, VAR Texture m)
<* FATAL Rd.Failure, Alerted, FloatMode.Trap, Lex.Error , "??");
{
  m.aR = Lex.Real(rd); Lex.Skip(rd);
  m.aG = Lex.Real(rd); Lex.Skip(rd);
  m.aB = Lex.Real(rd); Lex.Skip(rd);
  m.dR = Lex.Real(rd); Lex.Skip(rd);
  m.dG = Lex.Real(rd); Lex.Skip(rd);
  m.dB = Lex.Real(rd); Lex.Skip(rd);
  m.sR = Lex.Real(rd); Lex.Skip(rd);
  m.sG = Lex.Real(rd); Lex.Skip(rd);
  m.sB = Lex.Real(rd); Lex.Skip(rd);
  m.tR = Lex.Real(rd); Lex.Skip(rd);
  m.tG = Lex.Real(rd); Lex.Skip(rd);
  m.tB = Lex.Real(rd); Lex.Skip(rd);
  m.ir = Lex.Real(rd); Lex.Skip(rd);
  m.n = Lex.Real(rd); Lex.Skip(rd);
} /* ReadTexture */;


/*--- WRITING OF TYPES ------------------------------------------------------*/

PROCEDURE WriteTetrahedron(FILE *wr, READONLY Tetrahedron t, 
                           nat nbase) == 
{;
} /* WriteTetrahedron */;
  
PROCEDURE WriteVectors(FILE *wr, nat n, pos, Vectors3D vel, 
                       int base) == 
<* FATAL Wr.Failure, Alerted , "??");
{
  { /* with*/ d == Util.Digits(((double)n)) ) {
    for (i = 0;  i < n; i++) {
      WriteN(wr, base + i, d); fprintf(wr, ": ");
      WriteLR(wr, pos[i][0]); WriteSpace(wr); 
      WriteLR(wr, pos[i][1]); WriteSpace(wr);
      WriteLR(wr, pos[i][2]); WriteSpace(wr, 2);
      WriteLR(wr, vel[i][0]); WriteSpace(wr);
      WriteLR(wr, vel[i][1]); WriteSpace(wr);
      WriteLR(wr, vel[i][2]); WriteEOL(wr);
    };
  };
} /* WriteVectors */;

PROCEDURE WritePlainVectors(FILE *wr, nat n, pos, PlainVectors3D vel,
                            int base) == 
<* FATAL Wr.Failure, Alerted , "??");
{
  { /* with*/ d == Util.Digits(((double)n DIV 3)) ) {
    for (i = 0;  i <= n-1 BY 3;  i++) {
      WriteN(wr, base + i DIV 3, d); fprintf(wr, ": ");
      WriteLR(wr, pos[i+0]); WriteSpace(wr); 
      WriteLR(wr, pos[i+1]); WriteSpace(wr);
      WriteLR(wr, pos[i+2]); WriteSpace(wr, 2);
      WriteLR(wr, vel[i+0]); WriteSpace(wr);
      WriteLR(wr, vel[i+1]); WriteSpace(wr);
      WriteLR(wr, vel[i+2]); WriteEOL(wr);
    };
  };
} /* WritePlainVectors */;

void WriteFace(FILE *wr, READONLY Face f, nat base)
{
  WriteN(wr, base + f.u); WriteSpace(wr);
  WriteN(wr, base + f.v); WriteSpace(wr);
  WriteN(wr, base + f.w); WriteSpace(wr);
  WriteN(wr, f.tx);       WriteSpace(wr);
  WriteLR(wr, f.mu1);     WriteSpace(wr);
  WriteLR(wr, f.mu2);     WriteSpace(wr);
  WriteLR(wr, f.e);       WriteEOL(wr);
} /* WriteFace */;

void WriteFixed(FILE *wr, nat i, nat base)
{
  WriteN(wr, base + i); WriteEOL(wr);
} /* WriteFixed */;

void WriteKinetic(FILE *wr, ta, double tb, nat k, p, v, double xi)
<* FATAL Wr.Failure, Alerted , "??");
{
  WriteLR(wr, ta); WriteSpace(wr);
  WriteLR(wr, tb); WriteSpace(wr, 2);
  WriteN (wr, k DIV 3); WriteSpace(wr);
  Wr.PutChar(wr, VAL(ORD('x') + (k % 3), char)); WriteSpace(wr, 2);
  WriteLR(wr, p);  WriteSpace(wr);
  WriteLR(wr, v);  WriteSpace(wr);
  WriteLR(wr, xi); WriteSpace(wr);
} /* WriteKinetic */;

void WriteTexture(FILE *wr, READONLY Texture m)
{
  WriteR(wr, m.aR); WriteSpace(wr);
  WriteR(wr, m.aG); WriteSpace(wr);
  WriteR(wr, m.aB); WriteSpace(wr, 2);
  WriteR(wr, m.dR); WriteSpace(wr);
  WriteR(wr, m.dG); WriteSpace(wr);
  WriteR(wr, m.dB); WriteSpace(wr, 2);
  WriteR(wr, m.sR); WriteSpace(wr);
  WriteR(wr, m.sG); WriteSpace(wr);
  WriteR(wr, m.sB); WriteSpace(wr, 2);
  WriteR(wr, m.tR); WriteSpace(wr);
  WriteR(wr, m.tG); WriteSpace(wr);
  WriteR(wr, m.tB); WriteSpace(wr, 2);
  WriteR(wr, m.ir); WriteSpace(wr, 2);
  WriteR(wr, m.n);  WriteEOL(wr);
} /* WriteTexture */;

{; } /* SystemIO */.
