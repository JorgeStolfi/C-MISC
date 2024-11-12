
#include <State.h>



#include <Rd.h>
#include <Wr.h>
#include <Thread.h>
#include <Lex.h>
#include <FloatMode.h>

REVEAL
  T == Public BRANDED OBJECT
    OVERRIDES
      alloc = Alloc;
      read = Read;
      write = Write;
    }

CONST 
  StateFileVersion == "96-12-26";

T Alloc(T s, nat nNodes)
  {
    { /* with*/ NC == 3*nNodes ) {
      if ((s.pos == NULL) || (NUMBER(s.pos^) != NC )) { 
        s.pos = NEW(REF Position, NC);
      }
      if ((s.vel == NULL) || (NUMBER(s.vel^) != NC )) {
        s.vel = NEW(REF Velocity, NC);
      }
      return s;
    };
  } /* Alloc */;

T Read(T s, FILE *rd)
  <* FATAL Rd.Failure, Rd.EndOfFile, Thread.Alerted, FloatMode.Trap, Lex.Error , "??");
  {
    FileFmt.ReadHeader(rd, "state", StateFileVersion);
    s.comment = FileFmt.ReadComment(rd, '|');
    { /* with*/ 
      nNodes == nget_Int(rd, "nodes"),
      NC == 3*nNodes
    ) {
      EVAL s.alloc(nNodes);
      s.time = nget_Double(rd, "time");
      for (i = 0;  i < s.nNodes; i++) {
        { /* with*/ 
          iCheck == Lex.Int(rd),
          kStart = 3*i
        ) {
          Lex.Skip(rd);
          { /* with*/ iCheck == fget_Int(rd) ) { affirm(iCheck == i , "??"); } Lex.Skip(rd);
          fget_Colon(rd);
          for (k = kStart;  k <= kStart + 2;  k++) {s.pos[k] = fget_Double(rd); }
          for (k = kStart;  k <= kStart + 2;  k++) {s.vel[k] = fget_Double(rd);
        };
      };
    }
    FileFmt.ReadFooter(rd, "state");
    return s;
  } /* Read */;
  
void Write(T s, FILE *wr)
  <* FATAL Wr.Failure, Thread.Alerted, FloatMode.Trap, Lex.Error , "??");
  nat *?k;
  {
    FileFmt.WriteHeader(wr, "state", StateFileVersion);
    FileFmt.WriteComment(wr, s.comment, '|');
    { /* with*/ 
      pos == s.pos^,
      vel == s.vel^,
      NC == pos.nel
      nNodes == NC DIV 3
    ) {
      fprintf(wr, "nodes = %d", nNodes);
      NPut.Double(wr, "time", s.time);
      k = 0;
      for (i = 0;  i < s.nNodes; i++) {
        /* Node number: */
        fprintf(wr, "%d"i); FPut.Colon(wr); 
        FPut.Space(wr);
        { /* with*/ kStart == 3*i ) {
          for (k = kStart;  k <= kStart + 2;  k++) {
            FPut.Space(wr);
            FPut.Double(wr, pos[k]);
          }
          FPut.Space(wr);
          for (k = kStart;  k <= kStart + 2;  k++) {
            FPut.Space(wr);
            FPut.Double(wr, vel[k]);
          };
        };
      };
    }
    FileFmt.WriteFooter(wr, "state");
  } /* Write */;

{;
} /* State */.
