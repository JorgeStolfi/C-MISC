
#ifndef Candidate_H
#define Candidate_H



#include <Rd.h>
#include <Wr.h>

typedef
  T == struct ?? {
      bool vertexFace;
      nat i, j;      /* Box indices */;
    }
      
EXCEPTION EndOfInput;

void Read(FILE *rd, VAR T c);
/* Raises EndOfInput */
void Write(FILE *wr, T c);
;
} /* Candidate */.

#endif
