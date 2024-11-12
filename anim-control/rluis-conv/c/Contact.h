
#ifndef Contact_H
#define Contact_H



#include <Rd.h>
#include <Wr.h>
#include <Force.h>

typedef
  T == struct ?? {
      bool vertexFace;
      nat i, j;
      Force.T force;
      bool sign;  /* Spring direction for edge-edge collisions */
        /* 
          sign == TRUE means remove the spring when the edge-edge determinant
          becomes positive again */;
    }

EXCEPTION EndOfInput;

void Read(FILE *rd, VAR T c);
/* Raises EndOfInput */
void Write(FILE *wr, T c);
;
} /* Contact */.

#endif
