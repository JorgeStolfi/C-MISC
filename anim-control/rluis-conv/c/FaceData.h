
#ifndef FaceData_H
#define FaceData_H



#include <r3.h>
#include <Integrator.h>

typedef
  Data == struct ?? {
      pa, pb, pc, va, vb, vc,
      r3_t n, n1, n2, n3, dn, dn1, dn2, dn3;
    }
  T == struct ?? {_vec
      Data data;
      bool ready;
    }

void ComputeData(VAR Data data, a, b, nat c, pos, Vector vel);
PROCEDURE GetData(T data, nat f,
                  nat a, b, c; pos, Vector vel): Data;
Data GetReadyData(T data, nat f);
void Discard(T data);
;
} /* FaceData */.

#endif
