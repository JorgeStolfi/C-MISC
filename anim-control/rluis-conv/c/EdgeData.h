
#ifndef EdgeData_H
#define EdgeData_H



#include <r3.h>
#include <Integrator.h>

typedef
  Data == struct ?? {
      pa, pb, va, vb,
      r3_t n, dn;
    }
  T == struct ?? {_vec
      Data data;
      bool ready;
    }

void ComputeData(VAR Data data, a, b, c, nat d, pos, Vector vel);
PROCEDURE GetData(T data, nat l, a, b, c, nat d,
                  pos, Vector vel): Data;
Data GetReadyData(T data, nat l);
void Discard(T data);
;
} /* EdgeData */.

#endif
