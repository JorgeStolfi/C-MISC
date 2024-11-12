
#ifndef VertexData_H
#define VertexData_H



#include <r3.h>
#include <Integrator.h>

typedef
  Data == struct ?? {
      r3_t p, v;
    }
  T == struct ?? {_vec
      Data data;
      bool ready;
    }

void ComputeData(VAR Data data, nat a, pos, Vector vel);
PROCEDURE GetData(T data, nat p, nat a, 
                  pos, Vector vel): Data;
Data GetReadyData(T data, nat p);
void Discard(T data);
;
} /* VertexData */.

#endif
