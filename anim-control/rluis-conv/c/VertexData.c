
#include <VertexData.h>



#include <r3.h>
#include <Integrator.h>

void ComputeData(VAR Data data, nat a, pos, Vector vel)
{
  data.p = (r3_t){pos[a], pos[a+1], pos[a+2]}
  data.v = (r3_t){vel[a], vel[a+1], vel[a+2]};
} /* ComputeData */;

PROCEDURE GetData(T data, nat p, nat a, 
                  pos, Vector vel): Data == 
{
  if ((! data[p].ready )) {
    ComputeData(data[p].data, a, pos, vel);
    data[p].ready = TRUE;
  }
  return data[p].data;
} /* GetData */;

Data GetReadyData(T data, nat p)
{ return data[p].data; } /* GetReadyData */;

void Discard(T data)
{
  for (i = 0;  i <= ((data.nel - 1)|?|MAX_data);  i++) {
    data[i].ready = FALSE;
  };
} /* Discard */;

{; } /* VertexData */.
