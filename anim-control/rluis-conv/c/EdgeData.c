
#include <EdgeData.h>



#include <r3.h>
#include <LR3Extras.h>
#include <Integrator.h>

void ComputeData(VAR Data data, a, b, c, nat d, pos, Vector vel)
{
  data.pa = (r3_t){pos[a], pos[a+1], pos[a+2]}
  data.pb = (r3_t){pos[b], pos[b+1], pos[b+2]}
  data.va = (r3_t){vel[a], vel[a+1], vel[a+2]}
  data.vb = (r3_t){vel[b], vel[b+1], vel[b+2]}
  { /* with*/ 
    pc == (r3_t){pos[c], pos[c+1], pos[c+2]},
    pd == (r3_t){pos[d], pos[d+1], pos[d+2]},
    vc == (r3_t){vel[c], vel[c+1], vel[c+2]},
    vd == (r3_t){vel[d], vel[d+1], vel[d+2]},
    pm == r3_Scale(0.5D0, r3_Add(pc, pd)),
    vm == r3_Scale(0.5D0, r3_Add(vc, vd)),
    ba == r3_Sub(data.pb, data.pa),
    ma == r3_Sub(pm, data.pa),
    dba == r3_Sub(data.vb, data.va),
    dma == r3_Sub(vm, data.va)
  ) {
    data.n = LR3Extras.Cross(ba, ma);
    data.dn = r3_Add(LR3Extras.Cross(dba, ma), LR3Extras.Cross(ba, dma));
  };
} /* ComputeData */;

PROCEDURE GetData(T data, nat l, a, b, c, nat d, 
                  pos, Vector vel): Data == 
{
  if ((TRUE) || (! data[l].ready )) {
    ComputeData(data[l].data, a, b, c, d, pos, vel);
    data[l].ready = TRUE;
  }
  return data[l].data;
} /* GetData */;

Data GetReadyData(T data, nat l)
{ return data[l].data; } /* GetReadyData */;

void Discard(T data)
{
  for (i = 0;  i <= ((data.nel - 1)|?|MAX_data);  i++) {
    data[i].ready = FALSE;
  };
} /* Discard */;

{; } /* EdgeData */.
