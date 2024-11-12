
#include <FaceData.h>



#include <r3.h>
#include <LR3Extras.h>
#include <Integrator.h>

void ComputeData(VAR Data data, a, b, nat c, pos, Vector vel)
{
  data.pa = (r3_t){pos[a], pos[a+1], pos[a+2]}
  data.pb = (r3_t){pos[b], pos[b+1], pos[b+2]}
  data.pc = (r3_t){pos[c], pos[c+1], pos[c+2]}
  data.va = (r3_t){vel[a], vel[a+1], vel[a+2]}
  data.vb = (r3_t){vel[b], vel[b+1], vel[b+2]}
  data.vc = (r3_t){vel[c], vel[c+1], vel[c+2]}
  { /* with*/ 
    u == r3_Sub(data.pb, data.pa),
    v == r3_Sub(data.pc, data.pa),
    w == r3_Sub(data.pc, data.pb),
    du == r3_Sub(data.vb, data.va),
    dv == r3_Sub(data.vc, data.va),
    dw == r3_Sub(data.vc, data.vb),
    n == LR3Extras.Cross(u, v),
    dn == r3_Add(LR3Extras.Cross(du, v), LR3Extras.Cross(u, dv))
  ) {
    data.n = n;
    data.n1 = LR3Extras.Cross(n, w);
    data.n2 = LR3Extras.Cross(v, n);
    data.n3 = LR3Extras.Cross(n, u);
    data.dn = dn;
    data.dn1 = r3_Add(LR3Extras.Cross(dn, w), LR3Extras.Cross(n, dw));
    data.dn2 = r3_Add(LR3Extras.Cross(dv, n), LR3Extras.Cross(v, dn));
    data.dn3 = r3_Add(LR3Extras.Cross(dn, u), LR3Extras.Cross(n, du));
  };
} /* ComputeData */;

PROCEDURE GetData(T data, nat f, a, b, nat c, 
                  pos, Vector vel): Data == 
{
  if ((! data[f].ready )) {
    ComputeData(data[f].data, a, b, c, pos, vel);
    data[f].ready = TRUE;
  }
  return data[f].data;
} /* GetData */;

Data GetReadyData(T data, nat f)
{ return data[f].data; } /* GetReadyData */;

void Discard(T data)
{
  for (i = 0;  i <= ((data.nel - 1)|?|MAX_data);  i++) {
    data[i].ready = FALSE;
  };
} /* Discard */;

{; } /* FaceData */.
