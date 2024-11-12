# s/BuildMatrixGrad/BuildMatrixSLap/g
# s/GradDot/SLapDot/g
# s/<bool_t.h>/<bool.h>/g
/typedef.*vec_/s/nat/unsigned int/g
/_vec_new(/s/nat/unsigned int/g
/<affirm.h>/a\
#include <nat.h>

