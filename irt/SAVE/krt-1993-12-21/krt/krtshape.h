#ifndef KRTSHAPE_H
#define KRTSHAPE_H

#include <r3.h>

typedef void t_parms; /* Size and format is shape-dependent */

typedef struct 
  {
    double (*intp) (r3_t *org, r3_t *dir, t_parms *parms);
    void (*nrmp) (r3_t *hit, t_parms *parms, r3_t *nrm);
    void (*prnp) (FILE *f, t_parms *parms);
  }
  t_shape;

typedef struct
  {
    int neg;         /* TRUE if solid is complemented */
    t_shape *shape;  /* Shape of solid object */
    t_parms *parms;  /* Parameters of solid object */
    r3_t mag;        /* Scaling factors. */
    r3x3_t *rot;     /* Rotation matrix. */
    r3_t trn;        /* Translation vector (after scaling and rotation). */
  }
  t_solid;

t_solid krt_translate_solid (t_solid solid, r3_t trn);
  /* Translates a solid by the given vector. */

t_solid krt_rotate_solid (t_solid solid, r3x3_t *rot);
  /* Rotates a solid by the given matrix. */
  /* IMPORTANT: the matrix must be permanently allocated. */
  
t_solid krt_scale_solid (t_solid solid, r3_t mag);
  /* Scales a solid by the given scale factors. */

double krt_int_solid (r3_t *org, r3_t *dir, t_solid solid);

void krt_nrm_solid ( r3_t hit, t_solid solid, r3_t nrm );

void krt_prn_solid ( FILE *f, t_solid solid );

#endif
