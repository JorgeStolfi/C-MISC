#include "krtbox.h"
#include <math.h>
#include <malloc.h>
#include <ioprotos.h>
#include <stdio.h>
#include <krtshape.h>
#include <r3.h>
#include <js.h>

#define FAR_AWAY 99.99E+20

/* internal prototypes */

double krt_int_box( 
    r3_t *pos,       /* origin of ray */
    r3_t *ray,       /* ray vector */
    t_parms *parms   /* box parameters description */
  );

void krt_nrm_box(
    r3_t *pos,       /* point of intersection */
    t_parms *parms,  /* box parameters */
    r3_t *nrm        /* return surface normal */
  );

void krt_prn_box ( FILE *f, t_parms * parms );

/* box class: */

static t_shape boxshape = (t_shape) {krt_int_box, krt_nrm_box, krt_prn_box};

/* box parameters: */

typedef struct
  {
    t_shape *shape;      /* Shape (must be &boxshape) */
    double  sidehit;     /* side intersected     */
    double  xs, ys, zs;  /* size of sides        */
    double  x, y, z;     /* center position      */
  } 
  o_box;

/* create box object */

t_solid krt_make_box( 
  double x, double y, double z, 
  double xs, double ys, double zs
)
{
    o_box *box = (o_box *) malloc(sizeof(o_box));
    box->shape = &boxshape;
    box->x = x;
    box->y = y;
    box->z = z;
    box->xs = xs;
    box->ys = ys;
    box->zs = zs;
    
    return ( (t_solid) { &boxshape, (t_parms *) box} );
}

/* intersection calculation for ray and box */

double krt_int_box( 
    r3_t *org,       /* origin of ray */
    r3_t *dir,       /* ray vector */
    t_parms *parms   /* box parameters description */
)
{
  double s, ss, xhit, yhit, zhit;
  double xadj, yadj, zadj;
  o_box *box = (o_box *) parms;

  assert (box->shape == &boxshape, "krt_int_box: wrong parms!");

  ss = FAR_AWAY;
  /* translate ray origin to object's space */
  xadj = org->c[0] - box->x;
  yadj = org->c[1] - box->y;
  zadj = org->c[2] - box->z;

  /* check x faces */
  if( dir->c[0] != 0 )
    { 
      s = (box->xs - xadj) / dir->c[0];
      if((s > 0) && (s < ss))
        { yhit = fabs(yadj + s * dir->c[1]);
	  zhit = fabs(zadj + s * dir->c[2]);
	  if((yhit < box->ys) && (zhit < box->zs))
            { box->sidehit = 0;
	      ss = s;
	  }
      }
      s = (-box->xs - xadj) / dir->c[0];
      if((s > 0) && (s < ss))
        { yhit = fabs(yadj + s * dir->c[1]);
	  zhit = fabs(zadj + s * dir->c[2]);
	  if((yhit < box->ys) && (zhit < box->zs))
            { box->sidehit = 1;
	      ss = s;
	  }
      }
  }

  /* check y faces */
  if(dir->c[1] != 0)
    { s = (box->ys - yadj) / dir->c[1];
      if((s > 0) && (s < ss))
        { xhit = fabs(xadj + s * dir->c[0]);
	  zhit = fabs(zadj + s * dir->c[2]);
	  if((xhit < box->xs) && (zhit < box->zs))
            { box->sidehit = 2;
	      ss = s;
	    }
        }
      s = (-box->ys - yadj) / dir->c[1];
      if((s > 0) && (s < ss))
        { xhit = fabs(xadj + s * dir->c[0]);
	  zhit = fabs(zadj + s * dir->c[2]);
	  if((xhit < box->xs) && (zhit < box->zs))
            { box->sidehit = 3;
	      ss = s;
	    }
        }
    }

  /* check z faces */
  if(dir->c[2] != 0)
    { s = (box->zs - zadj) / dir->c[2];
      if((s > 0) && (s < ss))
        { xhit = fabs(xadj + s * dir->c[0]);
	  yhit = fabs(yadj + s * dir->c[1]);
	  if((xhit < box->xs) && (yhit < box->ys))
            { box->sidehit = 4;
	      ss = s;
	    }
        }

      s = (-box->zs - zadj) / dir->c[2];
      if((s > 0) && (s < ss))
        { xhit = fabs(xadj + s * dir->c[0]);
	  yhit = fabs(yadj + s * dir->c[1]);
	  if((xhit < box->xs) && (yhit < box->ys))
            { box->sidehit = 5;
	      ss = s;
            }
        }
    }

  if(ss == FAR_AWAY) return(0.0);
  return(ss);
}

/* normal calculation for box */

void krt_nrm_box(
    r3_t *hit,       /* point of intersection */
    t_parms *parms,  /* box parameters */
    r3_t *nrm        /* return surface normal */
)
{
    o_box  *box = (o_box *) parms;

    assert (box->shape == &boxshape, "krt_nrm_box: wrong parms!");

    nrm->c[0] = 0.0;
    nrm->c[1] = 0.0;
    nrm->c[2] = 0.0;
    switch((int)box->sidehit){

	case(0): nrm->c[0] = 1.0;
	         break;
	
	case(1): nrm->c[0] = -1.0;
	         break;
	
	case(2): nrm->c[1] = 1.0;
	         break;
	
	case(3): nrm->c[1] = -1.0;
	         break;
	
	case(4): nrm->c[2] = 1.0;
	         break;
	
	case(5): nrm->c[2] = -1.0;
	         break;
    }
    return;
}

void krt_prn_box (FILE *f, t_parms *parms)
{
    o_box  *box = (o_box *) parms;

    assert (box->shape == &boxshape, "krt_prn_box: wrong parms!");

    fprintf (f, "(box (%f %f %f) (%f %f %f))", 
      box->x, box->y, box->z,
      box->xs, box->ys, box->zs
    );
}
