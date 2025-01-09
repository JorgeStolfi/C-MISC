/* See {geomodel_pov.h}. */
/* Last edited on 2005-08-26 16:12:06 by stolfi */

#include <geomodel_pov.h>

#include <basic.h>
#include <geomodel.h>
#include <pov_utils.h>

#include <r3.h>
#include <affirm.h>
#include <bool.h>
#include <sign.h>

#include <values.h>
#include <stdio.h>
#include <math.h>

/* INTERNAL PROTOTYPES */

void pov_reflector(FILE *fpov, reflector_t *rf, double radius, char *txv, char *txe, char *txf);
  /* Writes a POV-Ray model of the reflector {rfl} to file {fpov}.
    Edges and vertices are plotted as cylinders and balls with the
    given {radius}. Vertices, edges and cells are painted with
    textures {txv}, {txe}, and {txf}, respectively. */

/* IMPLEMENTATIONS */
    
void pov_geomodel(FILE *fpov, geomodel_t *geo, double radius)  
  {
    int irf;
    for (irf = 0; irf < geo->nrf; irf++)
      { reflector_t *rf = &(geo->rf[irf]);
        char *txv = jsprintf("tx_rf_v[%d]", irf); 
        char *txe = jsprintf("tx_rf_e[%d]", irf); 
        char *txf = jsprintf("tx_rf_f[%d]", irf); 
        pov_reflector(fpov, rf, radius, txv, txe, txf);
        free(txv); free(txe); free(txf);
      }
  }

void pov_reflector(FILE *fpov, reflector_t *rf, double radius, char *txv, char *txe, char *txf)
  {
    int nx = rf->nv[0];
    int xinf = rf->bb[0].end[0], xsup = rf->bb[0].end[1];
    double dx = (xsup - xinf)/(nx-1);
    
    int ny = rf->nv[1];
    int yinf = rf->bb[1].end[0], ysup = rf->bb[1].end[1];
    double dy = (ysup - yinf)/(ny-1);
    
    int ix, iy;
    for (iy = 0; iy < ny; iy++)
      {
        r3_t S; /* Previous point on the previous row. */
        r3_t Q; /* Previous point on the same row. */
        for (ix = 0; ix < nx; ix++)
          { 
            r3_t P; /* Current point on this row. */
            P = (r3_t){{ xinf + ix*dx, yinf + iy*dy, rf->z[ix + nx*iy] }};
            if (iy > 0)
              { r3_t R; /* Current point on previous row. */
                R = (r3_t){{ P.c[0], P.c[1]-dy, rf->z[ix + nx*(iy-1)] }};
                if (ix > 0)
                  { pov_triangle(fpov, &S, &Q, &R, txf);
                    pov_triangle(fpov, &P, &Q, &R, txf);
                  }
                if (radius > 0.0)
                  { if (ix > 0) { pov_cylinder(fpov, &Q, &R, radius, txe); }
                    pov_cylinder(fpov, &R, &P, radius, txe);
                  }
                S = R;
              }
            if (radius > 0.0)
              { if (ix > 0) { pov_cylinder(fpov, &P, &Q, radius, txe); }
                pov_sphere(fpov, &P, radius, txv);
              }
            Q = P;
          }
      }
  }
