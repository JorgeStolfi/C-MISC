/* See tsp_cost.h */
/* Last edited on 2023-03-31 04:20:26 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include <stdlib.h>
#include <values.h>
#include <math.h>

#include <bool.h>
#include <sign.h>
#include <affirm.h>
#include <jsrandom.h>
#include <minu_brent.h>

#include <tsp_cost.h>
#include <tsp_lib.h>

#define SHOW_COSTS (TRUE)
/* Define this as TRUE to print the cost matrix. */

void gen_arc_cost_matrix(int32_t nv, Distr_t distr, int32_t d, double *c)
  {
    int32_t i, j, k;
    switch (distr)
      { case DISTR_EUC_UCUBE:
        case DISTR_EUC_UBALL:
        case DISTR_EUC_USPHE:
        case DISTR_EUC_GAUSS:
          { /* Generate {nv} random points in {R^d}: */
            double x[nv*d];
            for (i = 0; i < nv; i++) { gen_point_coords(distr, d, &(x[i*d])); }
            /* Print them if appropriate: */
            if (SHOW_COSTS) 
              { fprintf(stderr, "vertex positions\n");
                for (i = 0; i < nv; i++) 
                  { double *xi = &(x[i*d]);
                    for (k = 0; k < d; k++)
                      { fprintf(stderr, " %+7.4f", xi[k]); }
                    fprintf(stderr, "\n");
                  }
                fprintf(stderr, "\n");
              }
            /* Compute distance matrix: */
            for (i = 0; i < nv; i++)
              { c[i*nv+i] = 0.0;
                int32_t id = i*d;
                for (j = i+1; j < nv; j++) 
                  { int32_t jd = j*d;
                    double s2 = 0.0;
                    for (k = 0; k < d; k++)
                      { double dxk = x[id+k] - x[jd+k]; s2 += dxk*dxk; }
                    c[i*nv+j] = c[j*nv+i] = sqrt(s2);
                  }
              }
          }
          break;
        case DISTR_GEN_UUNIT:
        case DISTR_GEN_BIMOD:
        case DISTR_GEN_GAUSS:
          { /* Fill cost matrix with random numbers: */
            int32_t i, j;
            for (i = 0; i < nv; i++)
              { for (j = 0; j <= i; j++) 
                  { double v = gen_arc_cost(distr, d);
                    c[i*nv+j] = c[j*nv+i] = v;
                  }
              }
          }
          break;
        case DISTR_GEN_USPHE:
          { /* Fill the cost matrix with Gaussian random numbers: */
            int32_t i, j;
            double S2 = 0;
            for (i = 0; i < nv; i++)
              { for (j = 0; j <= i; j++) 
                  { double v = gen_arc_cost(DISTR_GEN_GAUSS, d);
                    c[i*nv+j] = v;
                    if (i != j) { S2 += v*v; c[j*nv+i] = v; }
                  }
              }
            /* Scale so that the elements off-diagonal have unit norm: */
            double S = sqrt(S2);
            if (S == 0.0) { S = 1.0; }
            for (i = 0; i < nv; i++)
              { for (j = 0; j < i; j++) 
                  { c[i*nv+j] /= S; c[j*nv+i] /= S; }
              }
          }
          break;
        default:
          affirm(FALSE, "bad distr");
      }

    /* Print costs if appropriate: */
    if (SHOW_COSTS)
      { fprintf(stderr, "cost matrix\n");
        for (i = 0; i < nv; i++)
          { for (j = 0; j <= i; j++) 
              { fprintf(stderr, " %+7.4f", c[i*nv+j]); }
            fprintf(stderr, "\n");
          }
        fprintf(stderr, "\n");
      }
  }

void gen_point_coords(Distr_t distr, int32_t d, double *x)
  { switch (distr)
      { case DISTR_EUC_UCUBE: 
          /* Uniform in cube {[-1,+1]^d}: */
          gen_ucube_point_coords(d, x);
          break;
        case DISTR_EUC_UBALL:
          /* Uniform in unit ball: */
          gen_uball_point_coords(d, x);
          break;
        case DISTR_EUC_USPHE:
          /* Uniform in unit ball: */
          { int32_t i;
            double S2;
            do
              { S2 = 0;
                gen_uball_point_coords(d, x);
                for (i = 0; i < d; i++) { S2 += x[i]*x[i]; }
              }
            while (S2 == 0);
            double S = sqrt(S2);
            for (i = 0; i < d; i++) { x[i] /= S; }
          }
          break;
        case DISTR_EUC_GAUSS:
          /* Gaussian with variance 1 in each coord: */
          gen_gauss_point_coords(d, x);
          break;
        default:
          affirm(FALSE, "bad distr");
      }
  }

double gen_arc_cost(Distr_t distr, int32_t d)
  { switch (distr)
      { case DISTR_GEN_UUNIT: 
          /* Uniform in {[0,1]}: */
          return drandom();
        case DISTR_GEN_GAUSS: 
          /* Gaussian: */
          return dgaussrand();
        case DISTR_GEN_BIMOD:
          /* Bimodal in {[0,1]} with {1/d} small values: */
          { double r = 1.999999*drandom()-0.9999995;
            double s = drandom();
            double z = 0.5*acos(r)/M_PI + (s*d < 1.0 ? 0.0 : 0.5);
            /* Make sure it lies in {[0,1]}: */
            while (z < 0.0) { z += 1.0; }
            while (z > 1.0) { z -= 1.0; }
            return z;
          }
        default:
          affirm(FALSE, "bad distr");
          return 0.0;
      }
  }

void gen_ucube_point_coords(int32_t d, double *x)
  {
    int32_t k;
    for (k = 0; k < d; k++)
      { double xk = 2*drandom()-1;
        x[k] = xk;
      }
  }

void gen_gauss_point_coords(int32_t d, double *x)
  {
    int32_t k;
    for (k = 0; k < d; k++) { x[k] = dgaussrand(); }
  }

void gen_uball_point_coords(int32_t d, double *x)
  {
    int32_t k;
    double r2 = 0.0;
    for (k = 0; k < d; k++) 
      { double xk = dgaussrand(); x[k] = xk; r2 += xk*xk; }
    double y = drandom();
    double f = (y <= 0.0 ? y : exp(log(y)/d))/sqrt(r2);
    for (k = 0; k < d; k++) { x[k] *= f; }
  }

void compute_perm_costs(int32_t nv, bool_t tour, double *c, int32_t np, vtx_t *p, double *Z)
  { int32_t k;
    for (k = 0; k < np; k++)
      { vtx_t *v = &(p[k*nv]); 
        /* Compute cost: */
        Z[k] = (tour ? tour_cost(nv, v, c) : perm_cost(nv, v, c));
      }
  }
    
double tour_cost(int32_t nv, vtx_t *v, double *c)
  { double Z = 0;
    int32_t i;
    for (i = 0; i < nv; i++)
      { int32_t x = v[i], y = v[(i+1)%nv];
        Z += c[x*nv+y];
      }
    return Z;
  }

double perm_cost(int32_t nv, vtx_t *v, double *c)
  { double Z = 0;
    int32_t i;
    for (i = 0; i < nv; i++)
      { int32_t x = i, y = v[i];
        Z += c[x*nv+y];
      }
    return Z;
  }

