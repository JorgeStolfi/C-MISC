/* See tup.h */
/* Last edited on 2023-02-21 11:40:25 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <values.h>

#include <affirm.h>
#include <jsfile.h>
#include <quad.h>

#include <stmap.h>
#include <stimage.h>

#include <tup.h>

dpair_vec_t st_vertex_phone_distances(Map *m, double maxDist, phone_vec_t *ph)
  { /* Work areas for {st_compute_distances}: */
    float d[m->nv];
    quad_arc_t e[m->nv];
    float c[2*m->ne];
    int32_t r[m->nv];
    int32_t nr;
    
    dpair_vec_t H = dpair_vec_new(2*m->nv);
    int32_t nH = 0;
    
    st_map_init_costs(m, d, e, c);
    for (uint32_t pj = 0;  pj < ph->ne; pj++)
      { phone_t *pph = &(ph->e[pj]);
        int32_t pi = pph->vertex;
        st_map_compute_costs(m, pi, (float)maxDist, r, &nr, d, e, c);
        for (uint32_t k = 0;  k < nr; k++)
          { int32_t ui = r[k];
            dpair_vec_expand(&H, nH);
            H.e[nH] = (dpair_t){ ui, pj, d[ui] };
            nH++;
          }
        st_map_reset_costs(m, r, nr, d, e, c);
      }
    dpair_vec_trim(&H, nH);
    return H;
  }

void st_recompute_phone_usage
  ( Map *m, 
    double_vec_t *dem, 
    double K,
    double partDist, 
    phone_vec_t *ph,
    double_vec_t *lost 
  )
  { /* Quotient for Gaussian exponent: */
    double S2 = partDist*partDist/log(2.0);
    /* Path-cost at which the captured fraction is negligible: */
    double tinyFrac = 0.005;
    float stopDist = (float)sqrt(-S2*log(tinyFrac));
    
    /* Compute distances from each vertex to each {stopDist}-reachable phone: */
    dpair_vec_t H = st_vertex_phone_distances(m, stopDist, ph);
    
    /* Pass 1: Compute vertex-phone weights, accum tot vertex weight {totw[ui]}. */
    /* Initialize {totw[ui]} with weight of `no phone'. */
    double totw[m->nv];
    { int32_t ui; for (ui = 0; ui < m->nv; ui++) { totw[ui] = 1.0; } }
    int32_t k;
    for (k = 0; k < H.ne; k++)
      { dpair_t *Hk = &(H.e[k]);
        int32_t ui = Hk->ui;
        /* int32_t pj = Hk->pi; */
        double upd = Hk->dist;
        /* Compute weight (attractiveness) {upw} of phone at {pi} to user at {ui}. */
        double upf = K*exp(-upd*upd/S2); /* Limit fraction of {ui} captured by {pi}. */
        double upw = upf/(1.0 - upf);    /* Attractiveness of {pi} to {ui}. */
        totw[ui] += upw;
        /* Replace distance in {H[k]} by weight: */
        Hk->dist = upw;
      }

    /* (Re) allocate the {lost} array if needed: */
    if ((lost->e != NULL) && (lost->ne != m->nv)) { free(lost->e); lost->e = NULL; }
    if (lost->e == NULL) { (*lost) = double_vec_new(m->nv); }

    /* Pass 2: for each phone, collect served demand of nearby vertices: */
    { int32_t pj; for (pj = 0; pj < ph->ne; pj++) { ph->e[pj].usage = 0.0; } }
    for (k = 0; k < H.ne; k++)
      { dpair_t *Hk = &(H.e[k]);
        int32_t ui = Hk->ui;
        int32_t pj = Hk->pi;
        double upw = Hk->dist;
        /* Assign fraction of {ui}'s demand to {pi}: */
        ph->e[pj].usage += upw*dem->e[ui]/totw[ui];
      }

    /* Compute lost demand {lost[ui]}: */
    { int32_t ui; 
      for (ui = 0; ui < m->nv; ui++) 
        { lost->e[ui] = dem->e[ui]/totw[ui]; } 
    }
    
    /* Reclaim matrix: */
    free(H.e);
  }

void st_write_phones(char *name, Map *m, phone_vec_t *ph)
  { FILE *wr = open_write(name, TRUE);
    fprintf(wr, "phones = %d\n", ph->ne);
    int32_t i;
    for(i = 0; i < ph->ne; i++)
      { phone_t *phi = &(ph->e[i]);
        int32_t vi = phi->vertex;
        VertexData *vdi = m->vd[vi];
        Point *pti = &(vdi->p);
        fprintf
          ( wr, 
            "%6d %6d  %5.0f %5.0f  %8.3f\n", 
            i, vi, pti->c[0], pti->c[1], phi->usage
          );
      }
    fclose(wr);
  }

vec_typeimpl(phone_vec_t,phone_vec,phone_t);
vec_typeimpl(dpair_vec_t,dpair_vec,dpair_t);
