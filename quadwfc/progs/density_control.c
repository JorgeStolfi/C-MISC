/* See {density_control.h}. */
/* Last edited on 2013-10-02 03:16:10 by stolfilocal */

#include <density_control.h>

#include <basic.h>
#include <wavefront.h>
#include <wave_prop.h>
#include <geomodel.h>
#include <interact.h>

void break_long_edges
  ( wavefront_t *wf, 
    double tol, 
    double time, 
    geomodel_t *geo
  );
  /* Inserts vertices in edges that are longer than {tol},
    until all edge lengths are {tol} or less. */
    
void contract_short_edges(wavefront_t *wf, double tol);    
  /* Contract edges that are shorter than {tol},
    as long as possible. */
    
void interpolate_samples
  ( sample_t *us,
    sample_t *vs,
    double frac,
    geomodel_t *geo,
    sample_t *ps
  );
  /* Computes a sample {ps} that lies {frac} of the way from {us} to {vs}.
    
    The procedure computes {ps.vel} by interpolating the velocities
    {us.vel} and {vs.vel}, preserving their length; and the position
    {ps.pos} from {us.pos} and {vs.pos}, taking into account the
    curvature of the front as implied by the velocities {us.vel} and
    {vs.vel}.
    
    The interpolation is performed as above only if both samples
    {us,vs} are alive, and {us.sgn == vs.sgn}. Otherwise the sample
    {ps} will be dead, with {pos,vel} undefined.
    
    Sometimes the interpolation fails, eg. if {us} and {vs} are so
    close to a reflector that {ps} would lie on the other side of it.
    Since the discrepancy may involve multiple reflectors, in this
    case the procedure just gives up, and leaves {ps} dead. */

segment_t *interpolate_segments
  ( segment_t *u, 
    segment_t *v, 
    double frac, 
    double time, 
    geomodel_t *geo
  );
  /* Creates a new ray segment {p} that ends at about {frac} of the
    way from ray segments {u} to {v}. The procedure first computes
    {p.prev} by interpolating {u.prev} and {v.prev}. Then the ray is
    traced forward from {p.prev} by {time} to obtain {p.curr}.
    
    If the interpolation of {u.prev} and {v.prev} fails for some
    reason, the procedure returns {NULL}. Otherwise, the segment {p}
    may have a live {prev} but dead {curr}. */
  
void split_quadrilateral(qarc_t a);
  /* Assumes that the face {f = LEFT(a)} is a regular (non-hole) face
    that has become a quadrilateral as a consequence of edge
    splitting. Inserts a new diagonal edge {d} between {p = quad_dst(a)}
    and {q = quad_org(quad_lprev(a))}. The old face record {f} remains as
    {LEFT(a) = LEFT(d)}, and a new face record is created for
    {RIGHT(d)}.  */

void density_control
  ( wavefront_t *wf, 
    double tol, 
    double time, 
    geomodel_t *geo
  )
  {
    break_long_edges(wf, tol, time, geo);
    contract_short_edges(wf, tol);
  }
  
void split_quadrilateral(qarc_t a)
  {
    
    qarc_t b = quad_lnext(a);
    qarc_t c = quad_lnext(b);
    qarc_t d = quad_lnext(c);
    
    /* Make sure that the left face is now a quadrilateral: */
    assert(quad_lnext(d) == a);
    
    /* The old left face of {a} will remain there: */
    face_t *oldf = LEFT(a); 
    
    /* Create a new left face for {b}: */
    face_t *newf = (face_t *)notnull(malloc(sizeof(face_t)), "no mem"); 
    newf->omit = FALSE;
    newf->num = 0; /* To be renumbered later. */
    
    /* Endpoints of the new diagonal edge: */
    segment_t *p = quad_org(b);
    segment_t *q = quad_org(d);
    
    /* Create the new diagonal edge: */
    qarc_t ediag = quad_make_edge();
    
    /* Attach it to {p}: */
    SET_quad_org(ediag, p); quad_splice(ediag, b);
    
    /* Attach it to {q}: */
    SET_quad_dst(ediag, q); quad_splice(quad_sym(ediag), d);
    
    /* Attach {oldf} to its new sides: */
    assert(LEFT(a) == oldf);
    SET_LEFT(ediag, oldf);
    assert(LEFT(d) == oldf);
    
    /* Attach {newf} to its sides: */
    SET_LEFT(quad_sym(ediag), newf); 
    SET_LEFT(c, newf); 
    SET_LEFT(b, newf);
  }

void break_long_edges
  ( wavefront_t *wf, 
    double tol, 
    double time, 
    geomodel_t *geo
  )
  {
    /* Queue of edges to check: */
    qarc_vec_t root = (qarc_vec_t){1, &(wf->a)};
    qarc_vec_t Q = renumber_edges(root);
    int nQ = Q.ne; /* Number of edges in queue. */
    int nC = 0; /* Number of edges already checked. */
    
    while (nC < nQ)
      { /* Here, the checked edges are {Q[0..nC-1]}, 
          the unchecked ones are {Q[nC..nQ-1]}. */
        /* Get next unchecked edge: */
        qarc_t e = Q.e[nC];
        /* Get endpoints and edge length: */
        segment_t *u = quad_org(e);
        segment_t *v = quad_dst(e);
        double len = r3_dist(&(u->curr.pos), &(v->curr.pos));
        if (sample_is_dead(&(u->curr)) || sample_is_dead(&(v->curr)))
          { /* One or both endpoints of {e} are dead, consider {e} OK: */ 
            nC++;
          }
        else if (len <= tol)
          { /* Edge {e} is short enough, OK: */ 
            nC++;
          }
        else
          { /* Edge {e} is too long, try to break it. */
            /* Try to interpolate the segment: */
            segment_t *p = interpolate_segments(u, v, 0.5, time, geo);
            
            if (p == NULL)
              { /* Interpolation failed, e.g. due to reflector interference. */
                /* Consider {e} OK: */
                nC++;
              }
            else
              { /* Checking:  */
                double lenu = r3_dist(&(u->curr.pos), &(p->curr.pos));
                double lenv = r3_dist(&(p->curr.pos), &(v->curr.pos));
                fprintf(stderr, "  %.3f -> %.3f + %.3f\n", len, lenu, lenv);

                /* Add {p} to sample set: */
                sref_vec_expand(&wf->st, wf->ns);
                wf->st.e[wf->ns] = p; wf->ns++;

                /* Get the two adjacent faces: */
                face_t *fl = LEFT(e);
                face_t *fr = RIGHT(e);

                /* Get key edges in the star of {e}: */
                qarc_t a = quad_oprev(e);
                qarc_t b = quad_lnext(a);
                qarc_t c = quad_lnext(e);
                qarc_t d = quad_lnext(c);
                /* Disconnect {e} at the destination end: */
                quad_splice(quad_sym(e), c);
                SET_quad_dst(e, p);
                /* Create the new edge and attach it to {p,v}: */
                qarc_t newe = quad_make_edge(); 
                /* Attach it to {p}: */
                SET_quad_org(newe, p); quad_splice(newe, quad_sym(e));
                /* Attach it to {v}: */
                SET_quad_dst(newe, v); quad_splice(quad_sym(newe), c);

                /* Split the face {LEFT(e)}, if appropriate: */
                if (fl->omit)
                  { /* Left face was a hole, do not split it. */
                    SET_LEFT(newe, fl);
                    assert(LEFT(c) == fl);
                    assert(LEFT(d) == fl);
                  }
                else
                  { split_quadrilateral(e); }  

                /* Split the face {RIGHT(e)}, if appropriate: */
                if (fr->omit)
                  { /* Right face was a hole, do not split it. */
                    SET_LEFT(quad_sym(newe), fr);
                    assert(LEFT(a) == fr);
                    assert(LEFT(b) == fr);
                  }
                else
                  { split_quadrilateral(a); }  

                /* Make sure that slots {Q[nQ..nQ+2]} exist: */
                quad_arc_vec_expand(&Q, nQ+2); 

                /* insert the three new edges in queue: */
                Q.e[nQ] = newe; nQ++;
                Q.e[nQ] = quad_onext(newe); nQ++;
                Q.e[nQ] = quad_oprev(newe); nQ++;

                /* Do not increment {nC}, so that {e} is checked again. */
              }
          }
      }
    free(Q.e);
  }  

void contract_short_edges(wavefront_t *wf, double tol)
  { 
    /* To be implemented. */
  }

segment_t *interpolate_segments
  ( segment_t *u, 
    segment_t *v, 
    double frac, 
    double time, 
    geomodel_t *geo
  )
  {
    if (sample_is_dead(&(u->prev)) || sample_is_dead(&(v->prev)))
      { 
        /* Starting points are not both alive: */
        return NULL;
      }
    else if (u->curr.sgn != v->curr.sgn)
      {
        /* Starting points are at different ray stages: */
        return NULL;
      }
    else
      { 
        /* Interpolate the {prev} samples: */
        sample_t pprev;
        interpolate_samples(&(u->prev), &(v->prev), frac, geo, &pprev);
        if (sample_is_dead(&pprev))
          { return NULL; }
        else
          { segment_t *p = (segment_t *)notnull(malloc(sizeof(segment_t)), "no mem");
            p->num = 0; /* Just in case. */

            /* Ray-trace from {prev} by {time} to obtain {curr}: */
            p->prev = pprev;
            ray_prop(&(p->prev), time, geo, &(p->curr));
            return p;
          }
      }
  }

void interpolate_samples
  ( sample_t *us,
    sample_t *vs,
    double frac,
    geomodel_t *geo,
    sample_t *ps
  )
  {
    if (sample_is_dead(us) || sample_is_dead(vs) || (us->sgn == vs->sgn))
      { /* Unable to interpolate: */ 
        kill_sample(ps);
      }
    else
      { 
        assert(us->imd == vs->imd);

        /* Interpolate the scalar speed, just in case they differ: */
        double carf = 1.0 - frac;
        double speedu2 = r3_norm_sqr(&(us->vel));
        double speedv2 = r3_norm_sqr(&(vs->vel));
        double speed = sqrt(carf*speedu2 + frac*speedv2);
        fprintf(stderr, "    speed = %.3g", speed);
        
        /* Interpolate velocities, taking care to preserve speed: */
        r3_t pdir; /* Direction of {p}'s velocity. */ 
        r3_mix(carf, &us->vel, frac, &vs->vel, &pdir);
        r3_dir(&pdir, &pdir);
        r3_scale(speed, &pdir, &(ps->vel));

        /* Interpolate the sample positions, taking curvature into account: */
        r3_mix(carf, &(us->pos), frac, &(vs->pos), &(ps->pos));
        /* Estimate sine of half the angle between {us->vel} and {vs->vel}: */
        double hpvsin = r3_sin(&(ps->vel), &(vs->vel));
        double hpusin = r3_sin(&(ps->vel), &(us->vel));
        double hvsin = 0.5*(hpvsin + hpusin);
        /* Compute half the distance from {u} to {v}: */
        double hdist = 0.5*r3_dist(&us->pos, &vs->pos);
        /* Adjust position of {p}: */
        double dp = hdist*hvsin;
        fprintf(stderr, " dp = %.3g", dp);
        r3_mix(dp, &pdir, 1.0, &(ps->pos), &(ps->pos));

        /* Find the model layer that contains the interpolated point: */
        int imd = containing_medium(&(ps->pos), geo);
        if ((imd == -1) || (imd != us->imd))
          { 
            /* Interpolated point lies in a different layer. Forget it: */
            fprintf(stderr, "  imd = %d  uimd = %d[tangled]", imd, us->imd);
            kill_sample(ps);
          }
        else
          {
            /* Copy into {ps} the situation of {us} in geophysical model: */
            ps->sgn = us->sgn;
            ps->imd = us->imd;
          }
        fprintf(stderr, "\n");
      }
  }
