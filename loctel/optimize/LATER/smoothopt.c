/* See smoothopt.h */
/* Last edited on 2007-12-26 13:01:22 by stolfi */

void st_smooth_opt_phones(Map *m, double *dem, double maxDist, Phone **phP, int *nphP)
  { 
    int nph = *nphP, maxph = nph;
    Phone *ph = *phP;
    
    /* Ideal positions and velocities of the phones: */
    Point q[nph];
    Point dq[nph];
    
    /* Time stepping: */
    int nsteps = 1000;     /* Number of time steps. */
    double  dt = 0.1;      /* Time step (sec). */
    
    /* Work vectors: */
    Point f[nph]; /* Force on each phone. */
    
    /* Constants: */
    double B = 1.00;  /* Viscosity coefficient. */
    double C = 1.00;  /* Friction coefficient. */
    
    /* Initialize state with current positions at rest: */
    int i;
    for (pi = 0; pi < nph; pi++)
      { Phone *phi = &(ph[pi]);
        int vi = phi->vertex;
        VertexData *vd = m->vd[vi];
        q[pi] = vd->p;
        dq[pi] = 0.0;
      }
      
    /* Evolve system: */
    double t = 0.0;
    int ti;
    for(ti = 0; ti < nsteps; ti++)
      { /* Compute force on each phone: */
        double tnew = t + dt;
        st_compute_forces(m, dem, q, dq, nph, maxDist, B, C, f);
        /* Update position and velocity, and plot: */
        for(pi = 0; pi < n; pi++) 
          { /* !!! Should use Runge-Kutta! */
            Point *qi = &(q[pi]);
            Point *dqi = &(dq[pi]);
            Point dqnew = r2_mix(*dqi, 1.0, f[pi], dt/M);
            Point dqmid = r2_mix(*dqi, 0.5, dqnew, 0.5);
            Point qnew = r2_mix(*qi, 1.0, dqmid, dt);
            if (r2_norm(f[i]) > fMax) 
              { fprintf(stderr, "phone %d diverged\n", pi);
                it = nsteps+1;
              }
            pswr_segment(psf, qi->c[0], qi->c[1], qnew.c[0], qnew.c[1]);
            pswr_set_fill_color(ps, 0.75, 0.00, 0.00);
            pswr_fill_dot(psf, qi->c[0], qi->c[1],  0.20,  TRUE, FALSE);
            q[pi] = qnew;
            dq[pi] = dqnew;
           }
        t = tnew;
      }

    /* Compute coverage of existing phones: */
    int vcover[m->nv];    /* Vertex coverage by TUP path-cost balls. */
    int ecover[m->ne];    /* Edge coverage by TUP path-cost balls. */
    int uph[nph];
    float dMax[nph];
    int i;
    for (i = 0; i < nph; i++)
      { /* Find vertex nearest to specified center: */
        uph[i] = ph[i].vertex;
        dMax[i] = (float)maxDist; 
      }
    st_compute_coverage(m, uph, dMax, nph, vcover, ecover);
    /* Make list of vertices sorted by increasing X coord: */
    int vseq[m->nv];
    int k;
    for (k = 0; k < m->nv; k++) { vseq[k] = k; }
    shp_heap_sort(vseq, m->nv, xyOrder, +1);
    /* Cover uncovered vertices: */
    for (k = 0; k < m->nv; k++)
      { int vi = vseq[k];
        if (vcover[vi] == 0)
          { int wi = st_cover_vertex(m, vi, maxDist, vcover, ecover); 
            phone_vec_expand(&ph, &maxph, nph);
            ph[nph] = (Phone){ wi, -1.0 };
            nph++;
          }
      }
    phone_vec_trim(&ph, &maxph, nph);
    (*phP) = ph;
    (*nphP) = nph;
  }

void st_compute_phone_usage(Map *m, double *dem, double partDist, Point *q, int nph)
  {
    /* 
      Demand partition model:
      
      Suppose that a phone is installed at vertex {pi}, and let {ui}
      be any other vertex. Let {uh} be the potential demand for calls
      generated in the neighborhood of {ui}.
      
      Even if {pi} was the only phone in town, it would consume only
      some fraction {upf*uh} of the potential demand {uh}. The balance
      represents calls that the user wanted to make, but didn't for
      any reason. Perhaps he didn't know that the phone existed, or
      decided that it was too far away, or went there but found the
      phone busy or broken. 
      
      The parameter {upf}, between 0 and 1, is the /limiting capture
      fraction/ of vertex {pi} for vertex {ui}. Considering all the
      reasons why a potential call may be lost, {upf} is always
      strictly less than 1, even when {ui = pi}. (Note, by the way,
      that an average user nominally located `at' the vertex {pi} is
      usually at some positive distance from the phone at {pi}, and
      may not even know of its existence.)
      
      We can also interpret this situation as a partition of the
      potential demand {uh} between the phone {pi} and a virtual
      `non-phone' that represents lost calls. We can always say that
      the partition is proportional to certain weights {upw} and
      {uNw}, that measure the /attractiveness/ of each choice
      for the user at {ui}.  If we assume that the attractiveness
      of the `no-phone' alternative is fixed, {uNw  = 1}, then 
      the attractiveness of the phone at {pi} is {upw = upf/(1 - upf)}.
      
      Suppose now that a second phone {qi} is installed within useful
      distance {uqd} from vertex {ui}. The the user at {ui} now has
      three available choices: `no-phone' (lost call), {pi}, and {qi}.
      We will assume that the attractiveness of each choice is not
      affected by the availability of the other choices, meaning that
      the potential demand {uh} is partitioned between them
      proportionally to the weights {uNw = 1}, {upw = upf/(1 - upf)},
      and {uqw = uqf/(1 - uqf)}, respectively.
      
      In other words, of the potential demand at node {ui}, {uh/uTW}
      calls will still be lost, {uh*upw/uTW} will be served by phone
      {pi}, and {uh*uqw/uTW} will be served by {qi}, where {uTW =
      uNw + upw + uqw}. Note that this model predicts that the new
      phone {qi} captures some of the calls that used to be lost, and
      also some of the calls that used to be satisfied by {pi}.
      Specifically, the usage of latter decreases from {uh*upw/(1 +
      upw)} to {uh*upw/(1 + upw + uqw)}.
      
      The same assumption will be used when there are three or more
      phones {p,q,r...}: the potential demand {uh} is assumed to be
      partitioned between `no-phone' and those phones proportionally
      to their weights {1, upw, uwq, urw...}.
      
      At present, we assume that the limiting capture fraction {upf}
      depends only on the distance {upd}, i.e. {upf = f(dist(ui,pi))},
      where {f(d)} is some function that decays to 0 as the distance
      {d} goes to infinity. Specifically, we assume that {f} is a
      Gaussian {g(d) = K*exp(-d^2/S2)}. The constant {S2} is such that
      the capture fraction is reduced by 50% as {d} increases from 0
      to {partDist}. The constant {K} (striclty less than 1.0)
      represents the faction of potential demand at {pi} that is
      captured by a phone located at {pi} itself. */
      
    
    /* Quotient for Gaussian exponent: */
    double S2 = partDist*partDist/log(2.0);
    /* Path-cost at which the capture factor is negligible: */
    double tinyFrac = 0.005;
    float stopDist = (float)sqrt(-S2*log(tinyFrac));

    st_user_phone_matrix(m, maxDist, ph, nph, &H, &nH);

    /* Pass 2: for each phone, collect served demand of nearby vertices: */


            double upf; /* Limiting capture fraction of {pi} at {ui}. */
            if (ui == pi)
              { upf = K; }
            else
              { double upd2 = upd*upd;
                double g = exp(-upd2/S2);
                upf = K*g;
              }
            double upw = upf/(1.0 - upf);
            totw[ui] += upw;

    st_map_init_costs(m, d, e, c);
    int j;
    for(j = 0; j < nph; j++)
      { Phone *pph = &(ph[j]);
        int pi = pph->vertex;
        pph->usage = 0.0;

        
        auto bool_t collect_demand(int ui, double upd);
        
        bool_t collect_demand(int ui, double upd)
          { if (d > stopDist)
              { return TRUE; }
            else
              { double upf; /* Limiting capture fraction of {pi} at {ui}. */
                if (ui == pi)
                  { upf = K; }
                else
                  { double upd2 = upd*upd;
                    double g = exp(-upd2/S2);
                    upf = K*g;
                  }
                double upw = upf/(1.0 - upf);
              }
          }


        Point *ppos = &(pd->p);
        /* Find vertices served by {pph}, add a fraction of their demands: */
        double usage = 0.0;
        st_map_compute_costs(m, pi, stopDist, r, &nr, d, e, c);
        int k;
        for (k = 0; k < nr; k++)
          { int ui = r[k];
            double upf; /* Limiting capture fraction of {pi} at {ui}. */
            if (ui == pi)
              { upf = K; }
            else
              { VertexData *ud = m->vd[ui];
                Point *upos = &(ud->p);
                double dx = ppos->c[0] - upos->c[0];
                double dy = ppos->c[1] - upos->c[1];
                double upd2 = dx*dx + dy*dy + 0.0001;
                double g = exp(-upd2/S2);
                upf = K*g;
              }
            double upw = upf/(1.0 - upf);
            totw[ui] += upw;
          }
          
            usage += frac*dem[ui];
        pph->usage = usage;
        st_map_reset_costs(m, r, nr, d, e, c);
      }

    /* Compute unsatisfied demand at each vertex: */
    for (ui = 0; ui < m->nv; ui++) { lost[ui] = dem[ui]/(1.0 + totw[ui]); }
  }

voif st_enum_served_vertices(Map *m, double partDist, int pi)
  { VertexData *pd = m->vd[pi];
    Point *ppos = &(pd->p);
    /* Find vertices {ui} served by {pph}, increment their {totw[ui]}: */
    double usage = 0.0;
    st_map_compute_costs(m, pi, stopDist, r, &nr, d, e, c);
    int k;
    for (k = 0; k < nr; k++)
      { int ui = r[k];
        visit(ui, d...)
      }
    st_map_reset_costs(m, r, nr, d, e, c);
  }
