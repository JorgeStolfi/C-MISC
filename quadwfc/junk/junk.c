/* Last edited on 2012-07-22 13:49:10 by stolfilocal */ 

//////////////////////////////////////////////////////////////////////

// From quadwfc.c - argument parsing

double parse_double(int *argn, int argc, char **argv, double lo, double hi)
  { char *rest;
    double x;
    if (*argn >= argc) 
      { arg_error("argument value is missing", (*argn)-1, argv[(*argn)-1]); }
    x = strtod(argv[*argn], &rest);
    if (*rest != '\0') 
      { arg_error("invalid number", *argn, argv[*argn]); } 
    if ((x < lo) || (x > hi)) 
      { arg_error("out of range", *argn, argv[*argn]); }
    ++(*argn);
    return x;
  }

void arg_error(char *msg, int argn, char *arg)
  {
    fprintf(stderr, "%s: **", PROG_NAME);
    if (argn >= 0) 
      { fprintf(stderr, " argv[%d]", argn); }
    if (arg != NULL) 
      { fprintf(stderr, " = %s", arg); }
    fprintf(stderr, " %s\n", msg);
    fprintf(stderr, "usage: %s \\\n%s", PROG_NAME, PROG_USAGE);
    exit(1);
  }     

//////////////////////////////////////////////////////////////////////
 
void generic_plot
  ( PSStream *fps,
    hr3_pmap_t *map,             /* Perspective projection matrix. */
    mesh_t *tri,                 /* Mesh to plot. */
    int N,                       /* Mesh refinement parameter. */
    triangle_visit_proc_t *proc  /* Called for each triangle fragment. */
  )
  { 
    auto void processTriangle(r3_t *P, r3_t *Q, r3_t *R);
      /* Subdivides the given triangle into {N^2} trianglets and calls {proc}
        on them. */

    r3_t S[N + 1]; /* Saved points. */
    
    void processTriangle(r3_t *P, r3_t *Q, r3_t *R)
      { double fN = (double)N;
        int i;
        for (i = 0; i <= N; i++)
          { int j;
            for (j = 0; j <= N-i; j++) 
              { int k = N - i - j;
                r3_t T = (r3_t)
                  {{(P->c[0]*i + Q->c[0]*j + R->c[0]*k)/fN,
                    (P->c[1]*i + Q->c[1]*j + R->c[1]*k)/fN,
                    (P->c[2]*i + Q->c[2]*j + R->c[2]*k)/fN
                  }};
                r3_dir(&T, &T); 
                if (i > 0)
                  { /* Plot triangles using points from previous row: */
                    /* Be sure to preserve orientation rel. to {P,Q,R}: */
                    r3_t *V = &(S[j]);  
                    r3_t *W = &(S[j+1]);
                    if (j > 0)
                      { r3_t *U = &(S[j-1]);
                        proc(V, U, &T);
                      }
                    proc(W, V, &T);
                  }
                /* Save point: */
                S[j] = T;
              }
          }
      }

    /* Process each part of the given triangulation */
    int NT = tri->side.nel;
    double epsilon = 0.0001; /* Relative perturbation to avoid edges */
    int i, j;
    /* Plot each triangle: */
    for (i = 0; i < NT; i++)
      { qarc_t e = tri->side.el[i];
        r3_t P = DEST(ONEXT(e))->curr;
        r3_t Q = ORG(e)->curr;
        r3_t R = DEST(e)->curr;

        /* Perturb points so that they lie slightly inside the triangle: */
        for (j = 0; j < 3; j++)
          { double Bj = (P.c[j] + Q.c[j] + R.c[j])/3.0;
            P.c[j] = (1.0-epsilon)*P.c[j] + epsilon*Bj;
            Q.c[j] = (1.0-epsilon)*Q.c[j] + epsilon*Bj;
            R.c[j] = (1.0-epsilon)*R.c[j] + epsilon*Bj;
          }
        r3_dir(&P, &P);
        r3_dir(&Q, &Q); 
        r3_dir(&R, &R);

        processTriangle(&P, &Q, &R);
      }
  }


//////////////////////////////////////////////////////////////////////
// instrumented version of {rec_delaunay} from triangulate.c:

void rec_delaunay(
    sref_vec_t *st,
    int sl, int sh,
    qarc_t *le, 
    qarc_t *re,
    int level
  )
  {
    auto void debug_arc(char *label, qarc_t e);
    
    bool_t debug = (sl == 39) && (sh == 44);
    { int i; 
      /* Set vertex numbers by order: */
      for (i = sl; i <= sh; i++) { st->el[i]->num = i; }
    }
      
    
    fprintf(stderr, "%*s  triangulating [%d..%d]", 2*level, "", sl, sh);
    fprintf(stderr, "  X = [%5g _ %5g]", st->el[sl]->curr.c[0], st->el[sh]->curr.c[0]);
    fprintf(stderr, "\n");
    
    if (debug)
      { int k;
        for (k = sl; k <= sh; k++)
          { sample_t *s = st->el[k];
            fprintf(stderr, "%*s    s[%03d] =", 2*level, "", k);
            fprintf(stderr, " ( %.5g %.5g )", s->curr.c[0], s->curr.c[1]);
            fprintf(stderr, "\n");
          }
      }
    if (sh == sl+1) 
      {
	/* Only two samples. */
        qarc_t a = quad_make_edge();
	SET_ORG(a, st->el[sl]); 
        SET_DEST(a, st->el[sl+1]);
	*le = a; *re = SYM(a);
      }
    else if (sh == sl+2) 
      {
	/* Only three samples. */
        qarc_t a = quad_make_edge();
	qarc_t b = quad_make_edge();
        sample_t *u = st->el[sl];
	sample_t *v = st->el[sl+1];
	sample_t *w = st->el[sl+2];
	double ct = orient(&(u->curr), &(v->curr), &(w->curr));
	quad_splice(SYM(a), b);
	SET_ORG(a, u); SET_DEST(a, v);
	SET_ORG(b, v); SET_DEST(b, w);
	if (ct == 0) 
	  { *le = a; *re = SYM(b); }
	else 
	  { qarc_t c = connect(b, a);
	    if (ct > 0) 
	      { *le = a; *re = SYM(b); }
	    else 
	      { *le = SYM(c); *re = c; }
	  }
      }
    else
      {
	qarc_t ldo, ldi, rdi, rdo;
	qarc_t basel, lcand, rcand;

        int sm = (sl+sh)/2;

        rec_delaunay(st, sl, sm, &ldo, &ldi, level+1);
	rec_delaunay(st, sm+1, sh, &rdi, &rdo, level+1);

	while (1) 
          {
	    if (leftof(&ORGP(rdi), ldi))
              { ldi = LNEXT(ldi); }
	    else if (rightof(&ORGP(ldi), rdi)) 
              { rdi = ONEXT(SYM(rdi)); }
	    else break;
	  }

	basel = connect(SYM(rdi), ldi);
	if (ORG(ldi) == ORG(ldo)) { ldo = SYM(basel); }
	if (ORG(rdi) == ORG(rdo)) { rdo = basel; }
        
	while (1) 
          {
	    if (debug) { debug_arc("basel =", basel); }
            
            lcand = ONEXT(SYM(basel));
	    if (rightof(&DESTP(lcand), basel))
              { 
                while (incircle(&DESTP(basel), &ORGP(basel), &DESTP(lcand), &DESTP(ONEXT(lcand)))) 
                  { qarc_t t = ONEXT(lcand); 
                    if (debug) { debug_arc("killed left", lcand); }
                    quad_destroy_edge(lcand); 
                    lcand = t;
                  }
              }

	    rcand = OPREV(basel);
	    if (rightof(&DESTP(rcand), basel))
	      { while (incircle(&DESTP(basel), &ORGP(basel), &DESTP(rcand), &DESTP(OPREV(rcand)))) 
                  { qarc_t t = OPREV(rcand); 
                    if (debug) { debug_arc("killed right", rcand); }
                    quad_destroy_edge(rcand);
                    rcand = t;
                  }
              }

	    if (!rightof(&DESTP(lcand), basel) && !rightof(&DESTP(rcand), basel)) break;

	    if ( !rightof(&DESTP(lcand), basel) ||
		 ( rightof(&DESTP(rcand), basel) && 
                   incircle(&DESTP(lcand), &ORGP(lcand), &ORGP(rcand), &DESTP(rcand))
                 )
               )
	      { basel = connect(rcand, SYM(basel)); }
	    else
	      { basel = connect(SYM(basel), SYM(lcand)); }
	  }
	*le = ldo; *re = rdo;
      }
      
    void debug_arc(char *label, qarc_t e)
      { 
        sample_t *u = ORG(e);
        sample_t *v = DEST(e);
        fprintf(stderr, "%*s    %s (%d,%d)\n", 2*level, "", label, u->num, v->num);
      }


    /* Paranoia: */
    create_face_records(*le);
    mesh_t *m = mesh_from_topology(*le);
    fprintf(stderr, "%*s  ns = %d", 2*level, "", sh-sl+1);
    fprintf(stderr, " nv = %d ne = %d nf = %d", m->out.nel, m->arc.nel/2, m->side.nel);
    fprintf(stderr, "\n");
    affirm(sh-sl+1 == m->out.nel, "lost vertices");
    { int i; 
      /* Recycle face records: */
      for (i = 0; i < m->side.nel; i++) { free(LEFT(m->side.el[i])); }
      /* Clean up face pointers: */
      for (i = 0; i < m->arc.nel; i++) { SET_LEFT(m->arc.el[i],NULL); }
      /* Restore vertex numbers: */
      for (i = sl; i <= sh; i++) { st->el[i]->num = i; }
      /* Free mesh storage: */
      free(m->arc.el); free(m->out.el); free(m->side.el);
      free(m);
    }
  }


//////////////////////////////////////////////////////////////////////
proj_map_t compute_persp_map(hr3_t *obs, r3_t *ctr, double rad, hr3_t *upp)
  { proj_map_t m;
    r3_t r, s, t, u, v;
    r4x4_t mt, mr, mc, mrt;
    int i, j;
    
    affirm(obs->h.c[0] >= 0.0, "bad obs"); /* Observer must be hither or infinite */
    affirm(upp->h.c[0] >= 0.0, "bad upp"); /* Zenith must be hither or infinite */

    /* Start with a translation from {ctr} to the origin: */
    for (i = 0; i < 4; i++){
      for (j = 0; j < 4; j++){
        if (i == j)
          { mt.c[i][j] = 1.0; }
        else if (i == 0)
          { mt.c[i][j] = -(ctr->c[j]); }
        else
          { mt.c[i][j] = 0.0; }
      }
    }
    
    /* Append the rotation matrix that moves {obs} to the Z-axis: */
    t = r3_hr3_dir(ctr, obs);
    u = r3_hr3_dir(ctr, upp);
    r3_decomp(&u, &t, &v, &s);
    /* Zenith reference point {upp} must not be on imagesys Z axis: */
    affirm(r3_norm_sqr(&s) >= 1.0e-12, "bad zenith"); 
    r3_dir(&s, &s);
    r3_cross(&s, &t, &r);
    r3_dir(&r, &r);

    mr.c[0][0] = 1.0;
    for (i = 1; i < 4; i++)
      { mr.c[0][i] = 0.0;
        mr.c[i][0] = 0.0;
        mr.c[i][1] = r.c[i-1];
        mr.c[i][2] = s.c[i-1];
        mr.c[i][3] = t.c[i-1];
      }
    
    /* Compose the two matrices: */
    r4x4_mul(&mt, &mr, &mrt);

    /* Do we need a conical projection step? */
    if (obs->h.c[0] == 0.0)
      { /* Observer is at infinity; cilindrical projection. */
        m.dir = mrt;
      }
    else
      { /* Observer is finite; add conical projection step. */
        double d = r3_hr3_dist(ctr, obs);
        double uno = -1.0;
        if (fabs(d) > 1.0)
          { uno = -1.0/d; d = 1.0; }
        for (i = 0; i < 4; i++)
          { for (j = 0; j < 4; j++)
              { if (i == j)
                  { mc.c[i][j] = d; }
                else
                  { mc.c[i][j] = 0.0; }
              }
          }
        mc.c[3][0] = uno;
        r4x4_mul(&mrt, &mc, &(m.dir));
      }

    /* Now combine {map} with a uniform scale of {1/rad}: */
    for (i = 0; i < 4; i++)
      for (j = 1; j < 4; j++)
        { m.dir.c[i][j] /= rad; }

    /* Compute inverse matrix: */
    r4x4_inv(&m.dir, &m.inv);
    return m;
  }

----------------------------------------------------------------------

typedef struct hr3_t  /* Point of R^3 in homogeneous coordinates. */
  { r4_t h;  /* {h.c[0..3]} are the homogeneous coords {[w,x,y,z]}. */
  } hr3_t;

typedef struct proj_map_t
  { r4x4_t dir; 
    r4x4_t inv;
  } proj_map_t;  
  /* {dir} is the map's matrix, {inv} is its inverse. */

    
r3_t r3_hr3_dir(r3_t *frm, hr3_t *tto);
  /* The unit vector pointing from {frm} towards {tto}. */

double r3_hr3_dist(r3_t *a, hr3_t *b);
  /* Euclidean distance from {a} to {b}, possibly {INF}. */

r3_t r3_hr3_dir(r3_t *frm, hr3_point_t *tto)
  { double fw = 1.0;
    double tw = tto->h.c[0];

    double fx = frm->c[1];
    double tx = tto->h.c[1];
    double dx = fw * tx - tw * fx;

    double fy = frm->c[2];
    double ty = tto->h.c[2];
    double dy = fw * ty - tw * fy;

    double fz = frm->c[3];
    double tz = tto->h.c[3];
    double dz = fw * tz - tw * fz;

    double length = hypot(hypot(dx, dy), dz);

    return (r3_t){{dx/length, dy/length, dz/length}};
  }

double r3_hr3_dist(r3_t *a, hr3_point_t *b)
  { double aw = 1.0;
    double bw = 1.0/b->h.c[0];

    double ax = a->c[1];
    double bx = b->h.c[1];
    double dx = ax*aw - bx*bw;

    double ay = a->c[2];
    double by = b->h.c[2];
    double dy = ay*aw - by*bw;

    double az = a->c[3];
    double bz = b->h.c[3];
    double dz = az*aw - bz*bw;

    return hypot(hypot(dx, dy), dz);
  }

typedef r3_t ProjPoint;  /* Persp-mapped point {(xplot,yplot,height)}. */

hr3_t hr3_from_r3(r3_t *c);
r3_t r3_from_hr3(hr3_t *p);
hr3_t hr3_map_point(hr3_t *p, proj_map_t *m);

hr3_t hr3_from_r3(r3_t *c)
  { return (hr3_t){(r4_t){{1.0, c->c[0], c->c[1], c->c[2]}}}; }
    
r3_t r3_from_hr3(hr3_t *p)
  { double w = p->h.c[0];
    affirm(w != 0.0, "null weight");
    return (r3_t){{p->h.c[1]/w, p->h.c[2]/w, p->h.c[3]/w}};
  }

hr3_t hr3_map_point(hr3_t *p, proj_map_t *m)
  { hr3_t q;
    r4x4_map_row(&(p->h), &(m->dir), &(q.h));
    return q;
  }

hr3_t hr3_inv_map_point(hr3_t *p, proj_map_t *m)
  { hr3_t q;
    r4x4_map_row(&(p->h), &(m->inv), &(q.h));
    return q;
  } 

hr3_point_t hr3_map_point(hr3_point_t *p, proj_map_t *m);
  /* Applies projective map {m} to point {p}. */

hr3_point_t hr3_inv_map_point(hr3_point_t *p, proj_map_t *m);
  /* Applies the inverse of projective map {m} to point {p}. */
  
//////////////////////////////////////////////////////////////////////

qarc_vec_t renumber_edges(qarc_vec_t a, unsigned int start);
  /* Enumerates undirected edges reachable from the arcs {a.e[..]}, and
    stores in their {num} fields consecutive integers beginning with {start}.
    Returns an array {res} where {res.e[i]} is one primally reachable
    arc from the edge with number {start+i}.

    An arc {b} is /primally reachable/ from an arc {a} iff it can be obtained 
    from {a} by a finite number of {Sym} and {Onext} operations (no {Rot}s
    or {Tor}s).  An edge {e} is /reachable/ iff some of its arcs is reachable. */

qarc_vec_t renumber_edges(qarc_vec_t a, unsigned int start)
  {
    auto void visit_edge(qarc_t a);
      /* If {Edge(a)} has not been visited before, renumber it. */

    auto bool_t edge_is_visited(qarc_t a);
      /* TRUE if {Edge(a)} has been visited before. */
    
    #define GUESS_NVISITED 1024
    qarc_vec_t visited = qarc_vec_new(GUESS_NVISITED);
    int nVisited = 0;
      /* {visited.e[0..nVisited-1]} contain one arc from each visited
        edge. An edge record {e} has been visited if
        {Edge(visited.e[e->num - start]) == e}. */

    bool_t edge_is_visited(qarc_t a)
      { unsigned int en = EDGE(a)->mark;
        return 
          (en >= start) && 
          (en < start + nVisited) && 
          (EDGE(visited.e[en-start]) == EDGE(a));
      }

    void visit_edge(qarc_t a)
      { if (! edge_is_visited(a))
          { /* Edge {Edge(a)} hasn't been visited yet. */
            qarc_vec_expand(&visited, nVisited);
            EDGE(a)->mark = start + nVisited;
            visited.e[nVisited] = a;
            nVisited++;
          }
      }

    int i;
    int nClosed = 0; 
      /* Edges {EDGE(visited.e[0..nClosed-1])} are the edges whose
        children were already visited.  The children of {e}, by 
        definition, are {Onext(e)} and {Onext(Sym(e))}. */

    /* Put the given arcs on the visited (minus repetitions): */
    for (i = 0; i < a.ne; i++) { visit_edge(a.e[i]); }
    /* Visit descendants of visited edges in BFS order: */
    while (nClosed < nVisited )
      { qarc_t s = visited.e[nClosed];
        visit_edge(quad_onext(s));
        visit_edge(quad_onext(quad_sym(s)));
        nClosed++;
      }
    /* Return the edges which were found: */
    qarc_vec_trim(&visited, nVisited);
    /* fprintf(stderr, "visited %d edges\n", nVisited); */
    return visited;
  }

