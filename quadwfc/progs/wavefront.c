/* See {wavefront.h}. */
/* Last edited on 2023-02-12 23:53:42 by stolfi */

#include <basic.h>
#include <wavefront.h>

/* TRAVERSAL */

void kill_sample(sample_t *s)
  {
    s->sgn = NULL;
    s->pos = (r3_t){{INF,INF,INF}};
    s->vel = (r3_t){{0,0,0}};
  }

bool_t sample_is_dead(sample_t *s)
  {  return ((s->sgn == NULL) || (*(s->sgn) == '\000'));  }

quad_arc_vec_t renumber_edges(quad_arc_vec_t a)
  {
    quad_arc_vec_t o = quad_arc_vec_new(a.ne);
    
    quad_edge_num_t nE = quad_renumber_edges(&a, &o);
    assert(nE == o.ne);
    
    return o;
  }

quad_arc_vec_t arcs_from_edges(quad_arc_vec_t edge)
  { quad_arc_vec_t arc = quad_arc_vec_new(2*edge.ne);
    int i;
    for (i = 0; i < edge.ne; i++)
      { arc.e[2*i] = edge.e[i];
        arc.e[2*i+1] = quad_sym(edge.e[i]);
      }
    return arc;
  }
  
quad_arc_vec_t renumber_vertices(quad_arc_vec_t arc)
  { 
    auto void visit_org(quad_arc_t a);
      /* Visits vertex {quad_org(a)}. */

    auto bool_t org_is_visited(quad_arc_t a);
      /* TRUE if {quad_org(a)} has been visited already. */
   
    quad_arc_vec_t r = quad_arc_vec_new(1);
    int nVertices = 0;
    int i;
   
    for (i = 0; i < arc.ne; i++)
      { if (! org_is_visited(arc.e[i])) { visit_org(arc.e[i]); } }

    quad_arc_vec_trim(&(r), nVertices);
    /* fprintf(stderr, "visited %d vertices\n", nVertices); */
    return r;

    bool_t org_is_visited(quad_arc_t a)
      { segment_t *v = quad_org(a);
        affirm(v != NULL, "null vertex");
        return (v->num < nVertices) && (quad_org(r.e[v->num]) == v);
      }

    void visit_org(quad_arc_t a)
      { segment_t *v = quad_org(a);
        affirm(v != NULL, "null vertex");
        quad_arc_t b = a;
        v->num = nVertices;
        do 
          { affirm(quad_org(b) == v, "inconsistent ORG");
            b = quad_onext(b);
          }
        while (b != a);
        quad_arc_vec_expand(&(r), nVertices);
        r.e[nVertices] = a;
        nVertices++;
      }
  }
  
quad_arc_vec_t renumber_faces(quad_arc_vec_t arc)
  { auto void visit_left(quad_arc_t a);
      /* Visits face {Left(a)}. */

    auto bool_t left_is_visited(quad_arc_t a);
      /* TRUE if {Left(a)} has been visited already. */
   
    quad_arc_vec_t r = quad_arc_vec_new(1);
    int nFaces = 0;
    int i;
   
    bool_t left_is_visited(quad_arc_t a)
      { face_t *f = LEFT(a);
        affirm(f != NULL, "null face");
        return (f->num < nFaces) && (LEFT(r.e[f->num]) == f);
      }

    void visit_left(quad_arc_t a)
      { face_t *f = LEFT(a);
        affirm(f != NULL, "null face");
        quad_arc_t b = a;
        f->num = nFaces;
        do 
          { affirm(LEFT(b) == f, "inconsistent Left");
            b = quad_lnext(b);
          }
        while (b != a);
        quad_arc_vec_expand(&(r), nFaces);
        r.e[nFaces] = a;
        nFaces++;
      }

    for (i = 0; i < arc.ne; i++)
      { if (! left_is_visited(arc.e[i])) { visit_left(arc.e[i]); } }

    quad_arc_vec_trim(&(r), nFaces);
    /* fprintf(stderr, "visited %d faces\n", nFaces); */
    return r;
  }
  
#define INIT_VISITED_SIZE 1024

void create_face_records(quad_arc_t e)
  {
    auto bool_t left_is_visited(quad_arc_t a);
      /* TRUE if face {LEFT(a)} has been visited before. */

    auto void visit_left(quad_arc_t a);
      /* If the left face of {a} has not been visited yet, 
        create its face record, set all pointers to it,
        and mark that face as visited. */

    quad_arc_vec_t side = quad_arc_vec_new(INIT_VISITED_SIZE);
    int nFaces = 0; 
      /* A face has been visited if all arcs around it have the {Left}
        field pointing to the same triangle record {f}, and
        {LEFT(side[f.num]) == f}. Outside of {visit_left}, either all
        arcs around the face satisfy this condition, or none of them
        does. */

    void visit_left(quad_arc_t a)
      { 
        if (! left_is_visited(a))
          { 
            face_t *f = (face_t *)malloc(sizeof(face_t));
            quad_arc_t b = a;
            int order = 0;
            f->num = nFaces;
            f->omit = FALSE;
            do 
              { SET_LEFT(b, f);
                b = quad_lnext(b);
                order++;
              }
            while (b != a);
            quad_arc_vec_expand(&side, nFaces);
            side.e[nFaces] = a;
            nFaces++;
          }
      }

    bool_t left_is_visited(quad_arc_t a)
      { face_t *f = LEFT(a);
        return 
          (f != NULL) && 
          (f->num < nFaces) && 
          (LEFT(side.e[f->num]) == f);
      }

    int nClosed = 0; 
      /* Faces {LEFT(side.e[0..nClosed-1])} are the visited faces whose 
        children were already visited.  The children of {LEFT(e)}, by 
        definition, are all faces that share a face or edge with 
        {LEFT(e)}. */

    /* Put the given face on the visited (minus repetitions): */
    visit_left(e);
      
    /* Visit descendants of visited faces in BFS order: */
    while (nClosed < nFaces )
      { quad_arc_t s = side.e[nClosed];
        /* Enumerate edges around {LEFT(s)} including {s}: */
        quad_arc_t e = s;
        do
          { /* Enumerate edges around {Dest(e)} except {e}: */
            quad_arc_t a = quad_dnext(e);
            while (a != e) { visit_left(a); a = quad_dnext(a); }
            e = quad_lnext(e);
          }
        while (e != s);
        nClosed = nClosed + 1;
      }
    /* Be nice, recycle: */
    free(side.e);
  }
  
/* VECTORS OF MESH ELEMENTS */
  
vec_typeimpl(sref_vec_t,sref_vec,sref_t);

vec_typeimpl(fref_vec_t,fref_vec,fref_t);

