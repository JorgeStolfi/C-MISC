/* See {mesh.h}. */
/* Last edited on 2013-10-02 03:15:47 by stolfilocal */

#include <mesh.h>
#include <basic.h>
#include <wavefront.h>

#include <r3.h>
#include <filefmt.h>
#include <fget.h>
#include <nget.h>

#define FILE_VERSION "2005-07-09"

void write_mesh(FILE *wr, mesh_t *tri)
  { 
    auto void write_point(r3_t *p);
    auto void write_arc(qarc_t a);
    auto void write_org(qarc_t a);
    auto void write_left(qarc_t a);
    auto void write_edge_data(qarc_t e);

    int NV = tri->out.ne;
    int NE = tri->arc.ne / 2;
    int NF = tri->side.ne;

    void write_point(r3_t *p)
      { int i;
        for (i = 0; i <= 2; i++)
          { if (i != 0){ fputc(' ', wr); }
            fprintf(wr, "%.16g", p->c[i]);
          }
      }

    void write_arc(qarc_t e)
      { fprintf(wr, "%d:%d", quad_edge(e)->mark, quad_tumble_code(e)); }

    void write_org(qarc_t a)
      { fputc('v', wr); fprintf(wr, "%d", quad_org(a)->num); }

    void write_left(qarc_t a)
      { fputc('f', wr); fprintf(wr, "%d", LEFT(a)->num); }

    void write_edge_data(qarc_t e)
      { write_arc(e); fputc(' ', wr);
        fputc(' ', wr);

        write_arc(quad_onext(e)); fputc(' ', wr);
        write_arc(quad_onext(quad_sym(e))); fputc(' ', wr);
        fputc(' ', wr);

        write_org(e); fputc(' ', wr);
        write_org(quad_sym(e)); fputc(' ', wr);
        write_left(e); fputc(' ', wr);
        write_left(quad_sym(e)); 
        fputc('\n', wr);
      }

    int i;

    filefmt_write_header(wr, "wfmesh", FILE_VERSION);

    fprintf(wr, "vertices = %d\n", NV);
    fprintf(wr, "edges = %d\n", NE);
    fprintf(wr, "faces = %d\n", NF);

    fprintf(wr, "vertices:\n");
    for (i = 0; i < NV; i++)
      { qarc_t a = tri->out.e[i];
        affirm(a != NULL_REF, "vertex with no arcs");
        segment_t *v = quad_org(a);
        affirm(v->num == i, "inconsistent vertex num");
        fprintf(wr, "%d ", i);
        write_arc(a); fputc(' ', wr);
        if (sample_is_dead(&(v->curr)))
          { fputc('*', wr); }
        else
          { write_point(&(v->curr.pos)); fputc(' ', wr);
            write_point(&(v->curr.vel)); 
          }
        fputc('\n', wr);
      }

    fprintf(wr, "edges:\n");
    for (i = 0; i < tri->arc.ne; i += 2)
      { qarc_t a = tri->arc.e[i];
        affirm(a != NULL_REF, "edge with no arcs");
        quad_edge_rec *e = quad_edge(a);
        affirm(e->mark == i/2, "inconsistent edge num");
        fprintf(wr, "%d ", i/2);
        write_edge_data(a);
      }

    fprintf(wr, "faces:\n");
    for (i = 0; i < tri->side.ne; i++)
      { qarc_t a = tri->side.e[i];
        affirm(a != NULL_REF, "face with no arcs");
        face_t *f = LEFT(a);
        affirm(f->num == i, "inconsistent face num");
        fprintf(wr, "%d ", i);
        write_arc(a); fputc(' ', wr);
        fputc((f->omit ? 'T' : 'F'), wr);
        fputc('\n', wr);
      }

    filefmt_write_footer(wr, "wfmesh");
    fflush(wr);
  }

mesh_t *read_mesh(FILE *rd)
  { 
    mesh_t *tri = (mesh_t*)notnull(malloc(sizeof(mesh_t)), "out of mem");
    filefmt_read_header(rd, "wfmesh", FILE_VERSION);
    int NV = nget_int(rd, "vertices"); fget_eol(rd);
    affirm(NV >= 0, "bad NV");
    int NE = nget_int(rd, "edges"); fget_eol(rd);
    affirm(NE >= 0, "bad NE");
    int NF = nget_int(rd, "faces"); fget_eol(rd);
    affirm(NF >= 0, "bad NF");
    tri->out = qarc_vec_new(NV);
    tri->arc = qarc_vec_new(2*NE);
    tri->side = qarc_vec_new(NF);
    {
      auto qarc_t read_arc(void);
      auto segment_t *read_vertex(void);
      auto face_t *read_face(void);
      auto void read_point(r3_t *c);
     
      sref_vec_t site = sref_vec_new(NV);
      qarc_vec_t edge = qarc_vec_new(NE);
      fref_vec_t face = fref_vec_new(NF);
     
      qarc_t read_arc(void)
        { int qnum, rnum; 
          qarc_t a;
          fget_skip_spaces(rd);
          qnum = fget_int(rd);
          affirm((qnum >= 0) && (qnum < NE), "bad quad number");
          fget_match(rd, ":");
          rnum = fget_int(rd);
          affirm((rnum >= 0) & (rnum <= 3), "inconsistent rot number");
          a = edge.e[qnum];
          while (quad_tumble_code(a) != rnum) { a = quad_rot(a); }
          return a;
        }

      segment_t *read_vertex(void)
        { int vnum;
          fget_skip_spaces(rd);
          fget_match(rd, "v");
          vnum = fget_int(rd);
          affirm((vnum >= 0) && (vnum < NV), "bad vertex number");
          return site.e[vnum];
        }

      face_t *read_face(void)
        { int fnum;
          fget_skip_spaces(rd);
          fget_match(rd, "f");
          fnum = fget_int(rd);
          affirm((fnum >= 0) && (fnum < NF), "bad vertex number");
          return face.e[fnum];
        }

      void read_point(r3_t *c)
        { int i;
          for (i = 0; i <= 2; i++)
            { fget_skip_spaces(rd); 
              c->c[i] = fget_double(rd);
            }
        }

      /* Create all vertices, edges, and faces: */
      int i;
      for (i = 0; i < NV; i++)
        { segment_t *v = (segment_t*)notnull(malloc(sizeof(segment_t)), "out of mem");
          site.e[i] = v; 
        }
      for (i = 0; i < NE; i++)
        { qarc_t a = quad_make_edge();
          edge.e[i] = a;
        }
      for (i = 0; i < NF; i++)
        { face_t *f = (face_t *)notnull(malloc(sizeof(face_t)), "out of mem");
          face.e[i] = f; 
        }

      /* Read vertex coordinates: */
      fget_match(rd, "vertices:"); fget_eol(rd);
      for (i = 0; i < NV; i++)
        { segment_t *v = site.e[i];
          int j = fget_int(rd);
          affirm(j == i, "inconsistent vertex number");
          fget_skip_spaces(rd);
          tri->out.e[i] = read_arc();
          fget_skip_spaces(rd);
          v->num = i;
          if (fget_test_char(rd, '*'))
            { kill_sample(&(v->curr)); }
          else
            { read_point(&(v->curr.pos));
              read_point(&(v->curr.vel));
            }
          fget_eol(rd);
        }

      /* Read edge pointers: */
      fget_match(rd, "edges:"); fget_eol(rd);
      for (i = 0; i < NE; i++)
        { qarc_t a = edge.e[i];
          int j = fget_int(rd);
          affirm(j == i, "inconsistent edge number");
          quad_edge(a)->mark = i;
          fget_skip_spaces(rd);
          { qarc_t a = read_arc();
            qarc_t b = read_arc(); 
            qarc_t c = read_arc();
            segment_t *o = read_vertex();
            segment_t *d = read_vertex();
            face_t *l = read_face();
            face_t *r = read_face();
            tri->arc.e[2*i] = a; 
            tri->arc.e[2*i+1] = quad_sym(a);
            quad_splice(a, quad_oprev(b));
            if (b != quad_sym(a)){ quad_splice(quad_sym(a), quad_oprev(c)); }
            quad_odata(a) = o;
            quad_ddata(a) = d;
            quad_ldata(a) = l;
            quad_rdata(a) = r;
          }
          fget_eol(rd);
        }

      /* Read face data: */
      fget_match(rd, "faces:"); fget_eol(rd);
      for (i = 0; i < NF; i++)
        { face_t *f = face.e[i];
          int j = fget_int(rd);
          affirm(j == i, "inconsistent face number");
          fget_skip_spaces(rd);
          tri->side.e[i] = read_arc();
          f->num = i; 
          f->omit = fget_bool(rd);
          fget_eol(rd);
        }
    }
    filefmt_read_footer(rd, "wfmesh");
    /* compute_face_geometry(tri); */
    return tri;
  }

mesh_t *mesh_from_topology(qarc_t e)
  { mesh_t *tri = (mesh_t *)notnull(malloc(sizeof(mesh_t)), "no mem for mesh");
    qarc_vec_t edge = renumber_edges((qarc_vec_t){1, &e}, 0);
    tri->arc = arcs_from_edges(edge);
    tri->out = renumber_vertices(tri->arc);
    tri->side = renumber_faces(tri->arc);
    /* compute_face_geometry(tri); */
    return tri;
  }

void enum_trianglets(mesh_t *tri, int N, double epsilon, triangle_visit_proc_t *proc)
  { 
    r3_t S[N + 1]; /* Saved points. */

    auto void process_triangle(r3_t *P, r3_t *Q, r3_t *R);
      /* Subdivides the given triangle into {N^2} trianglets and calls {proc}
        on them. */

    /* Process each part of the given triangulation */
    int NT = tri->side.ne;
    int k;
    /* Pov each triangle: */
    for (k = 0; k < NT; k++)
      { qarc_t e = tri->side.e[k];
        face_t *f = LEFT(e);
        if (! f->omit)
          { r3_t P = quad_dst(quad_onext(e))->curr.pos;
            r3_t Q = quad_org(e)->curr.pos;
            r3_t R = quad_dst(e)->curr.pos;
            if (epsilon != 0.0)
              { /* Perturb points so that they lie slightly inside the triangle: */
                int j;
                for (j = 0; j < 3; j++)
                  { double Bj = (P.c[j] + Q.c[j] + R.c[j])/3.0;
                    P.c[j] = (1.0-epsilon)*P.c[j] + epsilon*Bj;
                    Q.c[j] = (1.0-epsilon)*Q.c[j] + epsilon*Bj;
                    R.c[j] = (1.0-epsilon)*R.c[j] + epsilon*Bj;
                  }
              }
            process_triangle(&P, &Q, &R);
          }
      }
    
    void process_triangle(r3_t *P, r3_t *Q, r3_t *R)
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
                if (i > 0)
                  { /* Visit triangles using points from previous row: */
                    /* Be sure to preserve orientation rel. to {P,Q,R}: */
                    r3_t *V = &(S[j]);  
                    r3_t *W = &(S[j+1]);
                    if (j > 0) { r3_t *U = &(S[j-1]); proc(V, U, &T); }
                    proc(W, V, &T);
                  }
                /* Save point: */
                S[j] = T;
              }
          }
      }
  }
