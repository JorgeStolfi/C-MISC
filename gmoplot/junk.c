// Last edited on 2005-08-21 12:13:46 by stolfi

//  void compute_face_matrices(quad_arc_t e);
//    /* (Re)computes the barycentric-to-Cartesian matrix {f->b2c} and its
//      inverse {f->c2b} for the face {f = LEFT(e)}, assuming {e} is its
//      reference side. Assumes that all {Org} fields are set. */
//     
// void compute_face_geometry(plotlist_t *tri);
//   /* (Re)computes all derived geometry information, from the 
//     site coordinates. 
//     
//     Specifically, calls {compute_face_matrices(e)} for the reference edges
//     {e = tri->side[i]} of all faces in {tri}. Assumes that the {Org}
//     and {Left} fields are properly set, and the table {tri->side} is
//     up-to-date.
//     
//     This procedure is meant to be used after changing the site positions
//     in a mesh, provided that the topology hasn't changed. */

//  void compute_face_geometry(mesh_t *tri) 
//    { quad_arc_vec_t side = tri->side;
//      int i;
//  
//      for (i = 0; i < side.nel; i++) 
//        { quad_arc_t ei = side.el[i];
//          compute_face_matrices(ei);
//        }
//    }
//  
// void compute_face_matrices(quad_arc_t e)
//   { face_t *f = LEFT(e);
//     r3x3_t *b2c = &(f->b2c);
//     r3x3_t *c2b = &(f->c2b);
//     
//     sample_t u = DEST(ONEXT(e))->pos;
//     sample_t v = ORG(e)->pos;
//     sample_t w = DEST(e)->pos;
// 
//     *b2c = (r3x3_t){{
//       {u.c[0], v.c[0], w.c[0]},
//       {u.c[1], v.c[1], w.c[1]},
//       {u.c[2], v.c[2], w.c[2]}
//     }};
//     r3x3_inv(b2c, c2b);
//   }
