/* Last edited on 2014-07-23 19:49:58 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>

#include <gem.h>
#include <gem_print_graph.h>
#include <gem_print_face.h>

#include <gch_vector.h>
#include <gch_vector_list.h>
#include <gch_convexhull.h>

int main(int argc, char **argv)
  {
    if (argc != 4)
      { fprintf(stderr,"usage: %s <pointFile>  <graphfile>  <facefile> \n",argv[0]);
        return -1;
      }

    char *fnamePt = argv[1]; /* Name of input points file. */
    char *fnameGr = argv[2]; /* Name of output graph file. */
    char *fnameFc = argv[3]; /* Name of output 2-face file. */
    
    gch_vector_list_t *vl = gch_vector_list_read(fnamePt);
    fprintf(stderr, "Finding convex hull...\n");
    int dim;
    gem_ref_t gem = gch_convexhull_find(vl, &dim);

    if (dim <= 3)
      { fprintf(stderr, "Printing graph...\n");
        gem_print_graph_named(fnameGr,gem,dim);
      }
    if (dim >= 3)
      { fprintf(stderr, "Printing faces...\n");
        gem_print_face_2d_vertices_named(fnameFc, gem, dim);
      }
    gem_component_free(gem);

    gch_vector_list_free(vl);
    return 0;
  }
