#define PROG_NAME "povwfc"
#define PROG_DESC "generate a POV-Ray description of a geomodel and wavefront"
#define PROG_VERS "1.1"

/* Copyright © 2005 by the State University of Campinas (UNICAMP). */
/* See the copyright, authorship, and warranty notice at end of file. */
/* Last edited on 2007-12-27 02:06:49 by stolfi */

#define PROG_HELP \
  "  " PROG_NAME  " \\\n" \
  "    -triName <name> \\\n" \
  "    [ -geoName <name> ] \\\n" \
  "    -outName <name> \\\n" \
  "    [ -lineWidth NUM ]\n"

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  ???\n" \
  "\n" \
  "AUTHORS\n" \
  "  Created in aug/2005 by Lucas Freitas and Jorge Stolfi at IC-UNICAMP."

#include <basic.h>
#include <mesh.h>
#include <geomodel.h>
#include <pov_utils.h>
#include <mesh_pov.h>
#include <geomodel_pov.h>

#include <r3.h>

#include <affirm.h>
#include <argparser.h>
#include <bool.h>
#include <sign.h>
#include <jsfile.h>
#include <jsstring.h>

#include <values.h>
#include <stdio.h>
#include <math.h>

typedef struct options_t /* Parsed command line options */
  { char *triName;        /* Mesh file (minus ".tri" extension) */
    char *geoName;        /* Geophysical model file (minus ".geo" extension) or NULL */
    char *outName;        /* Output POV-Ray file (minus ".inc" extension). */
    double lineWidth;     /* Diameter of edges and vertices (0 = don't plot them). */
  } options_t;

void pov_everything
  ( FILE *fpov,              /* POV-Ray file. */
    mesh_t *tri,             /* mesh_t to draw, or NULL */
    geomodel_t *geo,         /* geophysical model to draw, or NULL */
    int N,                   /* Mesh subdivision parameter. */
    double lineWidth,        /* Diameter of edges and vertices (in model units). */
    bool_t verbose           /* TRUE to print messages along the way. */
  );
  /* Writes to {fpov} POV-ray description of the mesh {tri} and model {geo}.
    Edges and vertices are represented by bars and balls with 
    diameter {lineWidth}. */

options_t *parse_options(int argn, char **argc);
mesh_t *read_named_mesh(char *triName);
geomodel_t *read_named_geomodel(char *geoName);

int main(int argc, char **argv)
  { options_t *o = parse_options(argc, argv);

    mesh_t *tri = read_named_mesh(o->triName);
    geomodel_t *geo = read_named_geomodel(o->geoName);
      
    /* Open the figure stream: */
    char *fname = NULL;
    asprintf(&fname, "%s.inc", o->outName);
    FILE *fpov = open_write(fname, TRUE);

    /* Plot the triangulation: */
    pov_everything(fpov, tri, geo, 1, o->lineWidth, TRUE);
    
    fclose(fpov);
    free(fname);
    return 0;
  }
  
mesh_t *read_named_mesh(char *triName)
  { FILE *rd = open_read(txtcat(triName, ".tri"), TRUE);
    mesh_t *tri = read_mesh(rd);
    fclose(rd);
    return tri;
  }
    
geomodel_t *read_named_geomodel(char *geoName)
  { 
    if (geoName == NULL)
      { return NULL; }
    else
      { 
        char *fname = NULL;
        asprintf(&fname, "%s.geo", geoName);
        FILE *rd = open_read(fname, TRUE);
        geomodel_t *geo = notnull(malloc(sizeof(geomodel_t)), "no mem");
        (*geo) = read_geomodel(rd);
        free(fname);
        fclose(rd);
        return geo;
      }
  }
    
void pov_everything
  ( FILE *fpov,              /* POV-Ray file. */
    mesh_t *tri,             /* mesh_t to draw, or NULL */
    geomodel_t *geo,         /* geophysical model to draw, or NULL */
    int N,                   /* Mesh subdivision parameter. */
    double lineWidth,        /* Diameter of edges and vertices (in model units). */
    bool_t verbose           /* TRUE to print messages along the way. */
  )
  { 
    mumble("geophysical model...\n");;
    pov_geomodel(fpov, geo, lineWidth/2); 

    mumble("mesh...\n");;
    pov_mesh(fpov, tri, N, lineWidth/2); 
  }

options_t *parse_options(int argc, char **argv)
  { argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    options_t *o = notnull(malloc(sizeof(options_t)), "no mem");

    argparser_get_keyword(pp, "-triName");                           
    o->triName = argparser_get_next(pp);  

    if (argparser_keyword_present(pp, "-geoName"))                           
      { o->geoName = argparser_get_next(pp); }
    else
      { o->geoName = NULL; }

    argparser_get_keyword(pp, "-outName");                               
    o->outName = argparser_get_next(pp);  

    if (argparser_keyword_present(pp, "-lineWidth"))
      { o->lineWidth = argparser_get_next_double(pp, 0.0, INF); }
    else
      { o->lineWidth = 0.0; }

    argparser_finish(pp);
    return o;
  }

