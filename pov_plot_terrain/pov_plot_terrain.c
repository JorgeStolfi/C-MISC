#define PROG_NAME "pov_plot_terrain"
#define PROG_DESC "generate a POV_Ray triangular mesh model from terrain data"
#define PROG_VERS "1.0"

#define pov_plot_terrain_C_COPYRIGHT "Copyright © 2010 by the State University of Campinas (UNICAMP)"

/* Last edited on 2010-07-02 13:35:09 by stolfi */

#define PROG_HELP \
  PROG_NAME  " \\\n" \
  "  [ -clip { in | out } {XMIN} {XMAX} {YMIN} {YMAX} {ZMIN} {ZMAX} ] \\\n" \
  "  [ -verbose ] \\\n" \
  "  [ < ] {TERRAIN_FILE} > {INC_FILE}"

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  PROG_INFO_DESC "\n" \
  "\n" \
  "OPTIONS\n" \
  PROG_INFO_OPTS "\n" \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  pgm(5), ppm(5), pnm(5), pgmselect(1).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created in jun/2010 by Jorge Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "  Info and experiments by Rafael Saracchini, IC-UNICAMP/MVL-UWE." \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  2010-06-22 J. Stolfi: added \"-crop\" option.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " pov_plot_terrain_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS
  
#define PROG_INFO_DESC \
  "  The program reads a list of points {(X[k],Y[k],Z[k])} where {Z[k]} is" \
  " assumed to be a function of {(X[k],Y[k])}.  It outputs a mesh of" \
  " triangles that interpolates those points.  The projection of the" \
  " mesh on the {(X,Y)} plane is the Delaunay traingulation of the" \
  " points {(X[k],Y[k])}.\n" \
  "\n" \
  "  The mesh is optionally clipped against a user-given box."
  
#define PROG_INFO_OPTS \
  "  -clip { in | out } {XMIN} {XMAX} {YMIN} {YMAX} {ZMIN} {ZMAX}\n" \
  "    This argument specifies a box in three-space, and requests that the" \
  " output contains only those parts of the mesh which" \
  " lie respectively inside and outside that box.\n" \
  "\n" \
  "  -verbose\n" \
  "    Produces diagnostic output."

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <argparser.h>
#include <unistd.h>
#include <string.h>
#include <i2.h>
#include <r2.h>
#include <r3.h>
#include <tabela.h>
#include <jsfile.h>
#include <quad.h>
#include <delaunay.h>
#include <delaunay_plot_POV.h>

/* PROTOTYPES */

typedef struct options_t
  { double scale;           /* Multiplier for Z values. */
    double offset;          /* Additive offset for Z values. */
    bool_t clip_in;         /* True to output only parts inside the clip box. */     
    bool_t clip_out;        /* True to output only parts outside the clip box. */     
    interval_t clip_int[3]; /* Clipping interval on each axis. */
  } options_t;
  /* Command line arguments. */

int main(int argc, char* argv[]);

options_t *parse_options(int argc, char **argv);
  /* Parses the command line arguments. */

/* REAL CODE */

int main(int argc, char** argv)
  {
    /* Parse command line arguments. */
    options_t *o = parse_options(argc, argv);
    
    double scale = 1.0;
    double offset = 0.0;
    r3_vec_t pts = read_data(stdin);
    output_triangles(stdout, pts, o);
    return 0;
	
  }

void output_triangles(FILE *wr, r3_vec_t pts, options_t *o)
  {
    int nsites = pts.ne;
    delaunay_site_t *st = notnull(malloc(nsites*sizeof(delaunay_site_t)), "no mem");
    double *z = notnull(malloc(nsites*sizeof(double)), "no mem");
    int i;
    for (i = 0; i < nsites; i++)
      { st[i].index = i;
        r3_t *p = &(pts.e[i]);
        st[i].p.c[0] = p->c[0];
        st[i].p.c[1] = p->c[1];
        z[i] = o->scale * p->c[2] + o->offset;
      }
    quad_arc_t e = delaunay_build (st, nsites);
    plot_delaunay_pov(wr, e, nsites, st, z, o->clip_in, o->clip_out, o->clip_int);
    free(height);
    fflush(wr);
  }


void plot_delaunay_POV (quad_arc_t e, delaunay_site_t *st, int nsites, double height[] ,char *filename,double scale,double offset)
  { 
      
    /* Create Postscript document or EPS figure stream. */
    FILE* wr = open_write(filename,TRUE);
    float minX,minY,maxX,maxY;
    minX = minY = -1;
    maxX = maxY = +1;
    float maxZ = offset;
    float minZ = offset;
    int i;
    for(i = 0; i < nsites; i++){
	maxZ = fmax(maxZ,height[i]);
	minZ = fmin(minZ,height[i]);
    }
    fprintf(wr,"#declare minX = %f;\n",minX);
    fprintf(wr,"#declare maxX = %f;\n",maxX);
    fprintf(wr,"#declare minY = %f;\n",minY);
    fprintf(wr,"#declare maxY = %f;\n",maxY);
    fprintf(wr,"#declare maxZ = %f;\n",maxZ);
    fprintf(wr,"#declare minZ = %f;\n",minZ);
    fprintf(wr,"#declare zeroLevel = %f;\n",offset);
    fprintf(wr,"#declare graph = \n");
    fprintf(wr,"  mesh{\n");
    delaunay_plot_POV_triangles(wr, e,height,"graph_tx_out");
    fprintf(wr,"  }\n");
    fprintf(wr,"#declare skirt = \n");
    fprintf(wr,"  mesh{\n");
    delaunay_plot_POV_skirt(wr,e,height,"skirt_tx_out");
    fprintf(wr,"  }\n");
    fclose(wr);
  }


