#define PROG_NAME "wikipedia_vorticity"
#define PROG_DESC "generate data file for graphic visualization of vorticity"
#define PROG_VERS "1.0"

#define wikipedia_vorticity_C_COPYRIGHT "Copyright © 2010 by the State University of Campinas (UNICAMP)"

/* Last edited on 2023-02-20 04:52:18 by stolfi */

#define PROG_HELP \
  PROG_NAME  " \\\n" \
  "  \\\n" \
  "  -dim {DIM} \\\n" \
  "  -field {FNAME} -clock {CLOCK_INI} {CLOCK_STEP} {CLOCK_FIN} \\\n" \
  "  -outPrefix {OUT_PREFIX} \\\n" \
  "  [ -verbose ]"

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
  "  thenicemoonoutside(1).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created sep/2012 by Jorge Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  2020-09-28 J. Stolfi: created.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " wikipedia_vorticity_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS
  
#define PROG_INFO_DESC \
  "  The program outputs a set of POV-Ray \".inc\" files containing declarations for" \
  " a frame in a movie of particle paths in a velocity field.\n" \
  "\n" \
  "  The mesh is optionally clipped against a user-given box."
  
#define PROG_INFO_OPTS \
  "  -dim {DIM}\n" \
  "    This argument specifies the dimension of the flow.\n" \
  "\n" \
  "  -field {FNAME}\n" \
  "    This argument specifies the velocity field.\n" \
  "\n" \
  "  -clock {CLOCK_INI} {CLOCK_STEP} {CLOCK_FIN}\n" \
  "    This argument specifies the relative times of the frames in the" \
  " animation, typically between 0 (first frame) and 1 (last frame).  The" \
  " frame times will be {CLOCK_INI + k*CLOCK_STEP} for {k=0,1,...}, while" \
  " not exceeding {CLOCK_STEP}. \n" \
  "\n" \
  "  -outPrefix {OUT_PREFIX}\n" \
  "    Specifies the prefix for output file names.  The files" \
  " will be called \"{PREFIX}-{N.NNNNN}.inc\" where {N.NNNNN} is the" \
  " clock value, with five digits after the period.\n" \
  "\n" \
  "  -verbose\n" \
  "    Produces diagnostic output."

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <argparser.h>
#include <string.h>

#include <r3.h>
#include <vec.h>

/* PROTOTYPES */

typedef struct options_t
  { int dim;                /* Either 2 or 3. */
    char * field;           /* Velocity field name. */
    double clock_ini;       /* Clock of first frame. */
    double clock_step;      /* Clock increment between frame. */
    double clock_fin;       /* Clock of last frame must not exceed this. */
    char * outPrefix;       /* Output file prefix. */
    bool_t verbose;         /* True for diagnostics. */     
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
    
    r3_vec_t p0; = generate_particles(o->dim, density);
    for (i = 0; i <= +nc; i++)
      { generate_streakline(
    return 0;
	
  }

void euler_integral
  ( velocity_field_t *vel,
    r3_t p0,
    double ds,
    double dt,
    path_predicate_t *stop,
    double s0,
    double hs,
    double t0,
    double ht,
    path_predicate_t *plop );
  /* Integrates the velocity field {vel(p,t)} to obtain a path {P(t)}.
    
    Integration starts at {P(0) = p0} and extends in both directions
    until {stop(P(t),s(t),t)} is true.  Each step is at most {ds} 
    long and at takes at most {dt} time.
    
    In these formulas {t} is time and {s(t)} is the signed
    length along the path from {P(0)} to {P(t)}.
    
    The velocity field must be smooth and nonzero at all points and times 
    along the path.  Result is unpredictable otherwise.
    
    If {hs} is nonzero, calls {plop(P(t),s(t),t)} whenever {s(t)} is
    {s0 + k*hs} for some integer {k}; in that case {hs} must be greater
    than {ds}.  Likewise, if {ht} is nonzero, calls {plop(P(t),s(t),t)} 
    whenever {t} is {t0 + k*ht} for some integer {k}; in that case {ht} 
    must be greater than {dt}. */


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
    delaunay_plot_diagram_pov(wr, e, nsites, st, z, o->clip_in, o->clip_out, o->clip_int);
    free(height);
    fflush(wr);
  }


void delaunay_plot_diagram_POV (quad_arc_t e, delaunay_site_t *st, int nsites, double height[] ,char *filename,double scale,double offset)
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


