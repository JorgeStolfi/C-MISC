#define PROG_NAME "quadwfc"
#define PROG_DESC "seismic wavefront simulator based on quad-edge"
#define PROG_VERS "1.1"

/* Copyright © 2005 by the State University of Campinas (UNICAMP). */
/* See the copyright, authorship, and warranty notice at end of file. */
/* Last edited on 2023-02-12 11:19:43 by stolfi */

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    -outName NAME \\\n" \
  "    -maxIter NUM \\\n" \
  "    -timeStep NSECS \\\n" \
  "    -signature STRING \\\n" \
  "    -edgeLength NUM \\\n" \
  "    < INFILE.dat\n"
 
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
  "  Created in aug/2005 by Lucas Freitas and Jorge Stolfi at IC-UNICAMP.\n" \
  "  Based on \"simul.c\" written in jul/2005 by Lucas Freitas, IMECC-UNICAMP."
  
#include <basic.h>
#include <geomodel.h>
#include <wavefront.h>
#include <init_wf.h>
#include <output.h>
#include <wave_prop.h>

#include <argparser.h>

#include <string.h>

typedef struct options_t 
  { char *outName;      /* Prefix for file names. */
    int32_t maxIter;        /* Number of iterations to run. */
    double timeStep;    /* Time step in seconds. */
    double edgeLength;  /* Ideal edge length. */
    char *signature;    /* Signature of rays to be followed. */
  } options_t;
  
options_t *get_options(int32_t argc, char **argv);
double parse_double(int32_t *argn, int32_t argc, char **argv, double lo, double hi);
void arg_error(char *msg, int32_t argn, char *arg);
geomodel_t make_geomodel(void);

int32_t main (int32_t argc, char **argv)
{
  /* Parse command line arguments: */
  options_t *o = get_options(argc, argv);

  /* Build geophysical model: */
  geomodel_t geo = make_geomodel();
  output_model(&geo, o->outName);

  /* Build initial wavefront: */
  r3_t orig = (r3_t){{ 5000, 5000, 4000 }};
  double rad = 4.0*o->edgeLength; /* Radius of initial wavefront. */
  int32_t imd = 0; /* Initial layer index. */
  char *sgn = o->signature;
  wavefront_t wf = init_wf(&orig, rad, geo.md[0].v_P, o->edgeLength, imd, sgn);
  // wavefront_t wf = init_wf(&orig, 2500, geo.md[0].v_P, 100);
  
  output_wave(&wf, o->outName, -1);
  for (int32_t i = 1; i <= o->maxIter; i++)
    {
      wave_prop(&wf,o->timeStep, &(geo), o->edgeLength);
      plot_wave(&wf, o->outName, i);
      if ((i < 5) || (i % 5 == 0)) { output_wave(&wf, o->outName, i); }
    }
  output_wave(&wf, o->outName, INT32_MAX);
    
  return 0;
}

geomodel_t make_geomodel(void)
  {
    geomodel_t geo;
    
    int32_t nrf = 1; /* Number of reflectors */
    
    double xinf = 0, xsup = 10000;
    double yinf = 0, ysup = 10000;
    double zinf = 0, zsup = 5100;
    
    /* Domain bounding box: */ 
    geo.bb[0] = (interval_t){{ xinf, xsup }};
    geo.bb[1] = (interval_t){{ yinf, ysup }};
    geo.bb[2] = (interval_t){{ zinf, zsup }};

    /* Layers: */
    geo.nmd = nrf+1;
    geo.md = notnull(malloc((geo.nmd+1)*sizeof(medium_t)), "no mem");

    geo.md[0] = make_medium(0, 50, 40);
    geo.md[1] = make_medium(1, 60, 30);
    
    /* Reflectors (inter-layer boundaries): */
    geo.nrf = nrf;
    geo.rf = notnull(malloc((geo.nrf)*sizeof(reflector_t)), "no mem");
    
    geo.rf[0] = make_sinusoidal_reflector
      ( &(geo.md[0]), &(geo.md[1]), 
        17, 17, 
        xinf, xsup,
        yinf, ysup,
        2000, 2100,
        -1.50, 0.75
      );
    return geo;
  }

options_t *get_options(int32_t argc, char **argv)
  {
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    options_t *o = notnull(malloc(sizeof(options_t)), "no mem"); 
    
    argparser_get_keyword(pp, "-outName");                           
    o->outName = argparser_get_next(pp);  

    argparser_get_keyword(pp, "-maxIter");                           
    o->maxIter = (int32_t)argparser_get_next_int(pp, 1, 100000);  

    argparser_get_keyword(pp, "-timeStep");                           
    o->timeStep = argparser_get_next_double(pp, 1.0e-10, 1.0e10);  

    argparser_get_keyword(pp, "-edgeLength");                           
    o->edgeLength = argparser_get_next_double(pp, 1.0e-10, 1.0e10);  

    argparser_get_keyword(pp, "-signature");                           
    o->signature = argparser_get_next(pp);  

    argparser_finish(pp);

    return o;
  }
