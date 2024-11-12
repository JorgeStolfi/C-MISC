/* Last edited on 2023-02-12 07:52:15 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <argparser.h>
#include <argparser_geo.h>
#include <bool.h>
#include <jsfile.h>
#include <affirm.h>
#include <assert.h>
#include <stdlib.h>
#include <filefmt.h>
#include <nget.h>
#include <fget.h>
#include <jsmath.h>
#include <.h>
#include <.h>
#include <r3.h>

typedef struct config_t
  { int n;          /* Number of particles. */
    double *m;      /* Masses of particles. */
    r3_t *p;        /* Positions of particles. */ 
    r3_t *v;        /* Velocities of particles. */ 
    /* Domain: */
    r3_t sz;        /* Dimensions of cell. */
    bool_t per[3];  /* Periodicity of cell. */
    /* Global parameters: */
    double g;       /* Gravitational attraction. */
    double T;       /* Temperature of walls. */
  } config_t;
  /* A configuration of the system, consisting of {n} particles with
    masses {m[0..n-1]}, whose centers have positions {p[0..n-1]} and
    velocities {v[0..n-1]}, confined to a box with size {sz}.

    If {per[i]} is true, motions are unconstrained and the configuration
    is periodic along axis {i}; that is, the contents of the box
    are repeated indefinitely in that direction, and a particle whose
    center  leaves the box through a wall perpendicular to axis {i}
    imemdiately re-enters through the opposite wall with the same 
    velocity. 

    If {per[i]} is false, then particles bounce on the walls
    perpendicular to axis {i} as if they are hard barriers. */

typedef struct options_t 
  { 
    char *ini_name;   /* Name of initial configuration file. */
    double t0;        /* Nominal time at start of simulation. */
    double t1;        /* Nominal time at end of simulation. */
    double save_step; /* Interval for snapshot saving. */
    chat *out_pref;   /* Prefix for output file names. */
  } options_t;

int main (int argc, char **argv)
  {

    options_t *o = lqs_get_options(argc, argv);

    /* Get the initial configuration: */
    config_t *C = lqs_read_config(o->ini_name);

    /* Integrate over specified time interval: */
    lqs_integrate(C, o);

    /* Write the final configuration: */
    lqs_write_config(C, o->out_pref, "fin");

    /* Cleanup */
    return 0;
  }

#define config_TYPE "liqsim_config_t"
#define config_VERSION "2009-06-26"

config_t *lqs_read_config(char *f_name)
  {
    int i;
    FILE *rd = open_read(f_name, TRUE);
    config_t *C = notnull(malloc(sizeof(config_t)), "no mem");

    filefmt_read_header(rd, config_TYPE, config_VERSION);

    /* Domain: */
    C->sz = nget_r3_t(rd, "size"); fget_eol(rd);
    nget_name_eq(rd, "per"); 
    for (i = 0; i < 3; i++) { C->per[i] = fget_bool(rd); }
    fget_eol(rd);

    /* Global parameters: */
    C->g = nget_double(rd, "g");  fget_eol(rd);
    C->T = nget_double(rd, "T");  fget_eol(rd);

    /* Number of particles: */
    int n = nget_int32(rd, "n"); fget_eol(rd);
    C->n = n;

    /* Particle masses, positions, and velocities: */
    C->m = notnull(malloc(n*sizeof(double)), "no mem");
    C->p = notnull(malloc(n*sizeof(r3_t)), "no mem");
    C->v = notnull(malloc(n*sizeof(r3_t)), "no mem");
    for (i = 0; i < n; i++)
      { int ir = fget_int32(rd);
        demand(i == ir, "sequence error");
        C->m[i] = fget_double(rd, 0.0, INF);
        fget_rn(rd, 3, &(C->p[i].c));
        fget_rn(rd, 3, &(C->v[i].c));
        fget_eol(rd);
      }

    filefmt_read_footer(rd, config_TYPE);
    if (rd != stdin) { fclose(rd); }
    return C;
  }

void lqs_write_config(config_t *C, char *pref, char *tag)
  {
    int i;
    FILE *wr = open_write(f_name, TRUE);

    filefmt_write_header(wr, config_TYPE, config_VERSION);

    /* Domain: */
    fprintf(wr, "size = %22.16e %22.16e %22.16e\n", C->sz.c[0], C->sz.c[1], C->sz.c[2]);
    fprintf(wr, "per = %c %c %c\n", "FT"[C->per[0]], "FT"[C->per[1]], "FT"[C->per[2]]);
    
    /* Global parameters: */
    fprintf(wr, "g = %22.16e\n", C->g);
    fprintf(wr, "T = %22.16e\n", C->T);

    /* Number of particles: */
    int n = C->n;
    fprintf(wr, "n = %d\n", n);

    /* Particle masses, positions, and velocities: */
    int digs = digits(n-1);
    for (i = 0; i < n; i++)
      { fprintf(wr, "%*d\n", digs, i);
        fprintf(wr, " %22.16e", C->m[i]);
        r3_gen_print(wr, " %22.16e", &(C->p[i]), " ", "", "");
        r3_gen_print(wr, " %22.16e", &(C->v[i]), " ", "", "");
        fprintf(wr, "\n");
      }

    filefmt_write_footer(wr, config_TYPE);
    if (wr != stout) { fclose(wr); }
  }

void lqs_integrate(config_t *C, options_t *o)
  {
    

  }

options_t *lqs_get_options(int argc, char **argv)
  {
    options_t *o = notnull(malloc(sizeof(options_t)), "no mem");
    
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help_info(pp, PROG_HELP, PROG_INFO);
    argparser_parse_help_info_options(pp);
    
    
    argparser_skip_parsed(pp);
    
    argparser_finish(pp);
    return o;
  }
