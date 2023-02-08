/* See {irt_scene.h}. */
/* Last edited on 2008-07-14 20:33:21 by stolfi */

#define _GNU_SOURCE
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <malloc.h>

#include <frgb.h>
#include <frgb_ops.h>
#include <r3.h>
#include <ia.h>
#include <aa.h>
#include <flt.h>
#include <iaeval.h>
#include <flteval.h>
#include <affirm.h>
#include <pcode.h>
#include <jsfile.h>
#include <jsstring.h>

#include <irt_scene.h>
#include <irt_inter.h>

/*** INTERNAL PROTOTYES ***/

void irt_set_parm_defaults(scene_t *sc);
  /* 
    Assigns default values for all $*sc$ fields that can be 
    specified in the .parms file. */
    
void irt_read_parm_file (char *dir_name, char *scene_name, scene_t *sc);
  /* Reads the viewing, shading, and imaging parameters from
    the <scene_name>.parms file. */

void irt_read_pcode_file (char *dir_name, char *proc_name, pcode_proc_t *proc);
  /* 
    Reads the pseudocode of the scene's characteristic function
    from the <proc_name>.pcode file. */

void irt_compute_screen_axes(view_t *vw);
  /* 
    Computes vw->screen_x, vw->screen_y
    from vw->focus, vw->observer, vw->up, and vw->pixel_size */
    
void irt_compute_looks_flags(looks_t *lk);
  /* Computes the flags $lk->has_shiny$, $lk->has_mirror$, $lk->has_transp$ */

void irt_allocate_regs_and_stack(shape_t *sh);
  /* 
    Allocates the stack and registers needed for the evaluation of {sh->proc}
    with any of the supported arithmetic models. */
    
void irt_read_light_parms(FILE *parms_file, light_t *lt, int *nlines);
  /* Reads the parameters of a light source, starting at the next
    line of parms_file, up to and including the "end_light" directive.
    If the "ref_point" parameter is not specified, sets it to infinity. */

void irt_compute_light_parameters(light_t *lt, r3_t focus);
  /* 
    Computes $lt->ref_distance$ for light $lt$.
    The $focus$ is used as the default ref_point. */
    
/* 
  The following procedures continue parsing the current line by 
  calling $strtok(NULL, " ")$.  They bomb on errors. */

r3_t irt_read_point(void);
r3_t irt_read_vector(void);
  /* 
    Parses three real numbers (Cartesian coordinates). */
  
frgb_t irt_read_rgb (void);
  /* 
    Parses three real numbers (intensities, usually in [0__1]). */
  
char *irt_read_name(void);
  /* Parses an identifier. */
    
int irt_read_int(void);
  /* Parses an integer */

double irt_read_number(void);
  /* Parses a real number. */
    
pixel_num_t irt_read_pixel_num(void);
  /* Parses a pair of pixel indices (row and col). */

int irt_atoi(char *s);
double irt_atod(char *s);

/*** IMPLEMENTATIONS ***/

void irt_read_scene (char *dir_name, char *scene_name, scene_t *sc)
  { int i;
    sc->name = scene_name;
    irt_read_parm_file(dir_name, scene_name, sc);

    irt_compute_looks_flags(&(sc->looks));
    irt_compute_screen_axes(&(sc->view));
    for (i=0; i<sc->num_lights; i++)
      irt_compute_light_parameters(&(sc->light[i]), sc->view.focus);

    irt_read_pcode_file(dir_name, sc->shape.proc_name, &(sc->shape.proc));
    irt_allocate_regs_and_stack (&(sc->shape));
  }
  
void irt_set_parm_defaults(scene_t *sc)
  { 
    sc->print_ray = (pixel_num_t) {-1, -1};
    sc->plot_ray =  (pixel_num_t) {-1, -1};

    /* Shape: */
    sc->shape.proc_name = NULL;
    
    /* Looks: */
    sc->looks.ambient_color    = (frgb_t) {{0.05, 0.05, 0.05}};
    sc->looks.matte_color      = (frgb_t) {{0.05, 0.60, 0.05}};
    sc->looks.shiny_color      = (frgb_t) {{0.20, 0.20, 0.10}};
    sc->looks.shiny_spread     = 0.9;
    sc->looks.max_bounces      = 3;
    sc->looks.background_color = (frgb_t) {{0.10, 0.10, 0.60}};
    
    /* Viewing: */
    sc->view.image_width  = 128;
    sc->view.image_height = 128;
    sc->view.pixel_size = 0.3 / 500.0; 
    sc->view.observer = (r3_t) {{100.0, 10.0, 30.0}};
    sc->view.focus    = (r3_t) {{0.0, 0.0, 0.0}};
    sc->view.up       = (r3_t) {{0.0, 0.0, 1.0}};
    
    /* Lights: */
    sc->num_lights = 1;
    { /* Default light number 0 */
      sc->light[0].position = (r3_t) {{100.0, 50.0, 200.0}};
      sc->light[0].color   = (frgb_t) {{0.8, 0.8, 0.8}};
      sc->light[0].ref_point = (r3_t) {{Infinity, Infinity, Infinity}};
    }
  }

void irt_read_parm_file (char *dir_name, char *scene_name, scene_t *sc)
  { char *fname = NULL;
    asprintf(&fname, "%s/%s.parms", dir_name, scene_name);
    FILE *parm_file = open_read(fname, TRUE);

    char *s; /* Current line */
    int nlines = 0;
    int default_num_lights;
    
    irt_set_parm_defaults(sc);
    
    /* Discard default lights if user gives any lights: */
    default_num_lights = sc->num_lights;
    sc->num_lights = 0;
    
    fprintf(stderr, "----------------------------------------------\n");
    nlines=0;
    s = read_line(parm_file);
    while (s != NULL)
      {
        char *t;
        nlines++;
        fprintf(stderr, "%s\n",s);
        t = strtok(s, " ");
        if (t == NULL) 
          { /* Blank line, ignore: */ }
        else if (strcmp(t, "proc_name") == 0)
          { sc->shape.proc_name = irt_read_name(); }
        else if (strcmp(t, "arithmetic") == 0)
          { char *a = strtok(NULL, " ");
            if (strcmp(a, "IA") == 0) 
              { sc->arithmetic = arith_ia; }
            else if (strcmp(a, "IA/DIFF") == 0) 
              { sc->arithmetic = arith_ia_diff; }
            else if (strcmp(a, "AA") == 0)
              { sc->arithmetic = arith_aa; }
            else if (strcmp(a, "MIX") == 0)
              { sc->arithmetic = arith_mix; }
            else
              { fatalerror ("irt_read_parm_file: unknown arithmetic"); }
          }
        else if (strcmp(t, "print_ray") == 0)
          { sc->print_ray = irt_read_pixel_num(); }
        else if (strcmp(t, "plot_ray") == 0)
          { sc->plot_ray = irt_read_pixel_num(); }
        else if (strcmp(t, "ambient_color") == 0)
          { sc->looks.ambient_color = irt_read_rgb(); }
        else if (strcmp(t, "matte_color") == 0)
          { sc->looks.matte_color = irt_read_rgb(); }
        else if (strcmp(t, "shiny_color") == 0)
          { sc->looks.shiny_color = irt_read_rgb(); }
        else if (strcmp(t, "shiny_spread") == 0)
          { sc->looks.shiny_spread = irt_read_number(); }
        else if (strcmp(t, "max_bounces") == 0)
          { sc->looks.max_bounces = irt_read_int(); }
        else if (strcmp(t, "background_color") == 0)
          { sc->looks.background_color = irt_read_rgb(); }
        else if (strcmp(t, "image_width") == 0)
          { sc->view.image_width = irt_read_int(); }
        else if (strcmp(t, "image_height") == 0)
          { sc->view.image_height = irt_read_int(); }
        else if (strcmp(t, "pixel_size") == 0)
          { sc->view.pixel_size = irt_read_number(); }
        else if (strcmp(t, "observer") == 0)
          { sc->view.observer = irt_read_point(); }
	else if (strcmp(t, "focus") == 0)
          { sc->view.focus = irt_read_point(); }
	else if (strcmp(t, "up") == 0)
          { sc->view.up = irt_read_vector(); }
        else if (strcmp(t, "light_source") == 0)
          { affirm(sc->num_lights < MAXLIGHTS, "irt_read_parm_fil: too many lights");
            irt_read_light_parms(parm_file, &(sc->light[sc->num_lights]), &nlines);
            (sc->num_lights)++;
	  }
        else
          { fatalerror ("irt_read_parm_file: unrecognized keyword"); }
        
        /* Check for bogus input: */
        affirm (strtok(NULL, " ") == NULL, "irt_read_parm_file: superfluous parameters");
        s = read_line(parm_file);
      }
    if (sc->num_lights == 0) sc->num_lights = default_num_lights;
    if (sc->shape.proc_name == NULL) 
      { fatalerror ("irt_read_parm_file: missing parameter proc_name"); }
    fclose(parm_file);
    fprintf(stderr, "----------------------------------------------\n");
    free(fname);
  }

void irt_read_light_parms(FILE *parm_file, light_t *lt, int *nlines)
  {
    char *s;
    int done = 0;
    
    s = read_line(parm_file);
    
    lt->position = (r3_t) {{Infinity, Infinity, Infinity}};
    lt->ref_point = (r3_t) {{Infinity, Infinity, Infinity}};
    lt->color = (frgb_t) {{0.8, 0.8, 0.8}};

    while (s != NULL)
      { char *t;
        (*nlines)++;
        fprintf(stderr, "%s\n",s);
        t = strtok(s, " ");
        if (strcmp(t, "position") == 0)
          { lt->position = irt_read_point(); }
        else if (strcmp(t, "ref_point") == 0)
          { lt->ref_point = irt_read_point(); }
        else if (strcmp(t, "color") == 0)
          { lt->color = irt_read_rgb(); }
        else if (strcmp(t, "end_light") == 0)
          { done = 1; }
        else
          { fatalerror ("irt_read_light_parms: unrecognized keyword"); }

        /* Check for bogus input: */
        affirm (strtok(NULL, " ") == NULL, "irt_read_parm_file: superfluous parameters");
        if (done) 
          { if (lt->position.c[0] == Infinity)
	      { fatalerror("irt_read_light_parms: missing light position"); }
            return;
          }
        s = read_line(parm_file);
      }
    fatalerror("irt_read_light_parms: missing end_light"); 
  }
    
#define DEGREETORADIAN (M_PI/180.)

void irt_compute_screen_axes(view_t *vw)
  {
    r3_t sx, sy, sz;
    double dist, pix_dist;

    r3_sub(&vw->observer, &vw->focus, &sz);
    dist = r3_dir(&sz, &sz);

    r3_cross(&(vw->up), &sz, &sx);
    (void) r3_dir(&sx, &sx);
    
    r3_cross(&sz, &sx, &sy);
    (void) r3_dir(&sy, &sy);

    pix_dist = vw->pixel_size * dist;
    
    fprintf(stderr, "pixel size at focus = %f\n", pix_dist);
    fprintf(stderr, "image width at focus = %f\n", 
      ((double) vw->image_width) * pix_dist
    );
    fprintf(stderr, "image height at focus = %f\n", 
      ((double) vw->image_height) * pix_dist
    );
    
    r3_scale(pix_dist, &sx, &(vw->screen_x));
    r3_scale(pix_dist, &sy, &(vw->screen_y));
  }

void irt_compute_looks_flags(looks_t *lk)
  { lk->has_shiny = frgb_is_all_zeros(&(lk->shiny_color));
    if (lk->shiny_spread == 1) 
      { lk->has_mirror = lk->has_shiny; 
        lk->has_shiny = 0;
      }
    lk->has_transp = 0;
  }

void irt_compute_light_parameters(light_t *lt, r3_t focus)
  { if (lt->ref_point.c[0] == Infinity) lt->ref_point = focus;
    lt->ref_distance = r3_dist(&(lt->ref_point), &(lt->position));
  } 

void irt_read_pcode_file(char *dir_name, char *proc_name, pcode_proc_t *proc)
  { char *fname = NULL;
    asprintf(&fname, "%s/%s.pcode", dir_name, proc_name);
    FILE *pcode_file = open_read(fname, TRUE);
    *proc = pcode_parse(pcode_file);
    fclose(pcode_file);
    free(fname);
    fprintf(stderr, "----------------------------------------------\n");
    pcode_print(stderr, *proc);
    fprintf(stderr, "----------------------------------------------\n");
    affirm(proc->nin == 4, "irt_read_pcode_file: wrong number of arguments");
    affirm(proc->nout == 1, "irt_read_pcode_file: wrong number of results");
  }
  
void irt_allocate_regs_and_stack (shape_t *sh)
  { 
    int nr = sh->proc.nregs;
    int ns = sh->proc.nstack;
    
    /* Work areas for evaluating $proc$ along a ray: */
    sh->ia_regs  = (Interval *) malloc(nr * sizeof(Interval));
    sh->ia_stack = (Interval *) malloc(ns * sizeof(Interval));
    
    sh->id_regs  = (IntervalDiff *) malloc(nr * sizeof(IntervalDiff));
    sh->id_stack = (IntervalDiff *) malloc(ns * sizeof(IntervalDiff));
    
    sh->aa_regs  = (AAP *) malloc(nr * sizeof(AAP));
    sh->aa_stack = (AAP *) malloc(ns * sizeof(AAP));
    
    /* Work areas for computing the surface normal: */
    sh->nrm_regs =   (FloatDiff4 *) malloc(nr * sizeof(FloatDiff4));
    sh->nrm_stack =  (FloatDiff4 *) malloc(ns * sizeof(FloatDiff4));

    /* Work areas for evaluating $proc$ at a point (for debugging) */
    sh->flt_regs =   (Float *) malloc(nr * sizeof(Float));
    sh->flt_stack =  (Float *) malloc(ns * sizeof(Float));
  }

r3_t irt_read_point(void)
  { r3_t p;
    p.c[0] = irt_atod(strtok(NULL, " "));
    p.c[1] = irt_atod(strtok(NULL, " "));
    p.c[2] = irt_atod(strtok(NULL, " "));
    return(p);
  }
  
r3_t irt_read_vector(void)
  { r3_t v;
    v.c[0] = irt_atod(strtok(NULL, " "));
    v.c[1] = irt_atod(strtok(NULL, " "));
    v.c[2] = irt_atod(strtok(NULL, " "));
    return(v);
  }
  
int irt_read_int(void)
  { return(irt_atoi(strtok(NULL," "))); }
  
double irt_read_number(void)
  { return(irt_atod(strtok(NULL," "))); }
  
char *irt_read_name(void)
  { return(txtcat(strtok(NULL," "), "")); }
  
frgb_t irt_read_rgb (void)
 { frgb_t c;
   c.c[0] = irt_atod(strtok(NULL, " "));
   c.c[1] = irt_atod(strtok(NULL, " "));
   c.c[2] = irt_atod(strtok(NULL, " "));
   return(c);
 }

pixel_num_t irt_read_pixel_num(void)
  { pixel_num_t p;
    p.col = irt_atoi(strtok(NULL, " "));
    p.row = irt_atoi(strtok(NULL, " "));
    return(p);
  }
  
int irt_atoi(char *s)
  { affirm(s != NULL, "irt_atoi: integer expected");
    return(atoi(s));
  }

double irt_atod(char *s)
  { affirm(s != NULL, "irt_atoi: number expected");
    return(atof(s));
  }

