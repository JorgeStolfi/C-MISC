#define PROG_NAME "irt_main"
#define PROG_DESC "prototype of an interval ray tracer"
#define PROG_VERS "1.0"

/* Last edited on 2017-02-26 02:52:22 by stolfilocal */

#define irt_main_C_COPYRIGHT \
  "Copyright © 2005 by the State University of Campinas (UNICAMP)"
  /* See {argparser_help_info_STANDARD_RIGHTS} in {argparser.h}. */
  /* See {argparser_help_info_NO_WARRANTY} in {argparser.h}. */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>

#include <ia.h>
#include <aa.h>
#include <r3.h>
#include <hr3.h>
#include <r3.h>
#include <frgb.h>
#include <affirm.h>
#include <jstime.h>
#include <jsfile.h>
#include <jsstring.h>
#include <argparser.h>

#include <irt_info.h>

#include <irt_out.h>
#include <irt_scene.h>
#include <irt_evray.h>
#include <irt_shade.h>
#include <irt_inter.h>

#define DEBUG 0

typedef struct irt_options_t {
  char *scene_dir;
  char *scene_name;
} irt_options_t;

/*** INTERNAL PROTOTYES ***/

int main(int argc, char **argv);

irt_options_t *irt_get_options(int argc, char **argv);

r3_t irt_pixel_center(view_t *vw, int ip, int jp);

void irt_write_eval_count_image(char *scene_dir, char *scene_name, int width, int height, int count[]);
  /* Writes a ".pgm" image with the evaluation counts. */

PSStream *irt_open_postscript_stream(char *dir, char *name, int col, int row, char *tag);
  /* Opens a Postscript stream called "{dir}/out/{name}-{col}-{row}{tag}.ps"
    for the pixel {col,row}.  The picture layout has one picture per page
    with equal H and V dimensions. */

void irt_close_postscript_stream(PSStream *ps);
  /* Closes the Postscript stream {ps}, properly finishing the current plot. */

/*** MAIN PROGRAM ***/

int main(int argc, char **argv)
  {
    ia_init();
    aa_init();

    irt_options_t *o = irt_get_options(argc, argv);
    char *scene_dir = o->scene_dir;
    char *scene_name = o->scene_name;
    
    /* Read the scene description {sc}: */
    scene_t sc;
    irt_read_scene(scene_dir, scene_name, &sc);
    int width = sc.view.image_width;
    int height = sc.view.image_height;
    
    /* Get homogeneous coordinates {eye} of observer: */
    hr3_point_t eye = hr3_from_r3(&sc.view.observer);
    hr3_L_inf_normalize_point(&eye); /* To improve the numerics. */

    /* Open the output image file: */
    FILE *out_img_file = irt_open_color_image_file(scene_dir, scene_name, "", width, height);
    frgb_t out_img_row[width];

    /* Allocate the array {eval_count} of function evaluation counts per pixel: */
    int *eval_count = (int *)malloc(sizeof(int)*width*height);
    int irt_max_evals = 0; /* Mex eval count per ray. */

    double start_time = user_cpu_time_usec();
    int prev_evals;
    int ip, jp;
    for (ip = 0; ip < height; ip++)
      { for(jp = 0; jp < width; jp++)
          { r3_t pix;        /* Pixel center */
            r3_t dir;        /* Ray's direction */

            frgb_t color;
            int print_ray, plot_ray;

            print_ray = (ip == sc.print_ray.row) & (jp == sc.print_ray.col);
            plot_ray = (ip == sc.plot_ray.row) & (jp == sc.plot_ray.col);

            pix = irt_pixel_center(&(sc.view), ip, jp);
            r3_sub(&pix, &(sc.view.observer), &dir);
            (void) r3_dir(&dir, &dir);

            if (print_ray)
              { fprintf(stderr, ". . . . . . . . . . . . . . . . . . . . . . . . . . . . . .\n");
                fprintf(stderr, " col = %d row = %d\n", jp, ip);
              }
              
            PSStream *ps = 
              ( plot_ray ? 
                irt_open_postscript_stream(scene_dir,scene_name,jp,ip,"-debug") : 
                NULL
              );

            prev_evals = irt_num_evals;

            color = irt_compute_scene_color(&sc, &eye, &dir, sc.looks.max_bounces, print_ray, ps);
            
            if (plot_ray) { irt_close_postscript_stream(ps); }

            int nevals = irt_num_evals - prev_evals; /* Count of evals for this ray. */
            eval_count[(ip*width) + jp] = nevals;
            if (nevals > irt_max_evals) { irt_max_evals = nevals; }

            out_img_row[jp] = color;
          }

        irt_write_color_image_row(out_img_file, ip, width, out_img_row);
        fprintf(stderr, ".");
        if ((ip % 10 == 0) || (ip == height-1)) fprintf(stderr, "Did line %d.\n", ip);
      }
    double stop_time = user_cpu_time_usec();

    fprintf(stderr, "\n");
    irt_close_color_image_file(out_img_file);
    /* Print time statistics */
    {
      double tsec = (stop_time - start_time)/1.0e6;
      double usec_ray = (stop_time - start_time)/((double)irt_num_rays);
      double usec_eval = (stop_time - start_time)/((double)irt_num_evals);
      fprintf(stderr, "tot rays cast =  %d\n", irt_num_rays);
      fprintf(stderr, "tot proc evals = %d\n", irt_num_evals);
      double avg_evals_per_ray = ((double)irt_num_evals)/((double)irt_num_rays);
      fprintf(stderr, "evals per ray =  %8.1f avg, %d max\n", avg_evals_per_ray, irt_max_evals);
      fprintf(stderr, "\n");
      fprintf(stderr, "tot user time =  %8.1f seconds\n", tsec);
      fprintf(stderr, "time per ray =   %8.1f us\n", usec_ray);
      fprintf(stderr, "time per eval =  %8.1f us\n", usec_eval);
    }

    irt_write_eval_count_image(scene_dir, scene_name, width, height, eval_count);
    fprintf(stderr, "All done.\n");
    return (0);
  }

r3_t irt_pixel_center(view_t *vw, int ip, int jp)
  { r3_t dy, dx;
    r3_t pix = vw->focus;
    r3_scale( + jp - 0.5 * vw->image_width,  &vw->screen_x, &dx);
    r3_scale( - ip + 0.5 * vw->image_height,  &vw->screen_y, &dy);
    r3_add(&dx, &pix, &pix);
    r3_add(&dy, &pix, &pix);
    return(pix);
  }

void irt_write_eval_count_image(char *scene_dir, char *scene_name, int width, int height, int count[])
  { FILE *f = irt_open_count_image_file(scene_dir, scene_name, "-evals", width, height);
    int ip;
    for (ip = 0; ip < height; ip++)
      { int *ctrow = &(count[ip*width]);
        irt_write_count_image_row(f, ip, width, ctrow);
      }
    irt_close_count_image_file(f);
  }

irt_options_t *irt_get_options(int argc, char**argv)
  { 
    irt_options_t *o = (irt_options_t *)notnull(malloc(sizeof(irt_options_t)), "no mem");
    
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
/*     ??? /\* Note that level 40 requires seqs with at least {2^40} bases: *\/  */
/*     ??? if (argparser_keyword_present(pp, "-maxLevel")) */
/*     ???   { o->maxLevel = argparser_get_next_int(pp, 0, 40); } */
/*     ??? else */
/*     ???   { o->maxLevel = 40; }  */
    
/*     ??? argparser_get_keyword(pp, "-minBases"); */
/*     ??? o->minBases = argparser_get_next_int(pp, 0, INT_MAX); */
    
/*     ??? argparser_get_keyword(pp, "-scores"); */
/*     ??? o->sc = dm_score_args_parse(pp); */
    
/*     ??? if (argparser_keyword_present(pp, "-delta")) */
/*     ???   { o->delta = argparser_get_next_int(pp, 0, INT_MAX); } */
/*     ??? else */
/*     ???   { o->delta = 3; } */
    
/*     ??? if (argparser_keyword_present(pp, "-kappa")) */
/*     ???   { o->kappa = argparser_get_next_int(pp, 0, INT_MAX); } */
/*     ??? else */
/*     ???   { o->kappa = 6; } */
    
/*     ??? argparser_get_keyword(pp, "-maxMatches"); */
/*     ??? o->maxMatches = argparser_get_next_int(pp, 0, INT_MAX); */
    
/*     ??? argparser_get_keyword(pp, "-initFilter"); */
/*     ??? o->w0 = dm_weight_args_parse(pp, TRUE); */

/*     ??? argparser_get_keyword(pp, "-incrFilter"); */
/*     ??? o->w1 = dm_weight_args_parse(pp, TRUE); */

/*     ??? if (argparser_keyword_present(pp, "-contFrac")) */
/*     ???   { o->contFrac = argparser_get_next_double(pp, -1.0, +1.0); } */
/*     ??? else */
/*     ???   { o->contFrac = 0.99; } */
    
    argparser_skip_parsed(pp);
    
    o->scene_dir = argparser_get_next(pp);
    o->scene_name = argparser_get_next(pp);
    
/*     ??? o->seqNameA = argparser_get_next(pp); */
/*     ??? o->circA = argparser_get_next_bool(pp); */
    
/*     ??? o->seqNameB = argparser_get_next(pp); */
/*     ??? o->circB = argparser_get_next_bool(pp); */
    
/*     ??? o->outName = argparser_get_next(pp); */
    
/*     ??? o->dontRefine = FALSE; /\* DEBUG *\/ */
    
    argparser_finish(pp);
    
    return o;
  }

void irt_close_postscript_stream(PSStream *ps)
  {
    pswr_close_stream(ps);
  }

PSStream *irt_open_postscript_stream(char *dir, char *name, int col, int row, char *tag)
  {
    double size = 6.0 * 72.0;
    /* double xc = 4.25 * 72.0; */
    /* double hmin = xc - size / 2.0; */
    /* double hmax = xc + size / 2.0; */
    /* double yc = 6.50 *72.0; */
    /* double vmin = yc - size / 2.0; */
    /* double vmax = yc + size / 2.0; */

    char *fname = NULL;
    asprintf(&fname, "%s/out/%s-%04d-%04d%s", dir, name, col, row, tag);
    PSStream *ps = pswr_new_stream(fname, NULL, FALSE, "doc", "letter", FALSE, 0, 0);
    pswr_set_canvas_layout(ps, size, size, FALSE, 0, 0, 0, 1, 1);
    return ps;
  }
