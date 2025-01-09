#define PROG_NAME "irt_main"
#define PROG_DESC "prototype of an interval ray tracer"
#define PROG_VERS "1.0"

/* Last edited on 2023-02-22 20:13:25 by stolfi */

#define irt_main_C_COPYRIGHT \
  "Copyright © 2005 by the State University of Campinas (UNICAMP)"
  /* See {argparser_help_info_STANDARD_RIGHTS} in {argparser.h}. */
  /* See {argparser_help_info_NO_WARRANTY} in {argparser.h}. */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include <ia.h>
#include <aa.h>
#include <r3.h>
#include <hr3.h>
#include <r3.h>
#include <frgb.h>
#include <epswr.h>
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

typedef struct irt_options_t
  { int32_t imageSize_X;
    int32_t imageSize_Y;
    char *parmsFile;
    char *pcodeFile;
    char *outPrefix;

    /* Rendering options: */
    arith_t arith;            /* Which artihmetic to use in ray-tracing. */
    pixel_num_t printRay;     /* Print a trace of this ray's evaluation. */ 
    pixel_num_t plotRay;      /* Plot the shape function along this ray. */
  } irt_options_t;

/*** INTERNAL PROTOTYES ***/

int32_t main(int32_t argc, char **argv);

irt_options_t *irt_get_options(int32_t argc, char **argv);

r3_t irt_pixel_center(view_t *vw, int32_t ip, int32_t jp);

void irt_write_eval_count_image(char *outPrefix, int32_t width, int32_t height, int32_t count[]);
  /* Writes a ".pgm" image with the evaluation counts. */

epswr_figure_t *irt_open_eps_figure(char *outPrefix, int32_t col, int32_t row, char *tag);
  /* Opens an Encapsulated Postscript file called "{dir}/out/{name}-{col}-{row}{tag}.eps"
    for the pixel {col,row}.  The picture layout has one picture per page
    with equal H and V dimensions. */

void irt_close_eps_figure(epswr_figure_t *epsf);
  /* Closes the Postscript stream {epsf}, properly finishing the current plot. */

/*** MAIN PROGRAM ***/

int32_t main(int32_t argc, char **argv)
  {
    ia_init();
    aa_init();

    irt_options_t *o = irt_get_options(argc, argv);
    
    /* Read the scene description {sc}: */
    scene_t sc;
    irt_read_scene(o->pcodeFile, o->parmsFile, &sc);

    double pix_size = 0.3; /* Assumed screen pixel size (mm) */
    double D = 500;        /* Assumed screen viewing distance (mm) */
    int32_t width = o->imageSize_X;
    int32_t height = o->imageSize_Y;
    irt_set_image_parameters(&sc, width, height, pix_size, D);
    
    /* Get homogeneous coordinates {eye} of observer: */
    hr3_point_t eye = hr3_from_r3(&sc.view.observer);
    hr3_L_inf_normalize_point(&eye); /* To improve the numerics. */

    /* Open the output image file: */
    FILE *out_img_file = irt_open_color_image_file(o->outPrefix, "", width, height);
    frgb_t out_img_row[width];

    /* Allocate the array {eval_count} of function evaluation counts per pixel: */
    int32_t *eval_count = (int32_t *)malloc(sizeof(int32_t)*width*height);
    int32_t irt_max_evals = 0; /* Mex eval count per ray. */

    double start_time = user_cpu_time_usec();
    int32_t prev_evals;
    int32_t ip, jp;
    for (ip = 0; ip < height; ip++)
      { for(jp = 0; jp < width; jp++)
          { r3_t pix;        /* Pixel center */
            r3_t dir;        /* Ray's direction */

            frgb_t color;
            int32_t print_ray, plot_ray;

            print_ray = (ip == o->printRay.row) & (jp == o->printRay.col);
            plot_ray = (ip == o->plotRay.row) & (jp == o->plotRay.col);

            pix = irt_pixel_center(&(sc.view), ip, jp);
            r3_sub(&pix, &(sc.view.observer), &dir);
            (void) r3_dir(&dir, &dir);

            if (print_ray)
              { fprintf(stderr, ". . . . . . . . . . . . . . . . . . . . . . . . . . . . . .\n");
                fprintf(stderr, " col = %d row = %d\n", jp, ip);
              }
              
            epswr_figure_t *epsf = 
              ( plot_ray ? 
                irt_open_eps_figure(o->outPrefix, jp, ip, "-debug") : 
                NULL
              );

            prev_evals = irt_num_evals;

            color = irt_compute_scene_color(&sc, &eye, &dir, sc.looks.max_bounces, o->arith, print_ray, epsf);
            
            if (plot_ray) { irt_close_eps_figure(epsf); }

            int32_t nevals = irt_num_evals - prev_evals; /* Count of evals for this ray. */
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

    irt_write_eval_count_image(o->outPrefix, width, height, eval_count);
    fprintf(stderr, "All done.\n");
    return (0);
  }

r3_t irt_pixel_center(view_t *vw, int32_t ip, int32_t jp)
  { r3_t dy, dx;
    r3_t pix = vw->focus;
    r3_scale( + jp - 0.5 * vw->image_width,  &vw->screen_x, &dx);
    r3_scale( - ip + 0.5 * vw->image_height,  &vw->screen_y, &dy);
    r3_add(&dx, &pix, &pix);
    r3_add(&dy, &pix, &pix);
    return(pix);
  }

void irt_write_eval_count_image(char *outPrefix, int32_t width, int32_t height, int32_t count[])
  { FILE *f = irt_open_count_image_file(outPrefix, "-evals", width, height);
    int32_t ip;
    for (ip = 0; ip < height; ip++)
      { int32_t *ctrow = &(count[ip*width]);
        irt_write_count_image_row(f, ip, width, ctrow);
      }
    irt_close_count_image_file(f);
  }

#define irt_MAX_IMAGE_SIZE 20000
  /* Max image width and height, just paranoia. */

irt_options_t *irt_get_options(int32_t argc, char**argv)
  { 
    irt_options_t *o = (irt_options_t *)notnull(malloc(sizeof(irt_options_t)), "no mem");
    
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    argparser_get_keyword(pp, "-pcodeFile");
    o->pcodeFile = argparser_get_next_non_keyword(pp);
    
    argparser_get_keyword(pp, "-parmsFile");
    o->parmsFile = argparser_get_next_non_keyword(pp);

    argparser_get_keyword(pp, "-imageSize");
    o->imageSize_X = (int32_t)argparser_get_next_int(pp, 0, irt_MAX_IMAGE_SIZE);
    o->imageSize_Y = (int32_t)argparser_get_next_int(pp, 0, irt_MAX_IMAGE_SIZE);

    /* Rendering options: */
    argparser_get_keyword(pp, "-arith");
    char *arcode = argparser_get_next_non_keyword(pp);
    if (strcmp(arcode, "ia") == 0)
      { o->arith = arith_ia; }
    else if (strcmp(arcode, "id") == 0)
      { o->arith = arith_ia_diff; }
    else if  (strcmp(arcode, "aa") == 0)
      { o->arith = arith_aa; }
    else if  (strcmp(arcode, "mx") == 0)
      { o->arith = arith_mix; }
    else 
      { argparser_error(pp, "invalid arithmetic code"); }
      
    if (argparser_keyword_present(pp, "-printRay"))
      { o->printRay.col = (int32_t)argparser_get_next_int(pp, -1, irt_MAX_IMAGE_SIZE-1);
        o->printRay.row = (int32_t)argparser_get_next_int(pp, -1, irt_MAX_IMAGE_SIZE-1);
      }
    else
      { o->printRay = (pixel_num_t){ -1, -1 }; }

    if (argparser_keyword_present(pp, "-plotRay"))
      { o->plotRay.col = (int32_t)argparser_get_next_int(pp, -1, irt_MAX_IMAGE_SIZE-1);
        o->plotRay.row = (int32_t)argparser_get_next_int(pp, -1, irt_MAX_IMAGE_SIZE-1);
      }
    else
      { o->plotRay = (pixel_num_t){ -1, -1 }; }
    
    argparser_get_keyword(pp, "-outPrefix");
    o->outPrefix = argparser_get_next_non_keyword(pp);
    
    argparser_skip_parsed(pp);
   
    argparser_finish(pp);
    
    return o;
  }

void irt_close_eps_figure(epswr_figure_t *epsf)
  {
    epswr_end_figure(epsf);
  }

epswr_figure_t *irt_open_eps_figure(char *outPrefix, int32_t col, int32_t row, char *tag)
  {
    double maxSize = 160*epswr_pt_per_mm;
    
    double xmin = 0.0, xmax = 1.0;
    double ymin = 0.0, ymax = 1.0;

    char *fname = jsprintf("%s-%04d-%04d%s", outPrefix, col, row, tag);
    int32_t capLines = 0; 
    double capFontHeight = 10;
    bool_t epsf_verbose = FALSE;
    epswr_figure_t *epsf = epswr_new_captioned_figure
      ( NULL, fname, NULL, -1, NULL,
        xmin, xmax, ymin, ymax, maxSize, maxSize,
        capLines, capFontHeight, epsf_verbose
      );
    free(fname);
    return epsf;
  }
