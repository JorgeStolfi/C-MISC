#include <ioprotos.h>
#include <stdio.h>
#include <ia.h>
#include <aa.h>
#include <r3.h>
#include <h3.h>
#include <r3.h>
#include <rgb.h>
#include <irtout.h>
#include <irtscene.h>
#include <irtshade.h>

#define DEBUG 0

/*** INTERNAL PROTOTYES ***/

void main(int argc, char *argv[]);

r3_t irt_pixel_center(view_t *vw, int ip, int jp);

/*** MAIN PROGRAM ***/

void main(int argc, char *argv[])
  {
    scene_t sc;
    FILE *image_file;
    int *evals_image;
    FILE *evals_file;
    int ip, jp;
    double start_time, stop_time;

    h3_point_t eye; /* Homogeneous coordinates of observer */

    affirm(argc == 2, "irt_main: bad parameters\nusage: irt <scene name>");

    ia_init();
    aa_init();
    
    irt_read_scene(argv[1], &sc);
    h3_from_cart(&sc.view.observer, &eye);
    h3_inf_reduce(&eye); /* To improve numerical behavior */

    image_file = irt_open_output_image (sc.name, sc.view.image_width, sc.view.image_height);

    evals_image = (int *) malloc(sizeof(int)*sc.view.image_height*sc.view.image_width); 
     
    start_time = user_cpu_time_usec();
    for (ip = 0; ip < sc.view.image_height; ip++)
    {
      for(jp = 0; jp < sc.view.image_width; jp++)
	{
	  r3_t pix;        /* Pixel center */
	  r3_t dir;        /* Ray's direction */

	  rgb_t color;
          int print_ray, plot_ray, debug_flags;

          irt_evals_image = 0;
          
          print_ray = (ip == sc.print_ray.row) & (jp == sc.print_ray.col);
          plot_ray = (ip == sc.plot_ray.row) & (jp == sc.plot_ray.col);
          debug_flags = (plot_ray ? 2 : 0) | (print_ray ? 1 : 0);

	  pix = irt_pixel_center(&(sc.view), ip, jp);
	  r3_sub(&pix, &(sc.view.observer), &dir);
          (void) r3_normalize(&dir);
          
          if (print_ray)
            { fprintf(stderr, ". . . . . . . . . . . . . . . . . . . . . . . . . . . . . .\n");
              fprintf(stderr, "row = %d col = %d\n", ip, jp);
            }
              
          color = irt_compute_scene_color(&sc, &eye, &dir, sc.looks.max_bounces, debug_flags);

          irt_output_pixel (image_file, &color);
       
/* fprintf(stderr,"evals_image[%d] =  %d \n",(ip*sc.view.image_width) + jp,irt_evals_image);
*/          
          evals_image[(ip*sc.view.image_width) + jp] = irt_evals_image;
	}

	irt_end_image_row (image_file, ip);
	fprintf(stderr, ".");
        if (ip % 10 == 0) fprintf(stderr, "Did line %d.\n", ip);
      }

    irt_close_output_image (image_file);
    stop_time = user_cpu_time_usec();
    /* Print time statistics */
    { 
      double tsec = (stop_time - start_time)/1.0e6;
      double usec_ray = (stop_time - start_time)/((double)irt_num_rays);
      double usec_eval = (stop_time - start_time)/((double)irt_evals_image);
      fprintf(stderr, "tot rays cast =  %d\n", irt_num_rays);
      fprintf(stderr, "tot proc evals = %d\n", irt_evals_image);
      fprintf(stderr, "evals per ray =  %f\n", ((double)irt_evals_image)/((double)irt_num_rays));
      fprintf(stderr, "\n");
      fprintf(stderr, "tot user time =  %f seconds\n", tsec);
      fprintf(stderr, "time per ray =   %f us\n", usec_ray);
      fprintf(stderr, "time per eval =  %f us\n", usec_eval);
    }

    evals_file = fopen (txtcat(sc.name,".evals"),"w+t");
    for (ip = 0; ip < sc.view.image_height; ip++)
    {
      for(jp = 0; jp < sc.view.image_width; jp++)
            fprintf(evals_file,"%4d ",evals_image[(ip*sc.view.image_width) + jp]);
      fprintf(evals_file,"\n");
    }
    fclose(evals_file);

    fprintf(stderr, "All done.\n"); 
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

