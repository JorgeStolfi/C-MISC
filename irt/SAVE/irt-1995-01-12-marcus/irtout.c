/* Last edited on 2007-01-04 00:19:08 by stolfi */

/*  Escreve no formato ppm P6 com o valor maximo de cada
    componente rgb igual a 255  */

#include "irtout.h"
#include <ioprotos.h>
#include <stdio.h>
#include <math.h>
#include <rgb.h>

/*** INTERNAL PROTOTYPES ***/

int gammacorrect(double intensity);

/*** IMPLEMENTATIONS ***/

FILE *irt_open_output_image (char *scene_name, int width, int height)
  {
    FILE *image_file = open_write(txtcat(scene_name, ".ppm"), TRUE);

    /* Image format is is the "raw bits" ppm format. */
    /* See "man 5 ppm" in the Portable Bitmap package (pbmplus) */
    
    fprintf(image_file, "P6\n%d %d\n255\n", width, height);
    return (image_file);
  }

void irt_output_pixel (FILE *image_file, rgb_t *rgb)
  {
    fputc(gammacorrect(rgb->c[0]), image_file);
    fputc(gammacorrect(rgb->c[1]), image_file);
    fputc(gammacorrect(rgb->c[2]), image_file);
  }

void irt_end_image_row (FILE *image_file, int row)
  {
    fflush(image_file);
  }

void irt_close_output_image (FILE *image_file)
  {
    fclose(image_file);
  }

int gammacorrect (double y)
  {
    #define GAMMA 1.8
    int ival;

    if(y >= 1.0)
      { ival = 255; }
    else if(y <= 0.0)
      { ival = 0; }
    else
      { y = exp(log(y) / GAMMA);
        ival = (int)(255.0 * y + 0.5);
      }
    return(ival);
    #undef GAMMA
  }




