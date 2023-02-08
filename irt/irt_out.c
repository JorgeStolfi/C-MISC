/* See {irt_out.h}. */
/* Last edited on 2017-02-26 02:50:47 by stolfilocal */

/*  Escreve no formato ppm P6 com o valor maximo de cada
    componente rgb igual a 255  */

#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>

#include <affirm.h>
#include <jsfile.h>
#include <jsstring.h>
#include <jspnm.h>
#include <bool.h>
#include <frgb.h>

#include <irt_out.h>


/*** INTERNAL PROTOTYPES ***/

int gammacorrect(double intensity);

/*** IMPLEMENTATIONS ***/

FILE *irt_open_color_image_file(char *dir_name, char *scene_name, char *tag, int width, int height)
  { char *fname = NULL;
    asprintf(&fname, "%s/out/%s%s.ppm", dir_name, scene_name, tag);
    FILE *image_file = open_write(fname, TRUE);

    /* Image format is is the "raw bits" ppm format with 8 bits/sample. */
    /* See "man 5 ppm" in the Portable Bitmap package (pbmplus) */
    pnm_format_t format = RPPM_FORMAT;
    uint16_t maxval = 255;
    pnm_write_header(image_file, width, height, maxval, format);
    free(fname);
    return (image_file);
  }

void irt_write_color_image_row(FILE *image_file, int irow, int width, frgb_t rgb[])
  { uint16_t maxval = 255;
    uint16_t smp[3*width];
    int i, c;
    for (i = 0; i < width; i++)
      { for (c = 0; c < 3; c++) 
          { smp[i*3 + c] = gammacorrect(rgb[i].c[c]); }
      }
    pnm_write_pixels(image_file, smp, width, 3, maxval, TRUE, FALSE);
    fflush(image_file);
  }

void irt_close_color_image_file(FILE *image_file)
  {
    fclose(image_file);
  }

FILE *irt_open_count_image_file(char *dir_name, char *scene_name, char *tag, int width, int height)
  {
    char *fname = NULL;
    asprintf(&fname, "%s/out/%s%s.pgm", dir_name, scene_name, tag);
    FILE *image_file = open_write(fname, TRUE);
    
    /* Write header for ascii pgm file format */
    /* See "man 5 pgm" in the Portable Bitmap package (pbmplus) */
    
    pnm_format_t format = PGM_FORMAT;
    uint16_t maxval = PNM_FILE_MAX_MAXVAL;
    pnm_write_header(image_file, width, height, maxval, format);
    free(fname);
    return image_file;
  }

void irt_write_count_image_row(FILE *image_file, int irow, int width, int count[])
  { uint16_t maxval = PNM_FILE_MAX_MAXVAL;
    uint16_t smp[width];
    int i;
    for (i = 0; i < width; i++) { smp[i] = (count[i] < maxval ? count[i] : maxval); }
    pnm_write_pixels(image_file, smp, width, 1, maxval, FALSE, FALSE);
    fflush(image_file);
  }

void irt_close_count_image_file(FILE *image_file)
  {
    fclose(image_file);
  }

int gammacorrect (double y)
  { int ival;
    if(y >= 1.0)
      { ival = 255; }
    else if(y <= 0.0)
      { ival = 0; }
    else
      { y = exp(log(y) / irt_GAMMA);
        ival = (int)(255.0 * y + 0.5);
      }
    return(ival);
  }




