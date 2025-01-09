/* See {irt_out.h}. */
/* Last edited on 2023-02-22 20:09:39 by stolfi */

/*  Escreve no formato ppm P6 com o valor maximo de cada
    componente rgb igual a 255  */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <math.h>

#include <affirm.h>
#include <jsfile.h>
#include <jsstring.h>
#include <jspnm.h>
#include <bool.h>
#include <frgb.h>

#include <irt_out.h>

/*** INTERNAL PROTOTYPES ***/

int32_t gammacorrect(double intensity);

/*** IMPLEMENTATIONS ***/

FILE *irt_open_color_image_file(char *outPrefix, char *tag, int32_t width, int32_t height)
  { char *fname = NULL;
    char *fname = jsprintf("%s%s.ppm", outPrefix, tag);
    FILE *image_file = open_write(fname, TRUE);

    /* Image format is is the "raw bits" ppm format with 8 bits/sample. */
    /* See "man 5 ppm" in the Portable Bitmap package (pbmplus) */
    pnm_format_t format = RPPM_FORMAT;
    uint16_t maxval = 255;
    pnm_write_header(image_file, width, height, maxval, format);
    free(fname);
    return (image_file);
  }

void irt_write_color_image_row(FILE *image_file, int32_t irow, int32_t width, frgb_t rgb[])
  { uint16_t maxval = 255;
    uint16_t smp[3*width];
    int32_t i, c;
    for (i = 0; i < width; i++)
      { for (c = 0; c < 3; c++) 
          { smp[i*3 + c] = (uint16_t)gammacorrect(rgb[i].c[c]); }
      }
    pnm_write_pixels(image_file, smp, width, 3, maxval, TRUE, FALSE);
    fflush(image_file);
  }

void irt_close_color_image_file(FILE *image_file)
  {
    fclose(image_file);
  }

FILE *irt_open_count_image_file(char *outPrefix, char *tag, int32_t width, int32_t height)
  {
    char *fname = jsprintf("%s%s.pgm", outPrefix, tag);
    FILE *image_file = open_write(fname, TRUE);
    
    /* Write header for ascii pgm file format */
    /* See "man 5 pgm" in the Portable Bitmap package (pbmplus) */
    
    pnm_format_t format = PGM_FORMAT;
    uint16_t maxval = PNM_FILE_MAX_MAXVAL;
    pnm_write_header(image_file, width, height, maxval, format);
    free(fname);
    return image_file;
  }

void irt_write_count_image_row(FILE *image_file, int32_t irow, int32_t width, int32_t count[])
  { uint16_t maxval = PNM_FILE_MAX_MAXVAL;
    uint16_t smp[width];
    int32_t i;
    for (i = 0; i < width; i++) { smp[i] = (uint16_t)(count[i] < maxval ? count[i] : maxval); }
    pnm_write_pixels(image_file, smp, width, 1, maxval, FALSE, FALSE);
    fflush(image_file);
  }

void irt_close_count_image_file(FILE *image_file)
  {
    fclose(image_file);
  }

int32_t gammacorrect (double y)
  { int32_t ival;
    if(y >= 1.0)
      { ival = 255; }
    else if(y <= 0.0)
      { ival = 0; }
    else
      { y = exp(log(y) / irt_GAMMA);
        ival = (int32_t)(255.0 * y + 0.5);
      }
    return(ival);
  }




