#ifndef irt_out_H
#define irt_out_H

/* Last edited on 2017-02-26 02:50:21 by stolfilocal */
/* Image output facilities for the Interval Ray Tracer. */

#include <frgb.h>
#include <stdio.h>

#define irt_GAMMA 2.4

/* IMAGE FILES

  This interface offers procedures that open, write and close disk
  files containing images. The files are named
  "{dir_name}/out/{scene_name}{tag}.{ext}" where {ext} depends on the
  image type. */

/* COLOR IMAGES
  
  Color image files are in the "raw" PPM 8-bits-per-sample format,
  with {maxval = 255}. The pixels are represented in memory as
  {frgb_t} values in the {[0_1]} scale. They are gamma-corrected when
  written, with exponent {irt_GAMMA}. */

FILE *irt_open_color_image_file(char *dir_name, char *scene_name, char *tag, int width, int height);
void irt_write_color_image_row(FILE *image_file, int irow, int width, frgb_t rgb[]);
void irt_close_color_image_file(FILE *image_file);

/* COUNT IMAGES
  
  Count image files are in the "plain" PGM format with {maxval =
  65535}. The pixels are internally represented as {int} counts.
  The counts are clipped to the range {0..65535} and are written
  out without any rescaling or gamma correction. */

FILE *irt_open_count_image_file(char *dir_name, char *scene_name, char *tag, int width, int height);
void irt_write_count_image_row(FILE *image_file, int irow, int width, int count[]);
void irt_close_count_image_file(FILE *image_file);

#endif
