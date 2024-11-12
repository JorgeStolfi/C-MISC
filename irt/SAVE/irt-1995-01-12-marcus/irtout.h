#ifndef IRTOUT_H
#define IRTOUT_H

#include <rgb.h>
#include <stdio.h>

FILE *irt_open_output_image (char *scene_name, int width, int height);

void irt_output_pixel (FILE *image_file, rgb_t *rgb);

void irt_end_image_row (FILE *image_file, int row);

void irt_close_output_image (FILE *image_file);

#endif
