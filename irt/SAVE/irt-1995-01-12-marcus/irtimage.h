#ifndef IRTIMAGE_H
#define IRTIMAGE_H

#include <irtscene.h>

typedef struct {
    int width;
    int height;
    double max_v;
    double min_v;
    float *val;
  } irt_gray_image_t;
  
void irt_write_gray_image(char *image_name, irt_gray_image_t *image);
  /*
    Writes $image$ as <$image_name$>.pgm, after scaling it. */

#endif

