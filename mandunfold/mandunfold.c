#define PROG_NAME "mandunfold"
#define PROG_DESC "generates unfolded images of the Mandelbrot set"
#define PROG_VERS "1.0"

#define mandunfold_C_COPYRIGHT "Copyright © 2008 by the State University of Campinas (UNICAMP)"

/* Last edited on 2024-12-21 11:55:40 by stolfi */

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    -func {FUNC}.. \\\n" \
  "    -niter {NITER} \\\n" \
  "    -size {SIZE} \\\n" \
  "    -domain {CTR_RE} {CTR_IM} {RAD_RE} {RAD_IM} \\\n" \
  "    -range {CTR_RE} {CTR_IM} {RAD_RE} {RAD_IM} \\\n" \
  "    -outName {OUT_NAME} \\\n" \
  "    " argparser_help_info_HELP ""

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  The program generates a PGM image file \"{OUT_NAME}-W.pgm\" containing a picture" \
  " of a set {W} such that {h(W)} is a subset of the Mandelbrot set;" \
  " where {h} is a function identified by the {FUNC}" \
  " parameter.  It also generates a PGM image" \
  " file \"{OUT_NAME}-M.pgm\" containing a picture" \
  " of the set {h(W)}.\n" \
  "\n" \
  "  The set {W} consists of all complex numbers" \
  " {c}, in a specified rectangle of the complex plane," \
  " such that the iteration {z := z^2 + h(c)}, starting" \
  " with {z := 0}, does not diverge after {NITER}" \
  " iterations.  The iteration is assumed to" \
  " diverge if {\abs{z}} becomes greater than 2 at any point.\n" \
  "\n" \
  "OPTIONS\n" \
  "  -func {FUNC}\n" \
  "    This mandatory parameter specifies the numeric identifier {FUNC}" \
  " of the map {h}.  Function 0 is the identity.\n" \
  "\n" \
  "  -niter {NITER}\n" \
  "    This mandatory argument specifies the number" \
  " of iterations to be tried before deciding whether {c} is in the set.\n" \
  "\n" \
  "  -size {SIZE}\n" \
  "    This mandatory argument specifies that the output" \
  " image may have at most {SIZE} columns and {SIZE} rows." \
  " .\n" \
  "\n" \
  "  -domain {CTR_RE} {CTR_IM} {RAD_RE} {RAD_IM}\n" \
  "  -range {CTR_RE} {CTR_IM} {RAD_RE} {RAD_IM}\n" \
  "    These optional arguments specify the plot window" \
  " for the domain image (showing the set {W}) and range image" \
  " (showing the set {h(W)}, a part of {M}).  The nominal plot window is" \
  " the smallest axis-aligned rectangle in the complex plane that with corners" \
  " {(CTR_RE±RAD_RE) + I*(CTR_IM±RAD_IM)}.  (The actual" \
  " plot window may be slighltly larger, since it must" \
  " contain an integer number of square pixels.)  If this" \
  " parameter is omitted, or the given rectangle is" \
  " empty, uses a default rectangle.  The default" \
  " for {W} depends on {FUNC}; the default for" \
  " {M} is {[-2_+2]×[-2_+2]}. \n" \
  "\n" \
  "  -outName {OUT_NAME}\n" \
  "    This mandatory argument specifies the common prefix {OUT_NAME}" \
  " for all output image files.\n" \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  display(1), pgm(5).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created 2008-05-21 by Jorge Stolfi, IC-UNICAMP.\n" \
  "  Based on {m2.lua} created on 2008-05-20 by Luiz Henrique Figueiredo, IMPA.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  None yet.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " mandunfold_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include <bool.h>
#include <argparser.h>
#include <jsfile.h>

#include <munf_functions.h>

/* COMMAND-LINE OPTIONS */

typedef struct options_t
  { int func;                   /* Code of function. */
    int niter;                  /* Number of iterations to perform. */
    double complex Wctr, Wrad;  /* Center and radius of area of interest in {W} image. */
    double complex Mctr, Mrad;  /* Center and radius of area of interest in {M} image. */
    int size;                   /* Max count of image pixels in either direction. */
    char *outName;
  } options_t;

/* INTERNAL PROTOTYPES */

int main(int argc,char** argv);

options_t *parse_options(int argc, char **argv);
  /* Parses the command line arguments and packs them as an {options_t}. */

void parse_options_image(argparser_t *pp, double complex *ctr, double complex *rad);
  /* Parses the command line arguments after a "-domain" or "-range" keyword. */
  
void compute_image_params
  ( double complex rad, 
    int size,
    double *step, 
    int *ncols, 
    int *nrows
  );
  /* Given the radius (half-diagnal) {rad} of the desired plot window
    of an image, and the maximum {size} for the row and column counts,
    computes the raster parameters --- namely, the pixel size {*step}
    and the column and row counts {*ncols,*nrows}.
    
    The image dimensions will always be odd so that the center of the
    center pixel is the user-specified center point. The actual plot
    window may be slightly larger than the one implied by {rad}, since
    the image must contain an integer number of square pixels. */

int compute_image_size(double rad, double step, int max_size);
  /* Computes the size of an image along the {Re} or {Im} axis,
    given the desired half-extent {rad} of the plot window in that 
    direction, and the pixel size {step}.  
    
    The result {size} will be odd, and will not exceed {max_size}
    (which must be odd). If {size} is less than {max_size}, the
    implied half-extent {size*step/2} may be larger than {rad},
    as much as {rad+step}. */

int in_mandelbrot(double complex h, int niter, bool_t debug);
  /* Tests whether {h} is in the Mandelbrot set, with
    {niter} iterations. If it is, returns {-1}. Otherwise,
    returns the smallest {n} in {1..niter} such that iteration
    {n} produces a point outside the disk of radius 2. */
    
#define GAMMA (2.2)
  /* Power law exponent for intensity encoding. */

int pixel_from_lum(double lum);
  /* Converts a linear intensity {lum} in {[0_1]} to an integer pixel
    value in {0..PGM_MAXVAL}, with gamma correction. */

void bump_mandelbrot_image
  ( uint8_t *hit, 
    double complex ctr, 
    double step, 
    int ncols, 
    int nrows, 
    double complex h
  );
  /* If {h} falls inside the rectangle defined by {ctr,step,ncols,nrows},
    increment the count {hit[ire + iim*ncols]} of the pixel 
    {ire,iim} that contains {h}.  The count saturates at 255. */
  
/* IMPLEMENTATIONS */

#define PGM_MAXVAL 255
/* Maxval of output PGM image. */

int main(int argc, char** argv)
  {
    options_t *o = parse_options(argc, argv);
    
    /* Get h-function's attributes: */
    complex_map_t *h_proc;
    char *h_desc;
    double complex Wctr_default, Wrad_default;
    munf_get_function(o->func, &h_proc, &h_desc, &Wctr_default, &Wrad_default);
    fprintf(stderr, "function h%02d = {c -> %s}\n", o->func, h_desc);
    
    /* Compute parameters of {W} image: */
    if ((creal(o->Wctr - o->Wrad) == 0) || (cimag(o->Wctr - o->Wrad) == 0))
      { /* Given plot window is empty, provide plot window depending on {func}: */
        o->Wctr = Wctr_default;
        o->Wrad = Wrad_default;
      }
    double Wstep;        /* Pixel size of {W} image in {Re} and {Im}. */
    int Wncols, Wnrows;      /* Number of pixels of {W} image along each axis. */
    compute_image_params(o->Wrad, o->size, &Wstep, &Wncols, &Wnrows);
    
    /* Compute parameters of {M} image: */
    if ((creal(o->Mctr - o->Mrad) == 0) || (cimag(o->Mctr - o->Mrad) == 0))
      { /* Given plot window is empty, provide default plot window: */
        o->Wctr = (00.00000) + (00.00000)*I;
        o->Wrad = (02.01000) + (02.01000)*I;
      }
    double Mstep;        /* Pixel size of {M} image in {Re} and {Im}. */
    int Mncols, Mnrows;      /* Number of pixels of {M} image along each axis. */
    compute_image_params(o->Mrad, o->size, &Mstep, &Mncols, &Mnrows);
    
    /* Open the {W} file and output the PGM header: */
    char *Wname = NULL;
    char *Wname = jsprintf("%s-W.pgm", o->outName);
    FILE *Wfile = open_write(Wname, TRUE);
    fprintf(Wfile, "P2\n");
    fprintf(Wfile, "%d %d %d\n", Wncols, Wnrows, PGM_MAXVAL);
    
    /* Allocate the range hit array, clean it out: */
    int Mnpixs = Mncols*Mnrows;
    uint8_t *Mhit; /* {Mhit[Mire + Miim*Mncols]} counts range pixels that mapped to {Mire,Miim}. */
    Mhit = (uint8_t *)notnull(malloc(Mnpixs*sizeof(uint8_t)), "no mem");
    int i;
    for (i = 0; i < Mnpixs; i++) { Mhit[i] = 0; }
    
    /* Scan the domain rectangle, paint the {W} image: */
    int Wcol, Wrow; /* Row and column indices in {W} image. */
    for (Wrow = 0; Wrow < Wnrows; Wrow++) 
      { for (Wcol = 0; Wcol < Wncols; Wcol++)
          { /* bool_t debug = ((Wcol == Wncols/2) && (Wrow == Wnrows/2)); */
            bool_t debug = FALSE;
            
            /* Compute the pixel center coords {c}: */
            int Wdre = + (Wcol - (Wncols/2)); /* {Re} displacement from center pixel. */ 
            int Wdim = - (Wrow - (Wnrows/2)); /* {Im} displacement from center pixel. */ 
            double complex c = o->Wctr + Wstep*(Wdre + I*Wdim);

            /* Compute {h = h(c)}: */
            double complex h = h_proc(c);
            
            /* Perform the Mandelbrot test with {h}: */
            int nit = in_mandelbrot(h, o->niter, debug);
           
            if (debug) 
              { fprintf(stderr, "  c = ( %+15.12f %+15.12f )", creal(c), cimag(c)); 
                fprintf(stderr, "  h = ( %+15.12f %+15.12f )", creal(h), cimag(h));
                fprintf(stderr, "  nit = %5d\n", nit);
              }

            /* Decide pixel value {Xpix} of the {X} image: */
            double lum;
            if (nit < 0) 
              { /* Apparently {h(c)} is in Mandelbrot: */
                lum = 0;
              }
            else
              { /* Relative number of iterations before nominal divergence: */
                double relnit = ((double)nit)/((double)o->niter);
                /* Map to [0.5 _ 1.0]: */
                lum = 0.5 + 0.5*relnit;
              }
            /* Write {W} pixel: */
            fprintf(Wfile, "%d\n", pixel_from_lum(lum));
            
            /* If in mandelbrot, bump the {M} image pixel: */
            if (nit < 0) { bump_mandelbrot_image(Mhit, o->Mctr, Mstep, Mncols, Mnrows, h); }
          }
          
        /* Extra blank line between image rows: */
        fprintf(Wfile, "\n");
        fprintf(stderr, ".");
      }
    fprintf(stderr, "\n");  
    
    /* Close {W} file: */
    fclose(Wfile); free(Wname);
    
    /* Open the {M} file and output the PGM header: */
    char *Mname = NULL;
    char *Mname = jsprintf("%s-M.pgm", o->outName);
    FILE *Mfile = open_write(Mname, TRUE);
    fprintf(Mfile, "P2\n");
    fprintf(Mfile, "%d %d %d\n", Mncols, Mnrows, PGM_MAXVAL);
    
    /* Scan the range rectangle, paint the {M} image: */
    int Mcol, Mrow;
    for (Mrow = 0; Mrow < Mnrows; Mrow++) 
      { for (Mcol = 0; Mcol < Mncols; Mcol++)
          { /* bool_t debug = ((Mcol == Mncols/2) && (Mrow == Mnrows/2)); */
            bool_t debug = FALSE;
            
            /* Compute the pixel center coords {h}: */
            int Mdre = + (Mcol - (Mncols/2)); /* {Re} displacement from center pixel. */ 
            int Mdim = - (Mrow - (Mnrows/2)); /* {Im} displacement from center pixel. */ 
            double complex h = o->Mctr + Mstep*(Mdre + I*Mdim);

            /* Perform the Mandelbrot test with {h}: */
            int nit = in_mandelbrot(h, o->niter, debug);
           
            if (debug) 
              { fprintf(stderr, "  h = ( %+15.12f %+15.12f )", creal(h), cimag(h));
                fprintf(stderr, "  nit = %5d\n", nit);
              }

            /* Get log-count {Mdens} of domain pixels that mapped to this range pixel: */
            int hits = Mhit[Mcol + Mrow*Mncols];
            double Mdens = (hits == 0 ? 0 : 0.5 + 0.5*log(1 + hits)/log(256));
            
            /* Decide pixel value {Mpix} of the {M} image: */
            double lum;
            if (nit < 0) 
              { /* Apparently {h} is in Mandelbrot: */
                lum = 0.000 + 0.333*Mdens;
              }
            else
              { /* Apparently {h} is not in Mandelbrot: */
                lum = 0.667 + 0.333*Mdens;
              }
            /* Write {M} pixel: */
            fprintf(Mfile, "%d\n", pixel_from_lum(lum));
          }
          
        /* Extra blank line between image rows: */
        fprintf(Mfile, "\n");
        fprintf(stderr, ".");
      }
    fprintf(stderr, "\n");  
       
    /* Close {M} file: */
    fclose(Mfile); free(Mname);
    free(Mhit);

    return 0;

  }     

int in_mandelbrot(double complex h, int niter, bool_t debug)
  { int iter = 0;
    double complex z = 0;
    while (iter < niter)
      { z = z*z + h;
        iter++;
        double zabs2 = creal(z)*creal(z) + cimag(z)*cimag(z);
        if (debug) 
          { fprintf(stderr, "    iter = %5d\n", iter);
            fprintf(stderr, "  z = ( %+15.12f %+15.12f )", creal(z), cimag(z)); 
            fprintf(stderr, "  abs = %15.12f", zabs2);
          }
        if (zabs2 >= 4) { return iter; }
      }
    return -1;
  }

void bump_mandelbrot_image
  ( uint8_t *hit, 
    double complex ctr, 
    double step, 
    int ncols, 
    int nrows, 
    double complex h
  )
  { /* Check whether {h} is inside the image's domain: */
    double complex bot = ctr - (ncols + I*nrows)*step/2;
    double complex top = ctr + (ncols + I*nrows)*step/2;
    if ((creal(h) < creal(bot))|| (creal(h) > creal(top))) { return; }
    if ((cimag(h) < cimag(bot))|| (cimag(h) > cimag(top))) { return; }
    /* Compute the indices {col,row} of the pixel that contains {h}: */
    int col = (ncols/2) + (int)floor((creal(h) - creal(ctr))/step + 0.5);
    if (col < 0) { col = 0; } else if (col >= ncols) { col = ncols - 1; }
    int row = (nrows/2) + (int)floor((cimag(ctr) - cimag(h))/step + 0.5);
    if (row < 0) { row = 0; } else if (row >= nrows) { row = nrows - 1; }
    /* Increment that element of {hit}, saturating at 255: */
    int i = col + row*ncols;
    int v = hit[i]; 
    if (v < 255) { v++; hit[i] = v; }
  }

int pixel_from_lum(double lum)
  { /* Apply gamma correction: */
    double fpix = (lum <= 0 ? 0 : (lum >= 1 ? 1 : pow(lum, 1/GAMMA)));
    /* Convert to integer: */
    return (int)(fpix*PGM_MAXVAL + 0.5);
  }

void compute_image_params
  ( double complex rad, 
    int size,
    double *step, 
    int *ncols, 
    int *nrows
  )
  { /* Round {size} down to an odd number: */
    size = size - 1 + (size%2);
    /* Compute the pixel size {step}: */
    double rre = fabs(creal(rad));
    double rim = fabs(cimag(rad));
    (*step) = (2*fmax(rre, rim))/size;
    
    /* Compute the image size {ncols,nrows}: */
    (*ncols) = compute_image_size(rre, (*step), size);
    (*nrows) = compute_image_size(rim, (*step), size);
  }
  
int compute_image_size(double rad, double step, int max_size)
  {
    assert(max_size > 0);
    assert(max_size % 2 == 1);
    int hsize = (int)ceil(rad/step - 0.5); assert(hsize > 0);
    int size = 2*hsize + 1;
    if (size > max_size) { size = max_size; }
    return size;
  }

#define MAX_NITER 1000000
/* Max value of "-niter" parameter. */

#define MAX_SIZE 4096
/* Max value of "-size" parameter. */

#define MAX_CENTER 10.0
/* Max absolute coordinate of center of region of interest. */

#define MAX_RADIUS 10.0
/* Max absolute {re} or {Im} half-extent of region of interest. */

options_t *parse_options(int argc, char **argv)
  {
    /* Initialize argument parser: */
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    /* Allocate the command line argument record: */
    options_t *o = (options_t *)malloc(sizeof(options_t)); 
  
    /* Parse keyword parameters: */
    
    /* Code of {h}: */
    argparser_get_keyword(pp, "-func");
    o->func = argparser_get_next_int(pp, 0, MAX_NITER);
    
    /* Number of iterations: */
    argparser_get_keyword(pp, "-niter");
    o->niter = argparser_get_next_int(pp, 0, MAX_NITER);

    /* Max image dimension: */
    argparser_get_keyword(pp, "-size");
    o->size = argparser_get_next_int(pp, 1, MAX_SIZE);

    /* Domain image parameters: */
    argparser_get_keyword(pp, "-domain");
    parse_options_image(pp, &(o->Wctr), &(o->Wrad));

    /* Range image parameters: */
    argparser_get_keyword(pp, "-range");
    parse_options_image(pp, &(o->Mctr), &(o->Mrad));

    /* Output file name prefix: */
    argparser_get_keyword(pp, "-outName");
    o->outName = argparser_get_next(pp);

    /* Parse positional arguments: */
    argparser_skip_parsed(pp);
    
    /* Check for spurious arguments: */
    argparser_finish(pp);
    
    return o;
  }

void parse_options_image(argparser_t *pp, double complex *ctr, double complex *rad)
  {
    double ctr_re = argparser_get_next_double(pp, -MAX_CENTER, +MAX_CENTER);
    double ctr_im = argparser_get_next_double(pp, -MAX_CENTER, +MAX_CENTER);
    double rad_re = argparser_get_next_double(pp, -MAX_RADIUS, +MAX_RADIUS);
    double rad_im = argparser_get_next_double(pp, -MAX_RADIUS, +MAX_RADIUS);
    (*ctr) = ctr_re + I*ctr_im;
    (*rad) = rad_re + I*rad_im;
  }

