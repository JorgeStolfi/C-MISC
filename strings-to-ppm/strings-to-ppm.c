/* Reads a bunch of colored strings and produces a PPM image of them. */
/* Last edited on 2023-10-15 03:33:43 by stolfi */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <affirm.h>
#include <jsstring.h>

static char *usage = \
  " strings-to-ppm \\\n" \
  "   -numBits NUMBITS \\\n" \
  "   < INFILE > IMG.ppm " ;

/*
  The INFILE must contain records with the format STRING R G B where
  STRING is a string of NUMBITS binary digits, and R G B are
  non-negative real numbers. The same STRING may appear multiple times,
  but the total "R", "G", and "B" of all its entries must be in [0 _ 1].
*/

typedef struct color { float r; float g; float b; } color;
  
/* PROTOTYPES */

int main(int argc, char **argv);
void parse_args(int argc, char **argv, int *numBitsP);
void process_lines(int numBits, color *totColor);
int bin_from_string(char *str);
int pix_from_bin(int bin, int numBits);
void write_ppm(int maxBits, color *totColor);
int conv_int(float v);
void arg_error (char *msg);
void file_error (int NR, char *msg);

static double MAXLOGPIX; /* Set below. */

int main(int argc, char **argv)
  {
    int numBits, numPixels;
    int ip;
    color *totColor;

    parse_args(argc, argv, &numBits);
    
    /* Allocate image: */
    numPixels = 1 << numBits;
    totColor = (color*)malloc(sizeof(color)*numPixels);
    affirm (totColor != NULL, "memory exhausted");
    for(ip=0; ip < numPixels; ip++) { totColor[ip] = (color){0,0,0}; }

    /* Process input lines: */
    process_lines(numBits, totColor);

    /* Write image: */
    MAXLOGPIX = log(255.0 + 1.0);
    write_ppm(numBits, totColor);
    return 0;
  }

void process_lines(int numBits, color *totColor)
  {
    char wd[33];
    char eol;
    float r, g, b;
    int pCode;
    int res;
    color *col;
    int NR = 0;
    while((res = scanf("%32s %f %f %f%c", wd, &r, &g, &b, &eol)) != EOF)
      { NR++;
        while(eol == ' ') { eol = getchar(); }
        /* fprintf(stderr, "%s %f %f %f %o [%d]\n", wd, r, g, b, (int)eol, res); */
        if ((res != 5) || (eol != '\n'))  { file_error(NR, "bad record format"); }
        if (strlen(wd) > numBits)  { file_error(NR, "string too long"); }
        if ((r < 0) || (r > 1))  { file_error(NR, "bad red value"); }
        if ((g < 0) || (g > 1))  { file_error(NR, "bad green value"); }
        if ((b < 0) || (b > 1))  { file_error(NR, "bad blue value"); }
        pCode = pix_from_bin(bin_from_string(wd), numBits);
        if ( pCode == -1 ) { file_error(NR, "invalid char in word"); }
        /*Store in matrix: */
        col = &(totColor[pCode]);
        col->r += r;
        col->g += g;
        col->b += b;
      }
    fprintf(stderr, "%7d lines read\n", NR);
  }

int bin_from_string(char *wd)
  {
    char ch;
    int ip;
    int bin;
    /* Convert binary string to numeric: */
    bin = 0;
    for(ip=0; (ch=wd[ip]) != '\000'; ip++)
      { if (ch == '0') { bin = (bin << 1); }
        else if (ch == '1') { bin = (bin << 1) | 1; }
        else { return -1; }
      }
    return bin;
  } 

int pix_from_bin(int bin, int numBits)
  {
    int halfBits = numBits - (numBits >> 1);
    int k, gCode, xCode, yCode, temp, mask;
    if (bin < 0) { return -1; }
    
    /* Convert code to Gray code: */
    gCode = 0; temp = bin;
    while (temp > 0) { gCode ^= temp; temp >>= 1; }
    /* fprintf(stderr, "%08X = %08X\n", bin, gCode); */
    
    /* Split even/odd: */
    xCode = 0; yCode = 0;
    mask = (1 << halfBits) << halfBits;
    for (k=0; k<halfBits; k++) { 
      mask >>= 1; yCode = (yCode << 1) | ((mask & gCode) ? 1 : 0);
      mask >>= 1; xCode = (xCode << 1) | ((mask & gCode) ? 1 : 0);
    }
    /* fprintf(stderr, "%08X = [%08X,%08X]\n", gCode, yCode, xCode); */
    return (yCode << halfBits) | xCode;
  }

void parse_args(int argc, char **argv, int *numBitsP)
  {
    int numBits = -1;
    int i;
    for (i=1; i<argc; i++) 
      { /* fprintf(stderr, "argc[%d] = %s\n", i, argv[i]); */
        if (strcmp(argv[i],"-numBits") == 0)
          { i++; 
            if (i >= argc) { arg_error("no arg after \"-numBits\""); }
            if (numBits != -1) { arg_error("repeated \"-numBits\""); }
            if (sscanf(argv[i], "%d", &numBits) != 1) 
              { arg_error("bad format in \"-numBits\""); }
            if (numBits < 1)  
              { arg_error("bad \"-numBits\""); }
          }
        else
          { arg_error(txtcat3("bad option \"", argv[i], "\"")); }
      }
    if (numBits == -1) { arg_error("must specify \"-numBits\""); }
    (*numBitsP) = numBits;
  }

void write_ppm(int numBits, color *totColor)
  {
    int NY = 1 << (numBits/2);
    int NX = 1 << (numBits - (numBits/2));
    int numPixels = NX*NY;
    int ip;
    printf("P6\n");
    printf("%d %d 255\n", NX, NY);
    for (ip=0; ip<numPixels; ip++)
      { color *col=&(totColor[ip]);
        putchar((char)conv_int(col->r));
        putchar((char)conv_int(col->g));
        putchar((char)conv_int(col->b));
      }
    fclose(stdout);
  }

int conv_int(float v)
  /* Converts a float intensity to an integer pixel value [0..255]: */
  {
    if (v <= 0.0)
      { return 0; }
    else if (v >= 1.0 )
      { return 255; }
    else
      { v = log(v*255.0+1.0)/MAXLOGPIX;
        return (int)(255*v + 0.5);
      }
  }

void arg_error (char *msg)
  { fprintf (stderr, "*** %s\n", msg);
    fprintf (stderr, "usage %s\n", usage);
    exit(1);
  }
  
void file_error (int NR, char *msg)
  { fprintf (stderr, "*** line %d: %s\n", NR, msg);
    exit(1);
  }

