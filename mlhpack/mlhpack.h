/* Common declarations for mlhpack.c and mlhunpk.c */
/* Last edited on 2024-12-21 11:55:36 by stolfi */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define false 0
#define true 1 

typedef int bool_t;
typedef unsigned char byte; 
typedef signed char tiny;

#ifdef MSDOS 
/* #include <process.h> */
#endif

/* Internal prototypes: */

int main(int argc, char *argv[]);
void setup_digit_tables(void);
void check(bool_t cond, char *msg);

/* Max field length and number of fields: */
#define NC 35
#define NF 20
  
static char digit[NC];   /* Encodes ints in [0..NC-1] as chars */
static tiny digval[256];  /* converts digit chars back to ints */

long num_read = 0;       /* Counts bytes read from stdin */

void setup_digit_tables(void)
  /* digit table setup: */
  { int i;
    for (i=0; i<256; i++) digval[i] = -1;
    for (i=0; i<10; i++) digit[i] = '0' + i;
    for (i=10; i<36; i++) digit[i] = 'a' + (i-10);
    for (i=0; i<36; i++) digval[((byte)digit[i])] = i;
  }
      
void check(bool_t cond, char *msg)
  { if(!cond) 
      { fprintf(stderr, "** byte %ld: %s\n", num_read, msg);
        exit(1);
      }
  }
 
