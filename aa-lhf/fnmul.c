/* Last edited on 2023-01-14 12:22:49 by stolfi */

#include <memory.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <arith.h>
#include <fn.h>

/* PROTOTYPES */

char *fnname(void)
{ 
  return("fnmul"); 
} 
  
void fndescr(FILE *f)
{
  fprintf(f, "GIVEN (x, y, z) DO\n");
  fprintf(f, "  f = x * y - z;\n");
  fprintf(f, "RETURN (f)\n");
}

void fneval(AAform x, AAform y, AAform z, AAform f)
{
  aa_mul(f, x, y);
  aa_sub(f, f, z);
}
