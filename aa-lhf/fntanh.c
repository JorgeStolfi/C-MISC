/* Last edited on 2023-01-14 12:23:34 by stolfi */

#include <memory.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <arith.h>
#include <fn.h>

/* PROTOTYPES */

char *fnname(void)
{ 
  return("fntanh"); 
} 
  
void fndescr(FILE *f)
{
  fprintf(f, "GIVEN (x, y, z) DO\n");
  fprintf(f, "  f = tanh(x + y) - z;\n");
  fprintf(f, "RETURN (f)\n");
}

void fneval(AAform x, AAform y, AAform z, AAform f)
{
  aa_add(f, x, y);
  aa_tanh(f, f);
  aa_sub(f, f, z);
}
