
/* General interface for the target functions */

void fneval(AAform x, AAform y, AAform z, AAform f);
/* 
  Evaluates the target function */
  
char *fnname(void);
/*
  Returns name of function */
  
void fndescr(FILE *f);
/*
  Prints the function's formula */


