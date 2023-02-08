/* Last edited on 2006-03-19 12:46:13 by stolfi */

#include <SOTentFunction.h>
#include <stdio.h>

int main (int argc, char **argv)
{
  SOTentFunction *tf1, *tf2;
  char arqt1[20], arqt2[20];
  double grad_dot;
 
  printf(" #>> arq tenda1: ");
  scanf("%s", arqt1);
  printf("\n");
  printf(" DBG: ** %s ** \n\n", arqt1);

  FILE *file_t1 = open_read(arqt1, TRUE);

  printf(" #>> arq tenda2: ");
  scanf("%s", arqt2);
  printf("\n");
  printf(" DBG: ** %s ** \n\n", arqt2);

  FILE *file_t2 = open_read(arqt2, TRUE);

  tf1 = SOTentFunction_Read(file_t1, 2);
  tf2 = SOTentFunction_Read(file_t2, 2); 

  grad_dot = SOTentFunction_GradDotBoth(tf1, tf2);

  printf(" Teste! -> %.16g \n\n", grad_dot);

  FILE *file_out = open_write("testout.fun", TRUE);
  tf1->m->write(tf1, file_out);

  return 0;
}
