/*
[]--------------------------------------------------------------------------[]
| Este programa faz a validacao das operacoes em aritmetica afim verificando |
| se o resultados obtidos atraves das aritmeticas afim e ordinaria possuem   |
| pelo menos um ponto em comum.                                              |
[]--------------------------------------------------------------------------[]
*/

#include <ioprotos.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <flt.h>
#include <ia.h>
#include <aa.h>


#define N_TERMS 4          /* numero de termos em cada forma afim */
#define N_OP 11	           /* numero de operacoes a serem validadas */
#define N_VALIDATION 10    /* nro. vezes que o resultado da op. e' validado */
#define N_FA 20

#define DEBUG 1

AATerm eps[30];   /* vetores com os indices e valores dos e_i gerados aleatoriamente */
int    neps;

void gera_epsilons (AAP aa_x, AAP aa_y);

float gera_num(void);

void main()
{
  char *op_name;
  AAP aa_x, aa_y, aa_z;
  Interval ia_x, ia_y, ia_z, ia_aaz;
  int i,j,k;
  int st;



  flt_init();
  ia_init();
  aa_init();

  op_name = (char *) malloc (15);


  for (i=0; i < N_OP; i++)
  {

    /* utiliza a hora como semente para a funcao rand  */
/*    lt = time(NULL);
    st = ((unsigned int) lt) / 2;  */
    st = 314567;
    srandom (st);
    
    for (k=0; k < N_FA; k++)
      {
      do
	aa_x = aa_throw(N_TERMS);
      while (aa_is_full(aa_x));

      do
	aa_y = aa_throw(N_TERMS);
      while (aa_is_full(aa_y));

      switch (i)
      {
	case 0 : aa_z = aa_add(aa_x,aa_y); op_name="Adicao"; break;
	case 1 : aa_z = aa_sub(aa_x,aa_y); op_name="Subtracao"; break;
	case 2 : aa_z = aa_mul(aa_x,aa_y); op_name="Multiplicacao"; break;
	case 3 : aa_z = aa_div(aa_x,aa_y); op_name="Divisao"; break;
	case 4 : aa_z = aa_max(aa_x,aa_y); op_name="Maximo"; break;
	case 5 : aa_z = aa_min(aa_x,aa_y); op_name="Minimo"; break;
	case 6 : aa_z = aa_neg(aa_x); op_name="Negacao"; break;
	case 7 : aa_z = aa_sqr(aa_x); op_name="Elevado a 2"; break;
	case 8 : ia_x=aa_range(aa_x);
		 if (ia_x.lo >= 0) aa_z=aa_sqrt(aa_x);
		 else aa_z =aa_zero();
		 op_name="Raiz Quadrada"; break;
	case 9 : aa_z = aa_inv(aa_x); op_name="Inverso"; break;
	case 10 : aa_z = aa_abs(aa_x); op_name="Modulo"; break;
      }

     for (j=0; j < N_VALIDATION; j++)
       {

	gera_epsilons (aa_x, aa_y);

	if (!aa_is_full(aa_x))
	    ia_x = aa_range(aa_fix_eps(aa_x, neps, eps));

	if (!aa_is_full(aa_y))
	    ia_y = aa_range(aa_fix_eps(aa_y, neps, eps));

	if (!aa_is_full(aa_z))
	    ia_aaz = aa_range(aa_fix_eps(aa_z, neps, eps));

	switch(i)
	{
	    case 0 : ia_z = ia_add(ia_x,ia_y); break;
	    case 1 : ia_z = ia_sub(ia_x,ia_y); break;
	    case 2 : ia_z = ia_mul(ia_x,ia_y); break;
	    case 3 : ia_z = ia_div(ia_x,ia_y); break;
	    case 4 : ia_z = ia_max(ia_x,ia_y); break;
	    case 5 : ia_z = ia_min(ia_x,ia_y); break;
	    case 6 : ia_z = ia_neg(ia_x); break;
	    case 7 : ia_z = ia_sqr(ia_x); break;
	    case 8 : if (ia_x.lo >= 0) ia_z=ia_sqrt(ia_x);  break;
	    case 9 : ia_z = ia_inv(ia_x); break;
	    case 10 : ia_z = ia_abs(ia_x); break;
	}

	fprintf(stderr,"---------------------------------------------------------- \n");
	fprintf(stderr,"Na funcao %s ### teste j = %d \n",op_name,j
	);
	fprintf(stderr,"---------------------------------------------------------- \n");

	if (DEBUG)
	  {
	    fprintf(stderr,"aa_x = "); aa_print(stderr,aa_x); fprintf(stderr,"\n");
	    if (i < 4) { fprintf(stderr,"aa_y = "); aa_print(stderr,aa_y); fprintf(stderr,"\n"); }
	    fprintf(stderr,"aa_z = "); aa_print(stderr,aa_z); fprintf(stderr,"\n");
	    fprintf(stderr,"ia_aaz = "); ia_print(stderr,ia_aaz); fprintf(stderr,"\n");
	    fprintf(stderr,"\n ---------------------------------------------------------- \n");
	    fprintf(stderr,"---------------------------------------------------------- \n \n");
	    fprintf(stderr,"ia_x = "); ia_print(stderr,ia_x); fprintf(stderr,"\n");
	    if (i < 4) { fprintf(stderr,"ia_y = "); ia_print(stderr,ia_y); fprintf(stderr,"\n"); }
	    fprintf(stderr,"ia_z = "); ia_print(stderr,ia_z); fprintf(stderr,"\n");
	  }

	if ((ia_aaz.hi < ia_z.lo) || (ia_z.hi < ia_aaz.lo))
	  {
	    fprintf(stderr, " ********* INCORRETO ********** \n");
  	    getchar();
	  }
	else
	  fprintf(stderr, " \n \n ++++++++++ Validacao OK!!! +++++++++++ \n \n");
      }
    }
  }
}

void gera_epsilons (AAP aa_x, AAP aa_y)
  {
    AATermCount neps_x, neps_y;
    AATermP xp, yp;
    int i;

    if (!aa_is_full(aa_x)) neps_x= aa_x->nterms;
    else neps_x=0;
    if (!aa_is_full(aa_y)) neps_y= aa_y->nterms;
    else neps_y=0;
    xp = (AATermP) (aa_x+1);
    yp = (AATermP) (aa_y+1);
    i=0;
    while ((neps_x > 0) || (neps_y > 0))
    {
       if (neps_x == 0) eps[i].id = yp->id;
       else if (neps_y == 0) eps[i].id = xp->id;
            else if (xp->id < yp->id) eps[i].id = xp->id;
                 else eps[i].id = yp->id;
       eps[i].coef = gera_num();
       if ((neps_x > 0) && (xp->id == eps[i].id)) {xp++; neps_x--;}
       if ((neps_y > 0) && (yp->id == eps[i].id)) {yp++; neps_y--;}
       i++;
    }
    neps = i;
  }


float gera_num(void)
{
  float n, nn, num;

  n = flt_random();
  if ((n >= 0.25) && (n < 0.5)) num=-1;
  else
    if ((n >= 0.5) && (n < 0.75)) num=0;
    else
      if (n >= 0.75) num=1;
      else
      {
         nn = (random()&65535) / 65536.0;
         if (n < 0.125) nn *= -1;
         num = nn;
       }
  return(num);
}




