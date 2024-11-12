#include <ioprotos.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ia.h>
#include <aa.h>
#include <r3.h>
#include <flt.h>
#include <irtscene.h>
#include <irtinter.h>
#include <tst_iaeval.h>
#include <tst_aaeval.h>

#define DEBUG 0

/*** INTERNAL PROTOTYES ***/

void main(int argc, char *argv[]);

Interval tst_irt_seg_eval_ia (
    pcode_proc_t *proc,   /* The function's pseudo-code */
    Interval *regs,       /* Evaluation registers */
    Interval *stack       /* Evaluation stack */
  );
  /* Avalia a funcao sc->proc utilizando aritmetica intervalar ordinaria
     retornando o intervalo resultante desta avaliacao */

void compute_aa_WXYZ (
    h3_point_t *org,      /* Start of ray */
    h3_point_t *dst,      /* End of ray */
    Interval *tv,         /* Parameter interval (in [0__1]) */
    AAP *regs             /* Evaluation registers */
  );
  /* Avalia a funcao sc->proc utilizando aritmetica afim
     retornando a forma afim resultante desta avaliacao */

int merge_aa (AAP aa1, AAP aa2, AATerm *eps);

void merge_eps (AATerm eps1[], int size_eps1, AATerm eps2[], int size_eps2);

void fix_ia_WXYZ (  
    AAP *regs_aa,
    Interval *regs_ia
  );

  
AAP tst_irt_seg_eval_aa (
    pcode_proc_t *proc,   /* The function's pseudo-code */
    AAP *regs,            /* Evaluation registers */
    AAP *stack            /* Evaluation stack */
  );

float gera_num(void);
/* gera um numero aleatorio entre [-1 __ 1] com a seguinte distribuicao
   de probabilidade: 25% de chance de sair -1; 25% de chance de sair 0;
   25% de chance de sair 1 e 25% de chance de sair um outro numero diferente
   desses
*/


AATerm eps[50];
AATermCount neps;

/*** MAIN PROGRAM ***/

void main(int argc, char *argv[])
  {
    scene_t sc_ia, sc_aa;  /* scene for interval and affine arithmetic, respc. */
    
    Interval ia_f;         /* resultado da chamada de ia_eval(sc_ia->proc) */
    AAP aa_f;              /* resultado da chamada de aa_eval(sc_aa->proc) */
    Interval ia_aaf;       /* intervalo contendo todos os valores obtidos
                              a partir de aa_f fixando-se os epsilons com
                              valores aleatorios                           */

    Interval tv;           /* intervalo onde o raio sera' avaliado;
                              e' fornecido pelo usuario.                   */

    char str[15];
    double snum;

    h3_point_t eye;        /* Homogeneous coordinates of observer */
    r3_t dir;              /* Ray's direction */
    h3_point_t inf;        /* Ray's endpoint (at infinity) */

    AAP *regs_aa;
    int i,j;
    char *s;

    AATerm eps1[30], eps2[30];
    int size_eps1, size_eps2;
    MemP frame;
    AAP aa_temp;
      
    affirm(argc == 2, "irt_main: bad parameters\nusage: irt <scene name>");

    ia_init();
    aa_init();
  
    frame = aa_top();

    s = (char *) malloc(30);
    strcpy(s,argv[1]);
    sc_aa.name = strcat(s,"-aa");
    strcpy(s,argv[1]);
    sc_ia.name = strcat(s,"-ia");
        
    irt_read_scene(sc_aa.name, &sc_aa);
    irt_read_scene(sc_ia.name, &sc_ia);
    
    h3_from_cart(&sc_aa.observer, &eye);
    h3_inf_reduce(&eye); /* To improve numerical behavior */

    dir=sc_aa.focus;
    h3_infty(&dir,&inf);
    
    fprintf(stderr,"************ Intervalo onde o raio sera' avaliado *********** \n");
    fprintf(stderr,"lo = ");
    gets(str);
    snum = atof(str);
    tv.lo = (float) snum;
    fprintf(stderr,"hi = ");
    gets(str);
    snum = atof(str);
    tv.hi = (float) snum;

    fprintf(stderr," eye = "); h3_print_point(stderr,&eye); fprintf(stderr,"\n");
    fprintf(stderr," inf = "); h3_print_point(stderr,&inf); fprintf(stderr,"\n");
    fprintf(stderr,"tv ="); ia_print(stderr,tv); fprintf(stderr,"\n");
    
    srandom(1);
    
    compute_aa_WXYZ (&eye, &inf, &tv, sc_aa.seg_regs);

    regs_aa = (AAP *) sc_aa.seg_regs;

    size_eps1 = merge_aa (regs_aa[0], regs_aa[1], eps1);
    size_eps2 = merge_aa (regs_aa[2], regs_aa[3], eps2);
    
    merge_eps (eps1, size_eps1, eps2, size_eps2);

    for (j=0; j < 100; j++)
      {
        for (i=0; i < neps; i++) eps[i].coef = gera_num();
    
/*	for (i=0; i < neps; i++)
	   fprintf(stderr, "eps[%d].id = %ld  eps[%d].coef = %f \n",i,eps[i].id,i,eps[i].coef);    
*/
        fix_ia_WXYZ (regs_aa,sc_ia.seg_regs);

	aa_f = tst_irt_seg_eval_aa (&(sc_aa.proc), sc_aa.seg_regs, sc_aa.seg_stack);

        aa_temp = aa_fix_eps(aa_f,neps,eps);
	ia_aaf = aa_range (aa_temp);

	ia_f = tst_irt_seg_eval_ia (&(sc_ia.proc), sc_ia.seg_regs, sc_ia.seg_stack);

	fprintf(stderr,"\n");
	fprintf(stderr,"aa_f = "); aa_print(stderr,aa_f); fprintf(stderr,"\n");

	fprintf(stderr,"---------------------------------------------------------- \n");

	fprintf(stderr,"ia_aaf = "); ia_print(stderr,ia_aaf); fprintf(stderr,"\n");

	fprintf(stderr,"---------------------------------------------------------- \n");

	fprintf(stderr,"ia_f = "); ia_print(stderr,ia_f); fprintf(stderr,"\n");

	fprintf(stderr,"---------------------------------------------------------- \n");

	if ((ia_aaf.hi < ia_f.lo) || (ia_f.hi <= ia_aaf.lo))
	  {

	    fprintf(stderr, " \n\n ********* INCORRETO ********** \n");
	    getchar();
	  }
	else
	  fprintf(stderr, "\n +++++++++++ Validacao OK!!! +++++++++ \n");
      }
    aa_flush(frame);
  }
		     


Interval tst_irt_seg_eval_ia (
    pcode_proc_t *proc,   /* The function's pseudo-code */
    Interval *regs,       /* Evaluation registers */
    Interval *stack       /* Evaluation stack */
  )
  { 
    tst_ia_eval (regs, stack, proc->code);
    ROUND_NEAR;
    return(stack[0]);
  }


void compute_aa_WXYZ (
    h3_point_t *org,      /* Start of ray */
    h3_point_t *dst,      /* End of ray */
    Interval *tv,         /* Parameter interval (in [0__1]) */
    AAP *regs             /* Evaluation registers */
  )
  {
    AAP aat = aa_from_interval (*tv);
    AAP aar = aa_affine (aat, -1.0, 1.0, 1.0, 0.0); /* r = 1-t */
    int i;

    for (i=0; i<4; i++)
      { regs[i] = aa_affine_2 (aar, org->c[i], aat, dst->c[i], 1.0, 0.0, 0.0); 
fprintf(stderr,"regs[%d] = ",i); aa_print(stderr,regs[i]); fprintf(stderr,"\n");      
      }

  }


int merge_aa (AAP aa1, AAP aa2, AATerm *eps_aux)
  { 
    AATermCount neps1, neps2;
    AATermP ip1, ip2;
    int i;

    if (!aa_is_full(aa1)) neps1= aa1->nterms;
    else neps1=0;
    if (!aa_is_full(aa2)) neps2= aa2->nterms;
    else neps2=0;
    ip1 = (AATermP) (aa1+1);
    ip2 = (AATermP) (aa2+1);
    i = 0;
    while ((neps1 > 0) || (neps2 > 0))
    {

       if (neps1 == 0) eps_aux[i].id = ip2->id;
       else if (neps2 == 0) eps_aux[i].id = ip1->id;
	    else if (ip1->id < ip2->id) eps_aux[i].id = ip1->id;
		 else eps_aux[i].id = ip2->id;
       eps_aux[i].coef = 0;
       if ((neps1 > 0) && (ip1->id == eps_aux[i].id)) {ip1++; neps1--;}
       if ((neps2 > 0) && (ip2->id == eps_aux[i].id)) {ip2++; neps2--;}
       i++;
     }
     return(i);
  }


void merge_eps (AATerm eps1[], int size_eps1, AATerm eps2[], int size_eps2)
  { 
    int i1, i2, i;

    i1 = 0;
    i2 = 0;
    i = 0;
    while ((i1 < size_eps1) || (i2 < size_eps2))
    {
       if (i1 == size_eps1) eps[i] = eps2[i2];
       else if (i2 == size_eps2) eps[i] = eps1[i1];
	    else if (eps1[i1].id < eps2[i2].id) eps[i] = eps1[i1];
		 else eps[i] = eps2[i2];
       eps[i].coef = 0;
       if ((i1 < size_eps1) && (eps1[i1].id == eps[i].id)) {i1++; }
       if ((i2 < size_eps2) && (eps2[i2].id == eps[i].id)) {i2++; }
       i++;
    }
    neps = i;
  }

  
void fix_ia_WXYZ (  
    AAP *regs_aa,
    Interval *regs_ia
  )
  { int i;
    
    for (i=0; i < 4; i++) regs_ia[i] = aa_range(aa_fix_eps(regs_aa[i], neps, eps));
  }
  
AAP tst_irt_seg_eval_aa (
    pcode_proc_t *proc,   /* The function's pseudo-code */
    AAP *regs,            /* Evaluation registers */
    AAP *stack            /* Evaluation stack */
  )
  { MemP frame = aa_top();
  
    tst_aa_eval (regs, stack, proc->code);
    ROUND_NEAR;
    aa_flush(frame);
    return(stack[0]);
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




