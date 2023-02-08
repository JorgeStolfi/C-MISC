/* Last edited on 2007-12-23 22:47:03 by stolfi */
/* Test of libdgrid routines. */

#include <stdio.h>
#include <math.h>
#include <SOGrid.h>
#include <SOTent.h>
#include <SOBasic.h>
#include <SOFunction.h>
#include <SOIntegral.h>

#include <dg_grid.h>

#include <ps.h>

#define DIM 3

#define NUM_AXIS 2

#define next_axis(P) ((P+1)%NUM_AXIS)

typedef struct CellGeom{   /* Cell geometry relative to root or parent cell: */
  double min[MAX_PDIM];  /* Minimum coordinate on each axis. */
  double max[MAX_PDIM];  /* Maximum coordinate on each axis. */
} CellGeom;

typedef bool_t split_criterion(CellGeom *geom);

//double genf(double *varlist, dg_dim_t dim);

double length(double dx, double dy, double dz);

bool_t toobig(CellGeom* geom);

void dg_build_tree_test(SOGrid_Tree *t);

void dg_shatter_node(
  dg_tree_node_t *p, int maxlevel,
  split_criterion split,
  CellGeom* dim,
  dg_axis_t long_axis
  );

/* IMPLEMENTATIONS */

int main(int argc, char **argv)
{ int i;
  double I_Formula, C_Formula;
  //  double I_Gauss, C_Gauss;
  SOGrid_Tree *T; 

  SOTent_Pair tp;
  // SOTent t;

  //  SOTent_vec_t tv;
  SOTent_Pair_vec_t tpv; 
  CellGeom geom;

  dg_axis_t long_axis = 0;

  geom.min[0] = 0.0;
  geom.max[0] = 1.0;
  geom.min[1] = 0.0;
  geom.max[1] = sqrt(2.0)/2;
  
  T = SOGrid_Tree_new(2);

  dg_shatter_node(T->root, 6, toobig, &geom, long_axis);
  
  //dg_build_tree_test(T);

  tpv = SOTent_minimal_pairs(T);

  tp.cx[0] = tpv.el[0].cx[0]; tp.cx[1] = tpv.el[0].cx[1];
  tp.r[0]  = tpv.el[0].r[0];  tp.r[1]  = tpv.el[0].r[1];

  //I = dg_tent_esc_prod_gauss(&tp, 2, 4, 2);

  // I_Gauss = C_Gauss = 0.0;
  // SOTent_tt_dot_gauss(&tp, 2, 6, &I_Gauss, &C_Gauss);

  I_Formula = C_Formula = 0.0;
  SOTent_tt_dot(&tp, 2, &I_Formula, &C_Formula);

  dg_cell_index_t Index0, Index1;

  printf("Rank0 = %d, Rank1 = %d \n", tp.r[0], tp.r[1]);
  for(i = 0; i < 4; i++) 
  { Index0 = SOTent_get_brick(tp.cx[0], 2, tp.r[0], i); 
    Index1 = SOTent_get_brick(tp.cx[1], 2, tp.r[1], i);

    printf("Index0 = %ld, Index1 = %ld \n", Index0, Index1); 
  }  
    
  //    printf("Tent dot product (Gauss) = %20.15f , Correction: %20.15f \n", I_Gauss, C_Gauss);
    printf("Tent dot product (Formula) = %20.15f , Correction: %20.15f \n", I_Formula, C_Formula);

  //  fprintf("Tent dot product (Gauss) = %f \n", I_Gauss);
  //  fprintf("Tent dot product (Formula) = %f \n", I_Formula);

  return(0);

}

/*
double genf(double *var, dg_dim_t dim)
{ int i; 
  double result = 1;
   
  for(i = 0; i < dim; i++) result = result * (double)(pow(var[i],3)+5);
  return(result);
}
*/

double length(double dx, double dy, double dz)
{
  return sqrt(dx*dx + dy*dy + dz*dz);
}


bool_t toobig(CellGeom* geom)
/*
  An arbitrary splitting criterion. */
{
  double xctr = (geom->min[X] + geom->max[X])/2.0;
  double yctr = (geom->min[Y] + geom->max[Y])/2.0;
  double r = length((geom->max[X] - geom->min[X])/2.0, (geom->max[Y] - geom->min[Y])/2.0, 0.0);
  double d1 = length((xctr - 0.2)/9.0, (yctr - 0.3)/9.0, 0.003);
  double d2 = length((xctr - 0.6)/5.0, (yctr - 0.5)/5.0, 0.003);
  double d3 = length((yctr - (0.9 + 0.4*sin(2.0*xctr)))/4.0, 0.0, 0.003);
  double maxr = 3.0/(1.0/d1 + 1.0/d2 + 1.0/d3);
  return (r > maxr);
}


void dg_shatter_node(
  dg_tree_node_t *p, int maxlevel,
  split_criterion split,
  CellGeom* geom,
  dg_axis_t long_axis

 )
/*
  Given a node with no children, replaces it
  by a random subtree of depth at most maxlevel. */
{
  CellGeom d;

  d.min[X] = geom->min[X];
  d.max[X] = geom->max[X];
  d.min[Y] = geom->min[Y];
  d.max[Y] = geom->max[Y];

  if ((maxlevel > 1) && (split(geom)))
    {
      dg_tree_node_split(p);
     
      if (long_axis == 0)
	{
	  double xmid = (geom->max[X] + geom->min[X])/2.0;
	  d.max[X] = xmid;
	  dg_shatter_node(p->c[0], maxlevel-1, split, &d, 1);
	  d.max[X] = geom->max[X];
	  d.min[X] = xmid;
	  dg_shatter_node(p->c[1], maxlevel-1, split, &d, 1);
	}
      else
	{
	  double ymid = (geom->max[Y] + geom->min[Y])/2.0;
	  d.max[Y] = ymid;
	  dg_shatter_node(p->c[0], maxlevel-1, split, &d, 0);
	  d.max[Y] = geom->max[Y];
	  d.min[Y] = ymid;
	  dg_shatter_node(p->c[1], maxlevel-1, split, &d, 0);
	}
     }
  
}

void dg_build_tree_test(SOGrid_Tree *t)
{ int tam, i, inc = 0;
  dg_cell_index_t new_index, aux;
  dg_tree_node_t* p;

  tam = sizeof(new_index)*8;
  printf(" * Type the index of a new node: \n");
  scanf("%ld",&new_index);

  while(new_index!=0)
  {
    aux=0;   		/* Searches for the first bit 1 */
    i=tam;
    while(aux==0)
    {
      aux=new_index;
      i--;
      aux>>=i;
    }

    new_index<<=tam-i;       /* Clears the first bit 1 */
    new_index>>=tam-i;
    i--;

    p=t->root;
    while(i >= 0)
    {
      aux=new_index;
      aux<<=tam-i-1;
      aux>>=tam-1;

      if(aux)p=p->c[1];
      else p=p->c[0];
     i--;
    }
    dg_tree_node_split(p);
   
    printf(" * Type the index of a existing node: \n");
    scanf("%ld",&new_index);
    printf(" * Type the side to append it (0 | 1): \n");
    scanf("%d",&inc);
    new_index = new_index*2 + inc;
    printf(" *** The NEW node's index -> %ld \n" , new_index);
  }
}
