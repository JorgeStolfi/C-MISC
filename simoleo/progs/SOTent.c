/* See SOTent.h */
/* Last edited on 2007-01-04 04:36:20 by stolfi */

#include <SOGrid.h>
#include <SOBasic.h>
#include <SOTent.h>

#include <dg_grid.h>
#include <nat.h>

#include <math.h>

/* INTERNAL PROTOTYPES */

void SOTent_tf_cell_dot(
  int pos, 
  dg_cell_index_t k, 
  dg_rank_t rank, 
  dg_dim_t d, 
  int npoints, 
  SOIntegral_Func f,
  dg_dim_t fd,
  double *sum,
  double *corr
);

void SOTent_tt_cell_dot(
  int *pos, 
  dg_cell_index_t *k, 
  dg_rank_t *rank,
  dg_dim_t d, 
  int drank,
  double *sum,
  double *corr
);
  
void SOTent_tt_cell_grad_dot(
  int *pos, 
  dg_cell_index_t *cell, 
  dg_rank_t *rank, 
  dg_dim_t d, 
  int drank, 
  double *sum,
  double *corr
);

double SOTent_eval_brick(dg_dim_t d, int *s, double *xb);
  /* Evaluates a {d}-dimensional tent function {f(x)} given the 
    coordinates {xb} of the argument point relative to the 
    tent brick that contains {x}.  The bit vector {s}
    indicates where the pole of the tent is, relative to the 
    brick: namely, the pole is at {xb[i] = (double)s[i]}. */

bool_t SOTent_common_cells(dg_rank_t *r, dg_tree_node_t ***C, int *n);
  /* Given two lists of tent bricks, {C[0][0..n[0]-1]} and {C[1][0..n[1]-1]},
    of ranks {r[0]} and {r[1]}, respectively, returns {true} iff the
    two lists contain a common cell. */

/* IMPLEMENTATIONS */

double SOTent_eval(dg_dim_t d, dg_rank_t r, double *x)
  { double T = 1.0;
    int i;
    for(i=0; i < d; i++) 
      { if(x[i] > 1 || x[i] < -1) return 0.0;
        double fac = (i <= r ? 1.0 : 1.0 - fabs(x[i])); 
        T *= fac; 
      }
    return T;
  }

double SOTent_eval_brick(dg_dim_t d, int *s, double *x)
  { double T = 1;
    int i;
    for(i=0; i < d; i++) { T = T * (s[i] ? 1.0 - x[i] : x[i]); }
    return T;
  }

dg_cell_index_t SOTent_get_brick(dg_cell_index_t hbr, dg_dim_t d, dg_rank_t r, int k)
  { dg_grid_pos_t dp[dg_MAX_GRID_DIM];
    int i;
    for(i = 0; i < d; i++) { dp[i] = (k & 1) - 1; k >>= 1; }
    return dg_cell_translate(d, hbr, dp);
  }

/* Tent bases */

SOTent_vec_t SOTent_basis(SOGrid_Tree *t, bool_t inferior)
  { 
    nat_t ntv = 0;
    SOTent_vec_t tv = SOTent_vec_new(100); 
    dg_dim_t d = t->d;
    int_vec_t perp_axis = int_vec_new(d);

    auto void SOTent_basis_rec (dg_rank_t r, dg_dim_t s, dg_tree_node_t **C);
      /* Given the star {C[0..2^s-1]} of a face, appends to {tv} all
        tents whose center lies on that face. */

    void SOTent_basis_rec (dg_rank_t r, dg_dim_t s, dg_tree_node_t **C)
    {
      int i, mask;
      int NC = (unsigned char)(1 << s); /* Number of cells in star */
      dg_axis_t long_axis = dg_split_axis(d, r);

      if ((! inferior) && (s == d))
        { SOTent_vec_expand(&tv, ntv);
          tv.el[ntv] = (SOTent){C[NC-1]->index, r};
          ntv++;
        }
      
      bool_t any_leaves = FALSE;
      for( i = 0; i < NC; i++ ) 
        if((LOCH(*C[i]) == NULL) || (HICH(*C[i]) == NULL)) 
          { any_leaves = TRUE; } 

      if( any_leaves )
        { if( s == d )
            { SOTent_vec_expand(&tv, ntv);
              tv.el[ntv] = (SOTent){C[NC-1]->index, r};
              ntv++;
            } 
        }
      else // There are no leaves
        {
          bool_t split_face = TRUE; // Check whether face {f} will be split
          for( i = 0; i < s; i++ ) 
            if(long_axis == perp_axis.el[i]) 
              { split_face = FALSE;
                mask = 1; 
                mask = 1 << (s - 1 - i);
              }

          if( split_face )
            { 
              dg_tree_node_t *CX[2*NC]; /* Split of star */
              
              for(i = 0; i < NC; i++) { CX[i] = LOCH(*C[i]); }
              SOTent_basis_rec (r+1, s, CX);

              for(i = 0; i < NC; i++) { CX[i] = HICH(*C[i]); }
              SOTent_basis_rec (r+1, s, CX);

              perp_axis.el[s] = long_axis; 
              for( i=0; i<NC; i++ )
                { CX[2*i] = LOCH(*C[i]);
                  CX[2*i+1] = HICH(*C[i]);
                }
              SOTent_basis_rec (r+1, s+1, CX);
            }
          else
            { dg_tree_node_t *CX[NC]; /* Split of star */
              for(i = 0; i < NC; i++) { CX[i] = C[i]->ch[1 - (i & mask)]; }
              SOTent_basis_rec (r+1, s, CX);
            }
        }
    }
    SOTent_basis_rec (0, 0, &(t->root));
    SOTent_vec_trim(&tv, ntv);
    return tv;
  }

SOTent_Pair_vec_t SOTent_basis_pairs(SOGrid_Tree *t, bool_t inferior)
{ 
  nat_t ntpv = 0;
  SOTent_Pair_vec_t tpv = SOTent_Pair_vec_new(100); 
  dg_rank_t r_ini[2];
  nat_t s_ini[2];
  dg_tree_node_t ** C_ini[2], *cells_ini[2];
  int *perp_ini[2], perp, dup[2];
  dg_dim_t d = t->d;
    
  auto void SOTent_basis_pairs_rec(
    dg_rank_t *rank, 
    nat_t *S, 
    dg_tree_node_t ***cells, 
    int **Perp, 
    int *dup
  );
  
  void SOTent_basis_pairs_rec(
    dg_rank_t *rank, 
    nat_t *S, 
    dg_tree_node_t ***cells, 
    int **Perp, 
    int *dup
  )
  {
    int i, j, mask, divide, other, duplicated[2];
    int NC[2], *Next_perp[2];
    dg_rank_t r[2];
    nat_t s[2];
    dg_tree_node_t **cellsX, **Next_els[2];
    dg_axis_t long_axis;
    bool_t any_leaves[2], parallel_split;

    int *perp_axisI, *Next_perpI[2];
    dg_tree_node_t **cellsI, **Next_elsI[2];
    dg_rank_t rI[2];
    nat_t sI[2];
    
    for(i = 0; i < 2; i++) NC[i] = (unsigned char)(1 << S[i]);      
             
    if(SOTent_common_cells(rank, cells, NC))
    { 
      if((! inferior) && (S[0] == d) && (S[1] == d))
        { SOTent_Pair_vec_expand(&tpv, ntpv);
          tpv.el[ntpv].cx[0] = cells[0][0]->index;
          tpv.el[ntpv].cx[1] = cells[1][0]->index;
          tpv.el[ntpv].r[0] = rank[0];
          tpv.el[ntpv].r[1] = rank[1];
          ntpv++;     
        } 
      
      for(i = 0; i < 2; i++)
        { any_leaves[i] = FALSE;
          for( j = 0; j < NC[i]; j++ ) 
            if((LOCH(*cells[i][j]) == NULL) || (HICH(*cells[i][j]) == NULL))
              { any_leaves[i] = TRUE; }           
        }
       
      if(any_leaves[0] && any_leaves[1])
        { 
          if((S[0] == d) && (S[1] == d))
            { SOTent_Pair_vec_expand(&tpv, ntpv);
              tpv.el[ntpv].cx[0] = cells[0][NC[i]-1]->index;
              tpv.el[ntpv].cx[1] = cells[1][NC[i]-1]->index;
              tpv.el[ntpv].r[0] = rank[0];
              tpv.el[ntpv].r[1] = rank[1];
              ntpv++;     
            } 
          return;
        }  

      divide = 0; 
      if(rank[1] < rank[0]) divide = 1;  // The element with bigger cells shall be divided

      if(any_leaves[divide]) 
        { divide = 1 - divide; } // Unless the element with bigger cells has a leaf
      other = (divide + 1) % 2;

      perp_axisI = (int *)malloc((S[divide]+1)*sizeof(int));
      cellsX = (dg_tree_node_t **)malloc(NC[divide]*sizeof(dg_tree_node_t *));
      
      r[0] = rank[0]; r[1] = rank[1];
      s[0] = S[0]; s[1] = S[1];
    
      long_axis = rank[divide] % d;
      Next_els[other] = cells[other]; Next_els[divide] = cellsX;
      r[divide]++; // All new elements shall have increased ranks
      
      parallel_split = TRUE; // check parallelism

      for( i=0; i<s[divide]; i++ ) 
      { perp_axisI[i] = Perp[divide][i];
        if( long_axis == Perp[divide][i] )
        { 
          parallel_split = FALSE;
          mask = 1; 
          mask <<= (s[divide] - 1) - i;
        }
      }

      if( parallel_split ) // It's parallel
      {         
        cellsI = (dg_tree_node_t **)malloc(2*NC[divide]*sizeof(dg_tree_node_t *));
        rI[0] = rank[divide]+1; rI[1] = rank[divide]+1;
        sI[0] = S[divide]; sI[1] = S[divide]+1; 
        Next_perp[other] = Perp[other]; Next_perp[divide] = perp_axisI;
        Next_perpI[0] = perp_axisI; Next_perpI[1] = perp_axisI;
        perp_axisI[s[divide]] = long_axis;
        Next_elsI[0] = cellsX; Next_elsI[1] = cellsI;

        for( i=0; i<NC[divide]; i++ )
        { cellsI[2*i] = LOCH(*cells[divide][i]);
          cellsI[2*i+1] = HICH(*cells[divide][i]);
        }

        for( i=0; i<NC[divide]; i++ ) cellsX[i] = LOCH(*cells[divide][i]);
        duplicated[0] = 1; duplicated[1] = 1; 
	SOTent_basis_pairs_rec(r, s, Next_els, Next_perp, duplicated);
        duplicated[0] = 0; duplicated[1] = 0;
        if(!dup[divide])SOTent_basis_pairs_rec(rI, sI, Next_elsI, Next_perpI, duplicated);
        
        for( i=0; i<NC[divide]; i++ ) cellsX[i] = HICH(*cells[divide][i]);
        duplicated[0] = 1; duplicated[1] = 1;      
        SOTent_basis_pairs_rec(r, s, Next_els, Next_perp, duplicated);
        duplicated[0] = 0; duplicated[1] = 1;
        if(!dup[divide])SOTent_basis_pairs_rec(rI, sI, Next_elsI, Next_perpI, duplicated);
         
        Next_els[divide] = cellsI;
        s[divide]++;
        duplicated[divide] = 1; duplicated[other] = dup[other];
        SOTent_basis_pairs_rec(r, s, Next_els, Next_perp, duplicated);
        
        free(cellsI);
      }
      else
      {
        for( i=0; i<NC[divide]; i++ )
        { 
          if( (i & mask) != 0 ) cellsX[i] = LOCH(*cells[divide][i]);
          else cellsX[i] = HICH(*cells[divide][i]);
        }
        SOTent_basis_pairs_rec(r, s, Next_els, Perp, dup);       
      } 
      free(cellsX);
      free(perp_axisI);
    } 
  }
  
  if ((LOCH(*(t->root)) != NULL) && (HICH(*(t->root)) != NULL))
  {
    r_ini[0] = 1; r_ini[1] = 1;
    s_ini[0] = 0; s_ini[1] = 1;
    cells_ini[0] = LOCH(*t->root); cells_ini[1] = HICH(*t->root);
    C_ini[1] = cells_ini; 
    perp = 0; perp_ini[1] = &perp; // X axis

    C_ini[0] = &(LOCH(*t->root)); dup[0] = 0; dup[1] = 0;
    SOTent_basis_pairs_rec ( r_ini, s_ini, C_ini, perp_ini, dup);

    C_ini[0] = &(HICH(*t->root)); dup[0] = 0; dup[1] = 1;
    SOTent_basis_pairs_rec ( r_ini, s_ini, C_ini, perp_ini, dup);
  }
   
  SOTent_Pair_vec_trim(&tpv, ntpv);
  return tpv;
}

bool_t SOTent_common_cells(dg_rank_t *r, dg_tree_node_t ***C, int *NC)
{
  int i, j, shorter = 0, greater, dr;

  if(r[1] < r[0]) shorter = 1; 
  greater = (shorter + 1) % 2;
  
  dr =(int)r[greater]-r[shorter];
 
  for(i = 0; i < NC[shorter]; i++)
    for(j = 0; j < NC[greater]; j++)
      if(C[shorter][i]->index == ((C[greater][j]->index)>>dr)) return TRUE;
  return FALSE;
}

void SOTent_tf_dot(
  SOTent *t, 
  dg_dim_t d, 
  SOIntegral_Func f,
  dg_dim_t fd,
  int npoints,
  double *sum,
  double *corr
)
{ 
  int pos, max_cells;
  dg_cell_index_t k;
  
  max_cells = 1 << d;

  for(pos = 0; pos < max_cells; pos++)
  { 
    k = SOTent_get_brick(t->cx, d, t->r, pos);
    SOTent_tf_cell_dot(pos, k, t->r, d, npoints, f, fd, sum, corr);
  }
}

void SOTent_tf_cell_dot(
  int pos, 
  dg_cell_index_t k,
  dg_rank_t rank, 
  dg_dim_t d, 
  int npoints, 
  SOIntegral_Func f,
  dg_dim_t fd,
  double *sum,
  double *corr
)
{   
  static interval_t box[MAX_PDIM]; /* Coords of cell relative to root cell. */

  auto void integrand(double *var, double *fvar);
    /* Integrand of dot product: computes {h(x) = f(x)*t(x)},
      where {t} is the tent function, and {x} is a point within
      the cell. Note that {x} is given in cell-relative coordinates. */

  void integrand(double *x, double *fx)
    { double tx = 1;
      int i, j;
      /* Evaluate the tent function at {var}: */
      for(i = 0; i < d; i++)
      { /* Compute the position {s} of the cell within the tent's domain,
          along axis {i}: {s = 0} = low side, {s = 1} = high side. */
        int shift = ((rank+i)%d);
        int s = (pos>>shift) & 1;
        tx *= (s ? 1 - x[i] : x[i] );
        /* Convert {x} to root-relative coordinates: */
        x[i] = LO(box[i]) + x[i] * (HI(box[i]) - LO(box[i]));
      } 

      /* Evaluate the function {f} at {var}: */
      f(x, fx);

      /* Multiply fvar by the tent value: */
      for (j = 0; j < fd; j++) { fx[j] *= tx; }
    }

  /* Compute cell coordinates relative to root cell: */
  dg_cell_box_root_relative(d, k, box);
  /* Compute integral of {f(x)*t(x)} within the cell: */
  SOIntegral_Gauss(integrand, d, fd, sum, corr);
  
}


void SOTent_tt_dot(
  SOTent_Pair *tp,   
  dg_dim_t d,         /* Dimension of domain. */
  double *sum,      /* Accumulator for dot product. */
  double *corr      /* Low-order bits of {*sum}. */
)
{
  /* cell[0]: the smaller, with higher rank; 
     cell[1]: the bigger, with lower rank; 
     Same thing for rank[0] and rank[1] */

  int pos[2], lower, higher, drank, max_cells;
  dg_cell_index_t cell[2];
  dg_rank_t rank[2];

  max_cells = 1 << d;
  lower = 0; if(tp->r[1] < tp->r[0]) lower = 1;  // Gets the lower Rank
  higher = (lower + 1) % 2;
    
  for(pos[1] = 0; pos[1] < max_cells; pos[1]++)
    for(pos[0] = 0; pos[0] < max_cells; pos[0]++)
    { 
      cell[1] = SOTent_get_brick(tp->cx[lower], d, tp->r[lower], pos[1]);
      cell[0] = SOTent_get_brick(tp->cx[higher], d, tp->r[higher], pos[0]);
      rank[1] = tp->r[lower]; rank[0] = tp->r[higher];
      drank = tp->r[higher] - tp->r[lower];

      if((cell[0] >> drank) == cell[1]) // if lower is inside higher  
        SOTent_tt_cell_dot
         (pos, cell, rank, d, drank, sum, corr);           
    }
}

void SOTent_tt_cell_dot(
  int *pos, 
  dg_cell_index_t *cell, 
  dg_rank_t *rank, 
  dg_dim_t d, 
  int drank, 
  double *sum,
  double *corr
)
{ 
  /* cell[0]: the smaller, with higher rank; 
     cell[1]: the bigger, with lower rank; 
     Same thing for rank[0] and rank[1] */
  
  int a, i, j, kk, shift_l, shift_h, rel_diff, sl, sh; 
  double integral=1, I, l_start, l_size;
  int volcorr=0, fulldivs = rank[0]/d; /* Volume correction of integral relative to root cell */ 
  
  a = d - (drank % d);

  for(i = 0; i < d; i++)
  { volcorr += fulldivs + ((rank[0]%d) > i ? 1 : 0);
    
    rel_diff = (drank/d); if((drank%d)>i) rel_diff++;
    kk = 0; 
    for(j = 0; j < rel_diff; j++) 
      kk = kk + ((cell[0] & (1 << ((j * d)+i))) >> (((j * d)+i)-j)); 
           
    l_start = (double) kk / (1 << rel_diff);
    l_size = (double) 1 / (1 << rel_diff);
            
    shift_l = ((a + i)%d); shift_h = (i%d); 
    sl = 1 - ((pos[1]&(1<<shift_l))>>shift_l); 
    sh = 1 - ((pos[0]&(1<<shift_h))>>shift_h);        
    if(sl) kk = 1; else kk = -1;

    I = kk*(((double)1/(2-sh))*((double)l_size/3)+((double)l_start/2))+((1-sl)*(double)1/2);

    integral = integral * I;
  }       
  integral = (double)integral / (1 << volcorr);

  { double term = integral;
    /* Kahan's summation formula: */
    double tcorr = term - *corr;
    double newSum = *sum + tcorr;
    *corr = (newSum - *sum) - tcorr;
    *sum = newSum;
  }
}

void SOTent_tt_grad_dot(
  SOTent_Pair *tp,   
  dg_dim_t d,     /* Dimension of domain. */
  double *sum,      /* Accumulator for dot product. */
  double *corr      /* Low-order bits of {*sum}. */
)
{
  /* cell[0]: the smaller, with higher rank; 
     cell[1]: the bigger, with lower rank; 
     Same thing for rank[0] and rank[1] */

  int pos[2], lower, higher, drank, max_cells;
  dg_cell_index_t cell[2];
  dg_rank_t rank[2];

  max_cells = 1 << d;
  lower = 0; if(tp->r[1] < tp->r[0]) lower = 1;  // Gets the lower Rank
  higher = (lower + 1) % 2;
    
  for(pos[1] = 0; pos[1] < max_cells; pos[1]++)
    for(pos[0] = 0; pos[0] < max_cells; pos[0]++)
    { 
      cell[1] = SOTent_get_brick(tp->cx[lower], d, tp->r[lower], pos[1]);
      cell[0] = SOTent_get_brick(tp->cx[higher], d, tp->r[higher], pos[0]);
      rank[1] = tp->r[lower]; rank[0] = tp->r[higher];
      drank = tp->r[higher] - tp->r[lower];

      if((cell[0] >> drank) == cell[1]) // if lower is inside higher  
        SOTent_tt_cell_grad_dot
         (pos, cell, rank, d, drank, sum, corr);           
    }
}

void SOTent_tt_cell_grad_dot(
  int *pos, 
  dg_cell_index_t *cell, 
  dg_rank_t *rank, 
  dg_dim_t d, 
  int drank, 
  double *sum,
  double *corr
)
{ 
  /* cell[0]: the smaller, with higher rank; 
     cell[1]: the bigger, with lower rank; 
     Same thing for rank[0] and rank[1] */
  
  int a, i, j, kk, shift_l, shift_h, rel_diff, s[d][2]; 
  double l_start, l_size[d], I[d];
  int volcorr=0, fulldivs = rank[0]/d; /* Volume correction of integral relative to root cell */  


  a = d - (drank % d);

  for(i = 0; i < d; i++)
  { volcorr += fulldivs + ((rank[0]%d) > i ? 1 : 0);
    
    rel_diff = (drank/d); if((drank%d)>i) rel_diff++;
    kk = 0; 
    for(j = 0; j < rel_diff; j++) 
      kk = kk + ((cell[0] & (1 << ((j * d)+i))) >> (((j * d)+i)-j)); 
           
    l_start = (double) kk / (1 << rel_diff);
    l_size[i] = (double) 1 / (1 << rel_diff);
            
    shift_l = ((a + i)%d); shift_h = (i%d); 
    s[i][1] = 1 - ((pos[1]&(1<<shift_l))>>shift_l); // sl
    s[i][0] = 1 - ((pos[0]&(1<<shift_h))>>shift_h); // sh       
    if(s[i][1]) kk = 1; else kk = -1;

    I[i] = kk*(((double)1/(2-s[i][0]))*((double)l_size[i]/3)+
           ((double)l_start/2))+((1-s[i][1])*(double)1/2);

  }       

  for(i = 0; i < d; i++)
  {
    double integral=1;
    int sl, sh;  

    for(j = 0; j < d; j++)
      if(j != i) integral *= I[j];

    if(s[i][1]) sl = 1; else sl = -1;
    if(s[i][0]) sh = 1; else sh = -1;

    integral *= (double)(sl * sh * l_size[i]);

    integral = (double)integral / (1 << volcorr);

    { double term = integral;
      /* Kahan's summation formula: */
      double tcorr = term - *corr;
      double newSum = *sum + tcorr;
      *corr = (newSum - *sum) - tcorr;
      *sum = newSum;
    }
  }
}

/* Manipulation of vectors of tents */

vec_typeimpl(SOTent_vec_t,SOTent_vec,SOTent);

/* Manipulation of vectors of tent pairs */

vec_typeimpl(SOTent_Pair_vec_t,SOTent_Pair_vec,SOTent_Pair);


