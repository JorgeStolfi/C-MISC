/* SOTent.h -- polynomial splines defined on irregular dyadic grids */
/* Last edited on 2007-01-04 00:21:00 by stolfi */

#ifndef SOTent_H
#define SOTent_H

#include <SOGrid.h>
#include <SOIntegral.h>

#include <dg_tree.h>
#include <dg_grid.h>

#include <vec.h>

/* POLYNOMIAL SPLINES */

/* Let {G} be an arbitrary partition of some region {D} of {R^d} into
  pairwise disjoint cells. A `(polynomial) spline' for {G} is a
  function {f} from {R^d} to {R} such that {f(x)} is zero outside {D}
  and, for each cell {c} of {G}, the restriction {(f|c)(x)} is a
  polynomial on the {d} coordinates {x[0..d-1]} of the argument.
  
  Each partial function {f|c} is called a `piece' of the spline. The
  `degree' of the spline is the maximum exponent of any coordinate in
  any its pieces.  The set of all splines of degree {g} on a given 
  grid {G} constitute a linear vector space, {S^g(G)}.

  The procedures in this section are concerned with the sub-space
  {S^1_0(G)} of {S^1(G)}, consisting of the *continuous* splines of degree 1 
  for some finite *dyadic* grid {G}.  Note that since the union of all cells
  of {G} is the root cell {C_0}, any such spline is by definition zero outside 
  {G}, and (by continuity) also on the boundary of {G}.

  A polynomial of degree 1 is a multi-linear function - namely, a
  function which is affine (degree 1) along any axis. Note that if
  such a function is zero at two points of an axis-aligned line
  segment, it is zero on the whole segment. Thus, in particular, if a
  spline {f} from {S^1_0(G)} is zero at the two endpoints of some edge
  {e} of the grid, it is zero along the entire edge {e}. If {f} is
  zero on two opposite edges of some face {H} of the grid, it is zero
  on the whole face {H}. And so on. */ 

/* TENT FUNCTIONS */

/* A `tent function' is a spline {t} in {S^1_0(Q)}, where {Q}
  is a grid consisting of {2^d} axis-aligned boxes {b[0..2^d-1]}, of the same
  shape and size, that share a common corner {p}. By definition,
  within each of these boxes (the `bricks' of {t}), the tent function
  is a single non-zero multilinear function, which decreases from some
  central value {t(p)} to zero on the boundary of the tent's support {B} 
  (the union of its bricks).
  
  A `unit tent function' is a tent function that has value 1 at the
  central point. Observe that all the tent functions with the same
  set of bricks {b[i]} are multiples of the corresponding 
  unit tent function.
  
  A `dyadic tent' is a tent funtion whose bricks are cells of
  the dyadic cell tree. */
  
typedef struct SOTent  /* A unit dyadic tent function. */ 
  { dg_cell_index_t cx;
    dg_rank_t r;
  } SOTent;
  /* A diyadic tent function is completely determined by
    the index {cx} of the highest-numbred brick.  The rank {r} 
    is precomputed for programmer's convenience. */
    
dg_cell_index_t SOTent_get_brick(dg_cell_index_t hbr, dg_dim_t d, dg_rank_t r, int k);
  /* Given the grid index {hbr} of the highest brick in a
    {d}-dimensional tent domain, and an integer {k} in {0..2^d-1},
    returns the grid index of the {k}th brick of that tent. Namely, if
    bit {i} of {k} is 0 (resp. 1), then the requested brick lies in
    the lower (resp. upper) half of the tent, along coordinate axis
    {r+i}. */

/* TENT FUNCTION EVALUATION */

/* The evaluation of a unit tent function {t} at a point {p} of {R^d}
  is best divided in three steps. First, we compute the
  `support-relative' coordinates {x} of {p}, by the formula {x[i] =
  (p[r+i]-c[r+i])/h[r+i]}, where {r} is the rank of the bricks, {c} is
  the center of the tent's support, and {h[j]} is the extent along
  axis {j} of the bricks of {t}. (All indices are taken module {d}.)
  Note that the support-relative coordinate {x[0]} is the relative
  position of {p} along the longest axis of the bricks. 
  
  Then, if any of the {x[i]} lies outside the range {[-1 _ +1]}, the
  value of {t(p)} is zero.  Otherwise, the value is simply the product
  of {1 - |x[i]|}, for all {i}.
  
  In practice, the grid may be deformed by a complicated shape
  function {f}, which needs to be inverted before computing {t(p)}. 
  So, part of the computation of {x} is probaly best merged 
  with the computation of {f^{-1}(p)}. */

double SOTent_eval(dg_dim_t d, dg_rank_t r, double *x);
  /* Evaluates the {d}-dimensional tent function of rank {r}
    described by {t}, at the point whose support-relative 
    coordinates are {x[0..d-1]} and stores the result in {fp[0]}. */
 
/* ENUMERATING DYADIC TENTS IN A FINITE DYADIC GRID */

/* A dyadic tent {t} is said to be `compatible' with a grid {G} if each
   of its bricks is the union of one or more whole cells of {G}. */
   
vec_typedef(SOTent_vec_t,SOTent_vec,SOTent);
  /* A list of unit tents ({SOTent}s).  
    
    The procedure {SOTent_vec_new(n)} allocates a vector with space
    for {n} elements.
    
    The procedure {SOTent_vec_expand(&nv,ix)} makes sure that the
    vector {nv} has element {nv.e[ix]}. Reallocates {nv.e} if
    necessary, to have at least {index+1} elements -- but possibly
    more.
  
    The procedure {SOTent_vec_trim(nv,ne)} makes sure that the vector
    {nv} has exactly {ne} elements. Reallocates {nv.e} if
    necessary. */

SOTent_vec_t SOTent_basis(SOGrid_Tree *t, bool_t inferior);
  /* Returns the tent basis for the space {S^1_0(G)}, where {G} is the
    finite dyadic cell grid described by the tree {t}. Returns the
    inferior basis, consisting of all minimal tents, if {inferior=TRUE};
    else returns the superior basis, consisting of all maximal tents.
    
    A compatible tent {t} is `minimal' (resp.'maximal') if there is no
    other dyadic tent function with same center whose support is
    properly contained in (resp. properly contains) that of {t}.

    If {G} is a finite dyadic grid, it is easy to see that a dyadic tent
    {t} is minimal if, and only if, all its bricks are cells (leaf or
    non-leaf) of {G}, of the same rank {r}, and at least one of them is
    a leaf cell. It then follows that the support {B} of {t} is
    congruent to a dyadic cell of rank {r-d}, translated so that it is
    centered at the common vertex {p} of the bricks. (Note that {B}
    itself may not be a cell of the dyadic cell tree; it may
    straddle two or more such cells.)

    It is conjectured that, for any dyadic grid {G}, the set of all the
    unit minimal tent functions of {G} constitute a basis for the space
    {S^1_0(G)}. */

/* ENUMERATING PAIRS OF TENTS THAT OVERLAP */

typedef struct SOTent_Pair  /* For pairs of tents which have cells in comon. */ 
  { dg_cell_index_t cx[2];
    dg_rank_t r[2];
  } SOTent_Pair;
 
vec_typedef(SOTent_Pair_vec_t,SOTent_Pair_vec,SOTent_Pair);
  /* A list of minimal tent pairs. ({SOTent_Pair}s).  
    
    The procedure {SOTent_Pair_vec_new(n)} allocates a vector with
    space for {n} elements.
    
    The procedure {SOTent_Pair_vec_expand(&nv,ix)} makes sure that the
    vector {nv} has element {nv.e[ix]}. Reallocates {nv.e} if
    necessary, to have at least {index+1} elements -- but possibly
    more.
  
    The procedure {SOTent_Pair_vec_trim(nv,ne)} makes sure that the
    vector {nv} has exactly {ne} elements. Reallocates {nv.e} if
    necessary. */ 

SOTent_Pair_vec_t SOTent_basis_pairs(SOGrid_Tree *t, bool_t inferior);
  /* Returns a list of pairs of unit tent functions in the
    space {S^1_0(G)} which have overlapping domains; where {G} is the
    finite dyadic cell grid described by the tree {t}.  The tents 
    will betaken from the inferior basis if {inferior=TRUE}, else
    from the superior one. */

/* DOT PRODUCTS OF TENT FUNCTIONS */

/* These procedures add the computed integrals to the arguments {sum} and {corr}
  (which must be initialized by the client), using Kahan's summation formula.
  The final result will be in {sum}. */

void SOTent_tt_dot(
  SOTent_Pair *tp,   
  dg_dim_t d,     /* Dimension of domain. */
  double *sum,      /* Accumulator for dot product. */
  double *corr      /* Low-order bits of {*sum}. */
);
  /* Computes the dot product of two tent functions of dimension {d},
    using the exact algebraic formula.  The result is 
    accumulated in {*sum} and {*corr}. */

void SOTent_tt_grad_dot(
  SOTent_Pair *tp,   
  dg_dim_t d,     /* Dimension of domain. */
  double *sum,      /* Accumulator for dot product. */
  double *corr      /* Low-order bits of {*sum}. */
);

void SOTent_tf_dot( 
  SOTent *t, 
  dg_dim_t d,       /* Dimension of domain. */
  SOIntegral_Func f,  /* The function. */
  dg_dim_t fd,      /* Dimension of function's result {f(p)}. */ 
  int npoints,        /* Number of integration knots in each axis. */
  double *sum,        /* {sum[j]} accumulates component {j} of dot product. */
  double *corr        /* {corr[j]} low-order correction for {sum[j]}. */
);
  /* Computes the dot product of a tent function {t}, with domain
    dimension {d}, and a generic function {f} from {R^d} to {R^fd},
    using Gaussian quadrature with {npoints} samples along each axis.
    The result is accumulated in {sum[0..fd-1]} and
    {corr[0..fd-1]}. */

#endif
