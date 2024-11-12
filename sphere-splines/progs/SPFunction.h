/* SPFunction.h -- General functions on the sphere */
/* Last edited on 2006-02-28 11:50:36 by stolfi */

#ifndef SPFunction_H
#define SPFunction_H

#include <SPBasic.h>
#include <SPFuncMap.h>
#include <SPTriang.h>
#include <SPMatrix.h>
#include <SPH3.h>

#include <r3.h>
#include <vec.h>
#include <bool.h>
#include <nat.h>

#include <stdio.h>
 
/* An {SPFunction} object represents a real function defined on the sphere. */
  
#define OBJ void /* Actually {SPFunction} or a subclass */

typedef double SPFunction_EvalMth(OBJ *f, R3Point *p);     
typedef R3Gradient SPFunction_GradMth(OBJ *f, R3Point *p); 
typedef R3Hessian SPFunction_HessMth(OBJ *f, R3Point *p);  
typedef void SPFunction_WriteMth(OBJ *f, FILE *wr);        
typedef void SPFunction_MapleMth(OBJ *f, FILE *wr);        
typedef void SPFunction_ScaleMth(OBJ *f, double a);        
typedef void SPFunction_AddMth(OBJ *f, double a, OBJ *h);  
typedef OBJ *SPFunction_CopyMth(OBJ *f);                   
typedef void SPFunction_FreeMth(OBJ *f);        

typedef struct SPFunction_Methods    /* Methods for any SPFunction: */
  { SPFunction_EvalMth *eval;    /* Evaluates the function at {p}. */
    SPFunction_GradMth *grad;    /* Cartesian gradient at {p}. */
    SPFunction_HessMth *hess;    /* Cartesian Hessian at {p}. */
    SPFunction_WriteMth *write;  /* Writes the function to {wr}, for reading back. */
    SPFunction_MapleMth *maple;  /* Writes {f} to {wr} in Maple format. */
    SPFunction_ScaleMth *scale;  /* Scales the function {f} by {a}. */
    SPFunction_AddMth *add;      /* Adds {a*h} to {f} (same types only). */
    SPFunction_CopyMth *copy;    /* Make a clone of {f} (same meths, new data). */
    SPFunction_FreeMth *free;    /* Reclaims the function's private storage. */
  } SPFunction_Methods;
  /* Some of these methods may be missing for certain subclasses. */
  
typedef struct SPFunction_Data    /* Data fields for any SPFunction: */
  { } SPFunction_Data;            
  
typedef struct SPFunction /* Abstract class: a function defined on the sphere. */
  { char *type;             /* Type identifier */
    SPFunction_Methods *m;  /* Function methods. */
    SPFunction_Data *d;     /* Function parameters. */
  } SPFunction;
  
#define SPFunction_TypeId "SF."

SPFunction *SPFunction_Cast(OBJ *f);
  /* If {f} (according to {f->type}) is a subclass of {SPFunction},
    returns {f} cast to that type; otherwise returns NULL. */

S2Gradient SPFunction_SGrd(SPFunction *f, S2Point *p);  
  /* Computes the spherical gradient of {f} at {p},
    namely the component of {f->m->grad(f,p)} which is tangent
    to the sphere at {p}. Assumes {p} has unit norm. */

double SPFunction_SLap(SPFunction *f, S2Point *p);  
  /* Computes the spherical Laplacian of {f} at {p}. Assumes
    {p} has unit norm. */

SPH3_Plane SPFunction_GetSupportingPlane(SPFunction *f);
  /* Returns the supporting plane of {f} if it can find it;
    or {Omega} if {f} has no supporting plane. */
  
void SPFunction_M_Write(SPFunction *f, FILE *wr); 
  /* The default {write} method for a generic spherical function.
    Writes a standard {SPFunction} header line and the specific
    subclass ID, then calls the subclass-specific {write} method, and
    closes off with the generic {SPFunction} footer. */

SPFunction *SPFunction_Read(FILE *rd);  
  /* Reads a function's description from {rd}. */

/* FUNCTION SAMPLING */

S2Point_vec_t SPFunction_SamplePoints(int N);
  /* Generates N points, roughly uniformly distributed on the sphere. */

double_vec_t SPFunction_Sample(SPFunction *f, S2Point_vec_t p);
  /* Evaluates {f} at each point of {p}, returns the results. */

/* BASES FOR FUNCTION SPACES */

typedef SPFunction *SPFunctionRef;

vec_typedef(SPFunctionRef_vec_t,SPFunctionRef_vec,SPFunctionRef);

#define SPFunction_Basis SPFunctionRef_vec_t
#define Basis SPFunctionRef_vec_t 
  /* A {Basis} is a list of spherical functions,
    hopefully linearly independent. */

void SPFunction_WriteBasis(FILE *wr, Basis F);
  /* Writes a basis {F} to {wr}, in a format that can be 
    read back by {Read}.  */

Basis SPFunction_ReadBasis(FILE *rd);
  /* Reads from {rd} a basis for a space of
    spherical functions. Assumes the contents of
    {rd} was created by {write} above. */

Basis SPFunction_GetSubBasis(Basis F, int_vec_t ix);
  /* Returns a new {Basis} vector {G} of size {ix.ne} where
    {G[k] = F[ix[k]]}.  Note that the elements are shared, 
    not cloned. */

Basis SPFunction_CopyBasis(Basis F);
  /* Returns a new {Basis} vector {G} of same size as {F}, where
    each {G[i]} is a cloned copy of the corresponding {F[i]}. */

SPFunction *SPFunction_LinComb(double *a, Basis F);
  /* Returns the linear combination {SUM{a[i]*F[i], i=0..F.ne-1}}.
  
    The result is NULL if {F.ne == 0}; otherwise, it is a clone of
    the first basis element {F[0]}, modified with {add} and {scale}.
    All other functions {F[i]} with non-zero weights must be
    {add}-compatible with {F[0]}. */

/* INTEGRALS OF SPHERICAL FUNCTIONS */

/* The procedures below compute weighted integrals of a spherical
  function {F(p)}, defined as integrals of {w(p)*F(p)} over the
  sphere, where {w} is a given /weight/ function. All procedures
  assume {w(p) == 1.0} for all {p}, when {w == NULL}.
  
  A procedure {SPFunction_RawXXX} computes the same integral as
  {SPFunction_XXX}, but always uses {SPIntegral_OnSphere}, with the
  given {tri} argument, ignoring any underlying triangulation of {f}.
  
  A procedure {SPFunction_CustomXXX} computes the same integral as
  {SPFunction_XXX}, but uses {SPIntegral_Custom}, with sample points
  {sp} and sample weights {wt}. Namely, it computes
  {(A/N)*SUM(w(sp[k])*F(sp[k]), k=0..N-1)}, where {A} is the area of
  the sphere, and {N == sp.ne} 
  
  In most cases, the spherical function {F} is a given {SPFunction}
  {f} modified by a given pointwise operator {FMap} (a {FuncMap});
  i.e. {F(p) = FMap(f(p),p)}. The procedures assume {FMap(u,p) == u}
  when {FMap == NULL}. Procedures that require a second
  spherical-function operand {F}, such as {SPFunction_Dot}, accept a
  second {SPFunction} {g} and a corresponding point operator {GMap},
  with the same meaning and default behavior as {FMap}. However, some
  procedures -- especially {Raw} and {Custom} versions -- do not take
  {FMap}/{GMap} arguments, and always use {F(p) = f(p)}, {G(p) =
  g(p)}. */

double SPFunction_Integral
  ( SPFunction *f,
    FuncMap FMap, 
    SPFunction *w,
    Triangulation *tri,
    bool_t verbose
  );
  /* Computes the integral of {w(p)*FMap(f(p),p)} over the whole sphere.
    
    If the function {f} happens to be an instance of {SPSpline},
    and {tri} is either NULL or {f}'s underlying triangulation, then
    the procedure reduces to {SPSpline_IntegralPW} (q.v.).
    Otherwise, the integral is computed with {SPIntegral_OnSphere},
    with the same {tri} argument. */

double SPFunction_RawIntegral
  ( SPFunction *f,
    FuncMap FMap, 
    SPFunction *w,
    Triangulation *tri,
    bool_t verbose
  );
  /* {Raw} version of {SPFunction_Integral}. */

double SPFunction_CustomIntegral
  ( SPFunction *f, 
    FuncMap FMap, 
    SPFunction *w,
    S2Point_vec_t sp,
    double_vec_t wp
  );
  /* {Custom} version of {SPFunction_Integral}. */
    
/* CENTROIDS */ 

R3Point SPFunction_Centroid
  ( SPFunction *f,
    FuncMap FMap, 
    SPFunction *w,
    Triangulation *tri
  );
  /* The approximate weighted centroid of {F(p) = FMap(f(p),p)},
    namely {b} such that {b.c[i] == Int(p -> p.c[i]*F(p))}, for
    {i=0,1,2}, where {Int(f) = SPFunction_Integral(f,FMap,w,tri)}.
    Note that the result is not an S2Point, and is not a true average
    unless {Int(f) == 1}. */

R3Point SPFunction_RawCentroid
  ( SPFunction *f,
    FuncMap FMap, 
    SPFunction *w,
    Triangulation *tri,
    bool_t verbose
  );
  /* {Raw} version of {SPFunction_Centroid}. */ 
  
/* INNER PRODUCT OF SPHERICAL FUNCTIONS */

/* These procedures compute weighted dot product of spherical functions,
  defined as integrals of {w(p)*H(F(p),G(p))} over the sphere,
  where {w} is a weight function, and {H()} is some pointwise operator.
  The conventions are similar to those of the corresponding integration
  procedures above.  In particular, all these procedures 
  assume {w(p) == 1.0} for all {p}, when {w == NULL}. */

double SPFunction_Dot
  ( SPFunction *f,
    FuncMap FMap, 
    SPFunction *g,
    FuncMap GMap,
    SPFunction *w,
    Triangulation *tri,
    bool_t verbose
  );
  /* Computes the weighted dot product of {F(p)=FMap(f(p),p)} and
    {G(p)=GMap(g(p),p)}, that is the integral of {w(p)*F(p)*G(p)} over
    the whole sphere.
    
    If both functions {f,g} happen to be instances of {SPSpline}
    with the same triangulation, and {tri} is either NULL or that same
    triangulation, then the procedure reduces to
    {SPSpline_DotBothPW} (q.v.).
    
    Else, if either {f} or {g} happen to be an instance of {SPSpline}
    and {tri} is either NULL or that element's underlying triangulation,
    then the procedure reduces to {SPSpline_DotSinglePW} (q.v.).
    
    In all other cases, the integral is computed with {SPIntegral_OnSphere},
    with same {tri} argument. */

double SPFunction_RawDot
  ( SPFunction *f,
    FuncMap FMap, 
    SPFunction *g,
    FuncMap GMap,
    SPFunction *w,
    Triangulation *tri,
    bool_t verbose
  );
  /* {Raw} version of {SPFunction_Dot}. */
    
double SPFunction_CustomDot
  ( SPFunction *f, 
    FuncMap FMap, 
    SPFunction *g,
    FuncMap GMap,
    SPFunction *w,
    S2Point_vec_t sp,
    double_vec_t wp
  );
  /* {Custom} version of {SPFunction_Dot}. */

/* VELOCITY-GRADIENT INNER PRODUCTS */

/* These procedures integrate over the sphere the formula {w(p) *
  dot(v(p),SGrd(f)(p)) * g(p)}, where {SGrd} is the spherical gradient
  operator, and {v} is the velocity vector field corresponding to
  rigid rotation around the Z axis with unit angular speed -- that is,
  {v(x,y,z) = (-y,x,0)}. */

double SPFunction_VelSGrdDot
  ( SPFunction *f,
    SPFunction *g, 
    SPFunction *w,
    Triangulation *tri,
    bool_t verbose
  );
  /* Computes the dot product of the spherical gradients 
    {v ¤ SGrd(f)} and {g}, that is the integral of
    {w(p)*r3_dot(v(p) ¤ SGrd(f)(p))*g(p)} over the whole sphere.
    The integral is approximated as described for {SPFunction_Dot}. */

double SPFunction_RawVelSGrdDot
  ( SPFunction *f,
    SPFunction *g, 
    SPFunction *w,
    Triangulation *tri,
    bool_t verbose
  );
  /* {Raw} version of {SPFunction_VelSGrdDot}. */

double SPFunction_CustomVelSGrdDot
  ( SPFunction *f, 
    SPFunction *g,
    SPFunction *w,
    S2Point_vec_t sp,
    double_vec_t wp
  );
  /* {Custom} version of {SPFunction_VelSGrdDot}. */

/* SPHERICAL LAPLACIAN INNER PRODUCTS */

double SPFunction_SLapDot
  ( SPFunction *f,
    SPFunction *g, 
    SPFunction *w,
    Triangulation *tri,
    bool_t verbose
  );
  /* Computes the weighted scalar product of the fields {SLap(f)} and
    {g}, which, by Green's theorem, is equal to the spherical integral
    of {w(p)*r3_dot(SGrd(f)(p),SGrd(g)(p))}, where {SGrd} is the
    spherical gradient operator. The integral is approximated as
    described for {SPFunction_Dot}. */

double SPFunction_RawSLapDot
  ( SPFunction *f,
    SPFunction *g, 
    SPFunction *w,
    Triangulation *tri,
    bool_t verbose
  );
  /* {Raw} version of {SPFunction_SLapDot}. */

double SPFunction_CustomSLapDot
  ( SPFunction *f, 
    SPFunction *g,
    SPFunction *w,
    S2Point_vec_t sp,
    double_vec_t wp
  );
  /* {Custom} function of {SPFunction_SLapDot}. */

/* ORTHONORMALIZATION OF SPHERICAL FUNCTIONS - GENERAL METRICS  */

typedef double Metric(SPFunction *f, SPFunction *g);
  /* A general inner product (a symmetric bilinear functional, 
    usually positive definite). */

void SPFunction_GenMakePositive(SPFunction *f, Metric dot);
  /* Multiplies {f} by -1, if necessary, so that {dot(f,unit)}
    is non-negative --- where {unit} is the function that has
    value 1 over the whole sphere. */

double SPFunction_GenNormalize(SPFunction *f, Metric dot);
  /* Scales {f} by a positive factor so that {dot(f,f)} becomes 
    1, 0, or -1.  Returns the presumed new value of {dot(f,f)}. */

void SPFunction_GenOrthize
  ( SPFunction *f, 
    SPFunction *g, 
    Metric dot,
    double gg
  );
  /* Makes {f} orthogonal to {g}, under the metric {dot}. 
    
    More precisely, subtracts from {f} a multiple {alpha*g} of {g},
    such that {dot(f,g) == 0}. Assumes that {dot} is bi-linear, so
    {alpha = dot(f,g)/dot(g,g)}.
    
    The client must provide the value {gg = dot(g,q)}. The procedure
    fails if {dot(g,g) == 0} but {dot(f,g) != 0}. */

SPFunction *SPFunction_GenProject
  ( SPFunction *f, 
    SPFunction *g, 
    Metric dot,
    double gg
  );
  /* Returns the component of {f} that is parallel to {g}.
    
    More precisely, returns a copy {h} of {g}, scaled by a factor
    {alpha} such that {dot(f,g) == dot(h,g)} --- i.e., {alpha =
    dot(f,g)/dot(g,g)}. Fails if {dot(g,g) == 0} but {dot(f,g) != 0}.
    
    The client must provide the value {gg = dot(g,q)}. */

void SPFunction_GenOrthizeBasisAgainstBasis
  ( Basis f, 
    Basis g, 
    Metric dot,
    bool_t verbose,
    SPMatrix *B,
    SPMatrix *L
  );
  /* Modifies each {f[i]}, by linear combination with elements 
    of {g}, so that it becomes orthogonal to the space generated by {g}. 
    
    Namely, subtracts from {f[i]} some linear combination 
    {SUM { alpha[i,j]*g[j] : j}} of the {g[j]}, where the 
    coefficients {alpha[i,j]} are such that {dot(f[i], g[j])}
    becomes zero for all pairs {i,j}.
    
    The functions {f[i]} must be {add}-compatible with all {g[j]}, and
    all functions must be pairwise data-disjoint. May fail if {dot} is
    not positive definite in the space generated by the {g[j]}.
    
    The argument {B}, if not null, should point to the matrix of dot
    products for the basis {g}, namely {B[i,j] = dot(g[i],g[j]); and
    {L}, if not null, should be its lower triangular Cholesky factor.
    If {L} is null, the procedure computes it from {B}; if {B} is
    NULL, the procedure computes it from {g} and {dot}.) */

Basis SPFunction_GenProjectBasisOnBasis
  ( Basis f, 
    Basis g, 
    Metric dot,
    bool_t verbose,
    SPMatrix *B,
    SPMatrix *L
  );
  /* Computes the projection {h[i]} of each {f[i]} on the space generated
    by {g}. 
    
    Namely, sets {h[i] = SUM { alpha[i,j]*g[j] : j}}, where the coefficients
    {alpha[i,j]} are such that {dot(f[i], g[j]) == dot(h[i],g[j])} for all {i,j}. 
    
    The storage areas {fp->e} and {fo->e} must be allocated by the
    caller. If {fp == NULL}, the {fp} basis is not generated;
    otherwise the elements of {fp} will be newly allocated
    functions.  The same applies to {fo}. 
    May fail if {dot} is not positive definite in the space generated by
    the {g[j]}. The elements of {g} must be {add}-compatible with each other;
    if {fo} is requested, the {f[i]} must be {add}-compatible, too.
    
    The arguments {B} and {L} have the same meaning as in 
    {SPFunction_GenOrthizeBasisAgainstBasis}. */

void SPFunction_GenOrthonizeBasis(Basis bas, Metric dot, bool_t verbose);
  /* Makes the elements of {bas} orthogonal to each other, according
    to the scalar product {dot}, by the plain Gram-Schmidt method. May
    fail if {dot} is not positive definite. 
    
    Also normalizes all elements {b=bas[i]}, so that
    {dot(b,b)} is either 1, 0, or -1. */

void SPFunction_GenEigenFuncBasis
  ( Basis bas, 
    double dot(SPFunction *f, SPFunction *g),
    double mdot(SPFunction *f, SPFunction *g),
    double_vec_t ev,
    bool_t verbose
  );
  /* Replaces the elements of {bas} by a set of functions of their
    generated space which are orthonormal under the dot product
    {dot}, and are stationary points of the the metric induced by the
    alternate dot product {mdot}.
    
    The stationary functions are found by eigenvalue analysis of the
    matrices {B[i,j] = dot(bas[i]|bas[j])} and {G[i,j] =
    mdot(bas[i],bas[j])}. The corresponding eigenvalues are returned
    in the vector {ev}, which must have the same length as {bas}. */

int SPFunction_GenCheckOrthonization
  ( Basis bas, 
    int ini, int lim, 
    Metric dot, 
    int maxTests, 
    double maxErr, 
    int maxBad,
    bool_t verbose
  );
  /* Checks whether {dot(bas[i],bas[j]) = (i == j)}, for all {i,j} 
    in {ini..lim-1} with {i<=j}. Returns the number of pairs 
    that failed the test.
    
    Actually checks only {maxTests} pairs, chosen by some internal
    criterion. Complains if any pair has absolute error greater than
    {maxErr}. Gives up once it finds more than {maxBad} such pairs. */

/* ORTHONORMALIZATION OF SPHERICAL FUNCTIONS - SPECIFIC METRICS */

/* There follow versions of the general normalization procedures
  above, specialized for the standard integral-based inner product
  {<f|g> = SPFunction_Dot(f,NULL,f,NULL,w,tri,FALSE)}. */

void SPFunction_MakePositive
  ( SPFunction *f, 
    SPFunction *w, 
    Triangulation *tri
  );
  /* Multiplies {f} by -1, if necessary, so that its integral over the
    entire sphere (weighted by {w}) is non-negative. The integral is
    computed by {SPFunction_Integral(f,NULL,w,tri,FALSE)}. */

double SPFunction_StdNormalize
  ( SPFunction *f, 
    SPFunction *w, 
    Triangulation *tri
  );
  /* Scales {f} by a positive factor so that {<f|f>} 
    becomes 1 or 0. Returns the old value of {<f|f>}. */

void SPFunction_StdOrthize
  ( SPFunction *f, 
    SPFunction *g, 
    SPFunction *w, 
    Triangulation *tri,
    bool_t gNormal
  );
  /* Replaces {f} by {f + alpha*g} where {alpha} is such that {f}
    becomes orthogonal to {g} --- namely, {alpha = -<f,g>/<g,g>},
    where , where {<f,g> =
    SPFunction_Dot(f,NULL,g,NULL,w,tri,FALSE)}. If {gNormal} is
    true, assumes that {g} has unit norm, so that {alpha} is simply {-
    <f,g>} */

void SPFunction_StdOrthonizeBasis
  ( Basis bas, 
    SPFunction *w, 
    Triangulation *tri
  );
  /* Makes the basis elements {bas} orthogonal to each other (by the
    plain Gram-Schmidt method), under the scalar product 
    {<f,g> = SPFunction_Dot(f,NULL,g,NULL,w,tri,FALSE)}. */

void SPFunction_SLapDotEigenFuncBasis
  ( Basis bas, 
    double_vec_t ev,
    SPFunction *w, 
    Triangulation *tri,
    bool_t verbose
  );
  /* Calls {SPFunction_GenEigenFuncBasis} with the dot products
    {mdot(f,g) = <SGrd(f)|SGrd(g)>} and {dot(f,g) = <f|g>},
    where {SGrd} is the spherical gradient operator,
    and {<f,g> = SPFunction_Dot(f,NULL,g,NULL,w,tri,FALSE)}. */
  
int SPFunction_StdCheckOrthonization
  ( Basis bas, 
    SPFunction *w, 
    Triangulation *tri, 
    int maxTests, 
    double maxErr,
    int maxBad,
    bool_t verbose
  );
  /* Checks whether {<bas[i]|bas[j]> = (i == j)}, for all {i,j}
    with {i<=j}, where {<f,g> = SPFunction_Dot(f,NULL,g,NULL,w,tri,FALSE)}. 
    
    Returns the number of pairs that failed the test.
    The parameters {maxTests}, {maxErr}, and {maxBad} have
    the same meaning as for {SPFunction_GenCheckOrthonization}. */

#define Basis_new SPFunctionRef_vec_new
#define Basis_expand SPFunctionRef_vec_expand
#define Basis_trim SPFunctionRef_vec_trim

#endif
