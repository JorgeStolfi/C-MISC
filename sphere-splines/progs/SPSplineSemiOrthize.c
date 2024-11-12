/* See SPSplineSemiOrthize.h */
/* Last edited on 2009-01-09 22:07:14 by stolfi */

#include <SPSplineSemiOrthize.h>
    
#include <SPFunction.h>
#include <SPSpline.h>
#include <SPMatrix.h>
#include <SPBasisMatrices.h>
#include <SPSys.h>
#include <SPIntegral.h>
#include <SPTriang.h>
#include <SPQuad.h>
#include <SPH3.h>

#include <affirm.h>
#include <nat.h>
#include <bool.h>
#include <filefmt.h>
#include <nget.h>
#include <fget.h>
#include <r3.h>
#include <r4.h>
#include <r3x3.h>

#include <math.h>
#include <assert.h>
#include <limits.h>
#include <values.h>
#include <stdio.h>
#include <stdlib.h>

#define T SPSpline

/* INTERNAL PROTOTYPES */

int_vec_t SPSplineSemiOrthize_FindSubElems
  ( Basis bas, 
    SPFunction *s, 
    Triangulation *tri, 
    bool_t verbose
  );
  /* Returns all indices {ix}
    such that the support of {bas[j]} is (properly)
    contained in the support of {s}. */

void SPSplineSemiOrthize_Class
  ( Basis cls, 
    Basis prv, 
    SPMatrix prvM,
    Triangulation *tri,
    Metric dot,
    int gramSchmidt, 
    int eigenfuncs,
    Metric mdot,
    bool_t verbose
  );
  /* Assumes that {cls} is a list of splines with the same support
    {X}. Makes every element {cls[i]} in {cls} orthogonal to each
    other and to every element of {prv} whose support is contained in
    {X}. Assumes that {prv} is already semi-orthogonal and contains no
    element whose support is {X} or a superset thereof.

    Also assumes that {prvM} is the lower half of the stiffness matrix
    of {prv}. */

void SPSplineSemiOrthize_AgainstPrev
  ( Basis cls,
    Basis sub,
    SPMatrix subL,
    Triangulation *tri,
    Metric dot, 
    bool_t verbose
  );
  /* Makes every element of a support class {cls} orthogonal to every
    element of {bas} whose support is properly contained in the support
    {D} of {cls}. Assumes that {bas} has been semi-orthonized already. */

SPMatrix SPSplineSemiOrthize_GetSubMatrix(SPMatrix M, int_vec_t ix);
  /* Given the lower half {M} of a basis stiffness matrix, extracts the 
    (complete) stiffness matrix {S} for those elements with indices listed
    in the {ix} vector. */

void SPSplineSemiOrthize_ComputeMatrixRows
  ( Basis bas, 
    int cini, 
    int clim,
    Triangulation *tri,
    Metric dot,
    bool_t precise,
    MatEntry_vec_t *Ments,
    int *nM,
    double *maxOrthoError,
    double *maxNormError,
    double *minGenEntry,
    bool_t verbose
  );
  /* Appends to the lower half of the stiffness matrix {M} the dot 
    products of elements {bas[cini..clim-1]} with
    all previous elements and with themselves.  More precisely, computes
    {dot(bas[i], bas[j])} for all {i} in {cini..clim-1} 
    and all {j} in {0..i}. Assumes that the current entries in {M}
    are {Ments[0..nM-1]}, and updates {nM} accordingly. Expands {Ments}
    as needed. 

    The {precise} parameter has the same meaning as in
    {SPSplineSemiOrthize_Basis}.

    If {precise==TRUE}, the procedure updates {maxOrthoError} and {maxNormError}
    as in {SPSplineSemiOrthize_Check}. In any case, it also updates
    {minGenEntry} to record the minimum absolute value of
    {dot(bas[i],bas[j])} when neither support is contained in the
    other, provided that the dot product is nonzero.  The client must
    initialize {maxOrthoError} and {maxNormError} to zero, and {minGenEntry} to plus
    infinity. */

/* IMPLEMENTATIONS */

void SPSplineSemiOrthize_Basis
  ( Basis bas, 
    Triangulation *tri,
    bool_t weirdToo,
    Metric dot,
    bool_t precise,
    int gramSchmidt, 
    int eigenfuncs,
    Metric mdot,
    bool_t verbose
  )
  { 
    /* Make sure that basis elements are sorted in a subset-compatible order: */
    if (verbose) { fprintf(stderr, "Checking basis order...\n"); }
    SPSpline_CheckBasisOrder(tri, bas, 0, bas.ne);

    /* Lower half of stiffness matrix of the semi-orthized basis, built along the way: */
    MatEntry_vec_t Ments = MatEntry_vec_new(bas.ne); /* Will grow... */
    int nM = 0; /* Entries of {M} are {Ments[0..nM-1]}; note {nM <= Ments.e}. */

    /* The procedure builds a temporary matrix {M} that is the lower half 
      (including the diagonal) of the stiffness matrix of the 
      semi-orthogonalized basis, i.e. {M[i,j] = dot(bas[i],bas[j])}
      for all {i} and all {j <= i}. In particular, {M[i,j]} should be
      zero (i.e. missing) if the support of {j} is contained
      in the support of {i}, or vice-versa. */
    
    if (verbose) { fprintf(stderr, "semi-orthonormalizing basis...\n"); }
    int cini, clim; 
    cini = 0;
    double maxOrthoError = 0.0;  /* Max {fabs(M[i,j])} when {bas[i],bas[j]} should be orthogonal. */
    double maxNormError = 0.0;  /* Max {fabs(M[i,i] - 1)} when {bas[i]} should be normalized. */
    double minGenEntry = INFINITY;  /* Min nonzero {fabs(M[i,j])} other than those above. */
    while (cini < bas.ne)
      { /* Process next support class, starting at {bas[cini]}: */
        SPFunction *bini = bas.e[cini];
        
        /* Locate the last element {bas[clim-1]} of this support class: */
        for (clim = cini+1; clim < bas.ne; clim++)
          { SPFunction *blim = bas.e[clim];
            if (! SPSpline_SameSupport(bini, blim, tri)) { break; }
          }
        affirm((cini < clim) && (clim <= bas.ne), "duh?");
          
        /* Make basis descriptor {cls} for current class {bas[cini..clim-1]}: */
        Basis cls = (Basis){ /*nel*/ clim - cini, /*el*/ &(bas.e[cini])};

        /* Report class: */
        fprintf(stderr, "  cls=[%04d..%04d]", cini, clim-1);
        if (verbose) 
          { fprintf(stderr, " %4d elems ", cls.ne);
            SPSpline_DescribeSupport(stderr, cls.e[0], tri);
          }

        if (weirdToo || (! SPSpline_WeirdSupport(cls.e[0], tri)))
	  { /* Semi-orthogonalize class. */
            /* Make descriptor for all previous elements: */
            Basis prv = (Basis){ /*nel*/ cini, /*el*/ bas.e };
            /* Make descriptor for the stiffness matrix of those elements: */
            MatEntry_vec_t Pents = (MatEntry_vec_t) { /*nel*/ nM, /*el*/ Ments.e };
            SPMatrix prvM = (SPMatrix){ /*rows*/ cini, /*cols*/ cini, /*ents*/ Pents };
            /* Make all elements of {cls} orthogonal to {prv} and to themselves: */
            SPSplineSemiOrthize_Class
             ( cls, prv, prvM, tri, dot, gramSchmidt, eigenfuncs, mdot, verbose );
          }

        /* Make the class elements `mostly positive' and normalized: */ 
        if (verbose) { fprintf(stderr, "\n    [++]"); }
        int i; 
        for (i = cini; i < clim; i++)
          { if (verbose) { fprintf(stderr, "+"); }
            SPFunction_GenMakePositive(bas.e[i], dot);
            if (verbose) { fprintf(stderr, "|"); }
            SPFunction_GenNormalize(bas.e[i], dot);
          }

	/* Add to {M} the dot products of the class elements and previous elements: */
        if (verbose) { fprintf(stderr, "\n    [+M]"); }
        SPSplineSemiOrthize_ComputeMatrixRows
          ( bas, cini, clim, 
            tri, dot, precise, 
            &Ments, &nM, 
            &maxOrthoError, &maxNormError, &minGenEntry, 
            verbose 
          );

        fprintf(stderr, "\n");
        /* Prepare for next support class: */
        cini = clim;
      }

    fprintf(stderr, "\n");
    if (precise)
      { fprintf(stderr, "max orthogonalization error = %8.2e\n", maxOrthoError);
        fprintf(stderr, "max normalization error =     %8.2e\n", maxNormError);
      }
    fprintf(stderr, "min non-ortho dot prod =      %8.2e\n", minGenEntry);

    free(Ments.e);
  }

void SPSplineSemiOrthize_Class
  (
    Basis cls, 
    Basis prv, 
    SPMatrix prvM,
    Triangulation *tri,
    Metric dot,
    int gramSchmidt, 
    int eigenfuncs,
    Metric mdot,
    bool_t verbose
  )
  { int ngs = gramSchmidt, nef = eigenfuncs;
    int niter = (ngs > nef ? ngs : nef), iter;
    Basis sub;
    SPMatrix subL;
    bool_t weird = SPSpline_WeirdSupport(cls.e[0], tri);
    if (weird)
      { /* We will not orthogoalize against previous elements: */
        sub = SPFunctionRef_vec_new(0);
        subL = SPMatrix_Null(0,0);
      }
    else
      { /* Gather indices of elems in {prv} whose supp is proper subset of {supp(cls)}: */
        int_vec_t subIx = SPSplineSemiOrthize_FindSubElems(prv, cls.e[0], tri, verbose);
        /* Build basis {sub} with those elements: */
        sub = SPFunction_GetSubBasis(prv, subIx);
        /* Extract the (complete) dot product matrix {subM} of {sub}: */
        SPMatrix subM = SPSplineSemiOrthize_GetSubMatrix(prvM, subIx);
        /* Get the Cholesky factor {subL} of {prvM}: */
        SPMatrix_Cholesky(subM, 0.0, &subL);
        free(subIx.e);  free(subM.ents.e);
      }

    iter = 1;
    for (iter = 1; iter <= niter; iter++)
      { if (! weird)
          { /* Orthonormalize {cls} against all properly contained classes: */
            if (verbose) { fprintf(stderr, "\n    [+P]"); }
            SPSplineSemiOrthize_AgainstPrev(cls, sub, subL, tri, dot, verbose);
	  }
        if (cls.ne > 1)
          { /* Now orthogonalize the class {cls}: */
            if (iter <= ngs)
              { if (verbose) { fprintf(stderr, "\n    [+G]"); }
                SPFunction_GenOrthonizeBasis(cls, dot, verbose);
              }
            if (iter <= nef)
              { /* Eigenfunction orthogonalization: */
                double_vec_t ev = double_vec_new(cls.ne);
                if (verbose) { fprintf(stderr, "\n    [+E]"); }
                SPFunction_GenEigenFuncBasis(cls, dot, mdot, ev, verbose);
                if (verbose) 
                  { int i;
                    fprintf(stderr, " ev = (");
                    for (i = 0; i < cls.ne; i++)
                      { fprintf(stderr, " %10.6f", ev.e[i]); }
                    fprintf(stderr, " )");
                  }
                free(ev.e);
              }
          }
      }
    if (sub.ne > 0) { free(sub.e); free(subL.ents.e); }
  }

int_vec_t SPSplineSemiOrthize_FindSubElems
  ( Basis bas, 
    SPFunction *s, 
    Triangulation *tri, 
    bool_t verbose
  )
  {
    affirm(SPSpline_Cast(s) != NULL, "not PW function");
    affirm(! SPSpline_WeirdSupport(s, tri), "weird class");

    int_vec_t subIx = int_vec_new(10);
    int subN = 0; /* Elements found so far are {subIx[0..subN-1]}. */
    int j; 
    for (j = 0; j < bas.ne; j++)
      { SPFunction *bj = bas.e[j];
        affirm(! SPSpline_WeirdSupport(bj, tri), "weird basis");
        if (SPSpline_StrictSubSupport(bj, s, tri))
          { SPSpline *bpwj = SPSpline_Cast(bj);
            affirm(bpwj != NULL, "unsorted basis");
            if (verbose) { fprintf(stderr, "%s%d", (subN > 0 ? "," : "\n  sub=["), j); }
            /* Save in {subIx}: */
            int_vec_expand(&subIx, subN);
            subIx.e[subN] = j; subN++;
          }
      }
    int_vec_trim(&subIx, subN);
    if (verbose && (subN > 0)) { fprintf(stderr, "](%d)", subN); }
    return subIx;
  }

SPMatrix SPSplineSemiOrthize_GetSubMatrix(SPMatrix M, int_vec_t ix)
  {
    /* Sub-matrix (complete) being built: */
    MatEntry_vec_t Sents = MatEntry_vec_new(ix.ne); /* Will grow... */
    int nS = 0; /* Entries of {S} are {Sents[0..nS-1]}; note {nS <= Sents.e}. */
    
    int ki, kj;
    for (ki = 0; ki < ix.ne; ki++)
      { int i = ix.e[ki];
        for (kj = 0; kj < ix.ne; kj++)
	  { int j = ix.e[kj];
            double val;
            if (j <= i)
              { /* Get element from main matrix: */
                val = SPMatrix_GetValue(M, i, j);
              }
            else
              { /* Element has been extracted already, with transposed indices: */
                nat_t kx = SPMatrix_Find(Sents.e, nS, j, i);
                val = (kx < nS ? Sents.e[kx].va : 0.0); 
              }
	    if (val != 0.0)
              { MatEntry_vec_expand(&Sents, nS);
                Sents.e[nS] = (MatEntry){/*row*/ ki, /*col*/ kj, /*val*/ val };
                nS++;
              }
	  }
      }
    MatEntry_vec_trim(&Sents, nS);
    return (SPMatrix){ /*rows*/ ix.ne, /*cols*/ ix.ne, /*ents*/ Sents };
  }

void SPSplineSemiOrthize_AgainstPrev
  ( Basis cls,
    Basis sub,
    SPMatrix subL,
    Triangulation *tri,
    Metric dot,
    bool_t verbose
  )
  { /* Normalize elements of {cls}, should be good for numerics: */
    int i;
    if (verbose) { fprintf(stderr, "n"); }
    for (i = 0; i < cls.ne; i++) { SPFunction_GenNormalize(cls.e[i], dot); }
    
    /* Now make each class element {cls[i]} orthogonal to all
      {sub[k]}.  Since the elements of {sub} may not be orthogonal
      among themselves, we must solve for each {i} the linear system
      {B x = y} where {y[s] = dot(cls[i],sub[s])} and {B[r,s] =
      dot(sub[r],sub[s])}.  The system is solved through Cholesky
      factorization of {B}.  Then we make {cls[i] = cls[i] - SUM{
      x[j]*sub[j] : j }}.
  
      The matrix {B} may be large but hopefully sparse, so we use the 
      SPMatrix package rather than plain arrays.
    */

    double_vec_t y = double_vec_new(sub.ne); /* Right-hand side of system. */
    double_vec_t x = double_vec_new(sub.ne); /* Coefficients of projection. */ 
    double_vec_t t = double_vec_new(sub.ne); /* Work area for {CholeskySolve}. */ 
    double_vec_t subNorm = double_vec_new(sub.ne); /* Norm of {sub} vectors. */
    
    if (verbose) { fprintf(stderr, "m"); }
    for (i = 0; i < sub.ne; i++)
      { subNorm.e[i] = sqrt(fabs(dot(sub.e[i],sub.e[i]))); }
    if (verbose) { fprintf(stderr, "s"); }
    for (i = 0; i < cls.ne; i++)
      { SPFunction *ci = cls.e[i];
        double ciMag;
        int maxIter = 3, iter;
        for (iter = 1; iter <= maxIter; iter++)
          { double dciMax = 0.0;
            int j;
            /* Find coeffs {x[j]} of projection of {ci} on space of {sub}: */
            SPBasisMatrices_RecomputeVectorGen(ci, sub, dot, y, FALSE);
            fprintf(stderr, "x");
            SPSys_CholeskySolve(subL, y, t, x, FALSE);
            /* Subtract {x[j]*sub[j]} from {cls[i]}: */
            for (j = 0; j < sub.ne; j++)
              { double xj = x.e[j];
                if (xj != 0.0)
                  { ci->m->add(ci, -xj, sub.e[j]);
                    { double dci = fabs(xj)*subNorm.e[j]; 
                      if (dci > dciMax) { dciMax = dci; }
                    }
                  }
              }
            /* if (verbose) { fprintf(stderr, "(%8.2e)", xMax); } */
            if (dciMax == 0.0) { break; }
            if (iter == 1) { ciMag = sqrt(dot(ci,ci)); }
            if (dciMax <= 1.0e-13*ciMag) { break; }
            if (iter >= maxIter)
              { fprintf(stderr, "!");
                if (verbose)
                  { fprintf(stderr, "(%8.2e:%8.2e)", dciMax, ciMag); }
              }
          }
      }
    if (verbose) { fprintf(stderr, " "); }
    free(subNorm.e);
    free(x.e); free(t.e); free(y.e);
  }
  
void SPSplineSemiOrthize_ComputeMatrixRows
  ( Basis bas, 
    int cini, 
    int clim,
    Triangulation *tri,
    Metric dot,
    bool_t precise,
    MatEntry_vec_t *Ments,
    int *nM,
    double *maxOrthoError,
    double *maxNormError,
    double *minGenEntry,
    bool_t verbose
  )
  { int i, j;
    for (i = cini; i < clim; i++)
      { if (verbose) { fprintf(stderr, "m"); }
        SPFunction *bi = bas.e[i];
        for (j = 0; j <= i; j++)
          { SPFunction *bj = bas.e[j];
            bool_t nested = SPSpline_SubSupport(bj, bi, tri);
            double dprod;
            if (precise)
              { dprod = dot(bi, bj); }
            else 
              { if (j >= cini)
                  { dprod = (double)(j == i); }
                else if (nested)
                  { dprod = 0; }
                else
                  { dprod = dot(bi, bj); }
	      }
            /* Update errors: */
            if (i == j)
              { double dn = fabs(dprod - 1.0);
                if (dn > (*maxNormError)) { (*maxNormError) = dn; }
              }
            else if (nested)
              { double dd = fabs(dprod);
 	        if (dd > (*maxOrthoError)){ (*maxOrthoError) = dd; }
	      }
            else if (dprod != 0.0)
              { double dd = fabs(dprod);
 	        if (dd < (*minGenEntry)){ (*minGenEntry) = dd; }
	      }
    	    /* Store in array: */
	    if (dprod != 0.0)
              { MatEntry_vec_expand(Ments, (*nM));
                Ments->e[(*nM)] = (MatEntry){/*row*/ i, /*col*/ j, /*val*/ dprod };
                (*nM)++;
              }
          }
      }
  }

void SPSplineSemiOrthize_Check
  ( Basis bas, 
    int ini, int lim, 
    Triangulation *tri,
    bool_t weirdToo,
    Metric dot,
    double *maxOrthoError,
    double *maxNormError,
    double maxError,
    bool_t verbose
  )
  { 
    int i, j;
    (*maxOrthoError) = 0.0; double worstOrtho = 0.0;
    (*maxNormError) = 0.0;  double worstNorm = 1.0;
    if (verbose) { fprintf(stderr, "checking semi-orthonormality...\n"); }
    for (i = ini; i < lim; i++)
      { SPFunction *bi = bas.e[i];
        int nChk = 0;
        if ((! weirdToo) && SPSpline_WeirdSupport(bi, tri)) { break; }
        if (verbose) { fprintf(stderr, "[%04d]", i); }
        for (j = i; j < lim; j++)
          { SPFunction *bj = bas.e[j];
            bool_t check; char *rel;
             
            if ((! weirdToo) && SPSpline_WeirdSupport(bj, tri)) { break; }
            /* Check consistency of domain comparisons: */
            { bool_t sub = SPSpline_StrictSubSupport(bi, bj, tri);
              bool_t eql = SPSpline_SameSupport(bi, bj, tri);
              bool_t sup = SPSpline_StrictSubSupport(bj, bi, tri); 
              bool_t subeql = SPSpline_SubSupport(bi, bj, tri);
              bool_t supeql = SPSpline_SubSupport(bj, bi, tri);
              affirm(subeql == (sub | eql), "support cmp inconsistent (1)");
              affirm(supeql == (sup | eql), "support cmp inconsistent (2)");
              affirm(sub + eql + sup <= 1, "support cmp inconsistent (3)");
              check = sub | eql | sup;
              rel = (eql ? "=" : (sub ? "<" : (sup ? ">" : "?")));
            }

            if (check)
              { double dij = dot(bi, bj), err;
                if (verbose) 
                  { fprintf(stderr, "%s%s%d", (nChk == 0 ? " [" : ","), rel, j); }
                if (i == j)
                  { /* Check normalization: */
                    err = fabs(dij - 1.0);
                    if (err > (*maxNormError)) 
                      { (*maxNormError) = err; worstNorm = dij; }
                  }
                else
                  { /* Check orthogonality: */
                    err = fabs(dij);
                    if (err > (*maxOrthoError))
                      { (*maxOrthoError) = err; worstOrtho = dij; }
                  }
                if (err > maxError)
                  { fprintf
                      ( stderr, " ** <%04d|%04d> = %14.8e, aborted\n", i, j, dij );
		    return;
                  }
                nChk++;
              }
          }
        if (verbose) { fprintf(stderr, "](%d)\n", nChk); }
      }
    fprintf(stderr, "worst <b[i]|b[i]> = %14.8e err = %8.2e\n", worstNorm, (*maxNormError));
    fprintf(stderr, "worst <b[i]|b[j]> = %14.8e err = %8.2e\n", worstOrtho, (*maxOrthoError));
  }
