/* Last edited on 2005-10-27 01:45:22 by stolfi */

/* SPHessTest */

r3_t PickOrthoDir(S2Point *p)
  { int i, imax = 0, imin = 0;
    r3_t u;
    for (i = 0; i < 3; i++)
      { double pi = fabs(p->c[i]); 
        if (pi > fabs(p->c[imax])) { imax = i; }
        if (pi <= fabs(p->c[imin])) { imin = i; }
      }
    affirm(imax != imin, "huh?");
    u = (r3_t){{0.0, 0.0, 0.0}};
    u.c[imax] = -(p->c[imin]);
    u.c[imin] = p->c[imax];
    r3_dir(&u, &u);
    return u;
  }

/* SPSpline.h */

r3_t r3_pick_ortho_vec_t(r3_t *u)
  { r3_t v = (r3_t){{0.0, 0.0, 0.0}};
    int mi = 0, i;
    double mc = fabs(u->c[0]);
    for(i = 1; i < 3; i++) 
      { double uc = fabs(u->c[i]); 
        if (uc > mc) { mi = i; mc = uc; }
      }
    { int xi = (mi + 1) % 3;
      v.c[xi] = u->c[mi];
      v.c[mi] = - u->c[xi];
    }
    return v;
  }

void SPSpline_SemiOrthonizeBasis
  ( Basis bas, 
    Triangulation *tri,
    bool_t weirdToo,
    Metric dot,
    int gramSchmidt,
    int eigenfuncs,
    Metric mdot,
    bool_t verbose
  );
  /* Applies a partial orthonormalization procedure to the elements of
    basis {bas}.
    
    The elements must have been previously sorted in an order
    compatible with domain containment.
    
    The elements are processed in increasing order; each normal
    element {bas[i]}, with support {D}, is made {dot}-orthogonal to all
    elements {bas[j]} whose support is strictly contained in {D} (and
    which therefore must have been already processed).
    
    Then, the elements in each support class are made {dot}-orthogonal to
    each other:
    
      If {gramSchmidt > 0}, the elements in the class are made
      {dot}-orthogonal by the Gram-Schmidt method. (This step is
      repeated {gramSchmidt} times for each class.)
    
      Them if {eigenfuncs > 0}, the elements in the class are replaced
      by an orthonormal set of stationary functions of the metric
      induced by {mdot}. (This step is repeated {eigenfuncs} times for
      each class.)
    
    Any weird elements in {bas} are moved to the end of the subset. If
    {weirdToo} is true, they are orthogonalized among themselves, as
    described above. They are never combined with normal ones.
    
    In any case, all elements (including the weird ones) are
    made `mostly positive' (with {SPFunction_GenMakePositive(f,dot)})
    and normalized with {SPFunction_GenNormalize(f,dot)}.
    
    All elements in {bas} are modified by the procedure via the
    {scale} and {add}, and thus must be {add}-compatible with
    each other. */

void SPSpline_CheckSemiOrthonization
  ( Basis bas, 
    int ini, int lim, 
    Triangulation *tri,
    bool_t weirdToo,
    Metric dot,
    double maxErr,
    bool_t verbose
  );
  /* Checks whether the elements {bas[ini..lim-1]} are
    semi-orthonormalized. Namely, whether {fabs(dot(bas[i],bas[i])) ==
    1} for every {i}, and {dot(bas[i],bas[j]) = 0} whenever {i != j}
    and the domain of {bas[i]} is a subset or superset of that of
    {bas[j]}.  Bombs out if the absolute error in any of these
    conditions is greater than {maxErr}. */

void SPSpline_SemiOrthizeAgainstPrev
  ( Basis cls,
    Basis bas,
    Triangulation *tri,
    Metric dot, 
    bool_t verbose
  );
  /* Makes every element of a support class {cls} orthogonal to every
    element of {bas} whose support is properly contained in the support
    {D} of {cls}. Assumes that {bas} has been semi-orthonized already. */

void SPSpline_SemiOrthonizeBasis
  ( Basis bas, 
    Triangulation *tri,
    bool_t weirdToo,
    Metric dot,
    int gramSchmidt, 
    int eigenfuncs,
    Metric mdot,
    bool_t verbose
  )
  { int cini, clim; 
    
    /* Gather all support classes and sort them in some subset-compatible order: */
    if (verbose) { fprintf(stderr, "Checking basis order...\n"); }
    SPSpline_CheckBasisOrder(tri, bas, 0, bas.ne);

    if (verbose) { fprintf(stderr, "semi-orthonormalizing basis...\n"); }
    cini = 0;
    while (cini < bas.ne)
      { /* Process next support class, starting at {bas[cini]}: */
        SPFunction *bini = bas.e[cini];
        bool_t weird = SPSpline_WeirdSupport(bini, tri);
        Basis cls;
        
        /* Locate the last element {bas[clim-1]} of this support class: */
        for (clim = cini+1; clim < bas.ne; clim++)
          { SPFunction *blim = bas.e[clim];
            if (! SPSpline_SameSupport(bini, blim, tri)) { break; }
          }
          
        /* Make basis descriptor {cls} for class {bas[cini..clim-1]}: */
        { int cn = clim - cini;
          affirm(cn > 0, "empty class?");
          cls = (Basis){ /*nel*/ cn, /*el*/ &(bas.e[cini])};
        }

        /* Report class: */
        fprintf(stderr, "  [%04d..%04d]", cini, clim-1);
        if (verbose) 
          { fprintf(stderr, " %4d elems ", cls.ne);
            SPSpline_DescribeSupport(stderr, cls.e[0], tri);
          }

        /* Semi-orthogonalize class (if asked to), then normalize: */
        { int ngs = gramSchmidt, nef = eigenfuncs;
          int niter = (ngs > nef ? ngs : nef), iter;
          iter = 1;
          for (iter = 1; iter <= niter; iter++)
            { if (! weird)
                { /* Orthonormalize {cls} against all properly contained classes: */
                  /* Make descriptor for previous elements: */
                  Basis prv = (Basis){ /*nel*/ cini, /*el*/ bas.e };
                  if (verbose) { fprintf(stderr, "\n    [+P]"); }
                  SPSpline_SemiOrthizeAgainstPrev(cls, prv, tri, dot, verbose);
                }
              if ((cls.ne > 1) && ((!weird) || weirdToo))
                { /* Now orthogonalize the class {cls}: */
                  if (iter <= ngs)
                    { if (verbose) { fprintf(stderr, "\n    [+G]"); }
                      SPFunction_GenOrthonizeBasis(cls, dot);
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
        }

        /* Finally, make the class elements `mostly positive' and normalized: */ 
        { int i; 
          for (i = cini; i < clim; i++)
            { SPFunction_GenMakePositive(bas.e[i], dot);
              SPFunction_GenNormalize(bas.e[i], dot);
            }
        }

        if (verbose) { fprintf(stderr, "\n"); }
        fprintf(stderr, "\n");
        /* Prepare for next support class: */
        cini = clim;
      }
    
    /* Now for the salutary paranoia: */
    SPSpline_CheckSemiOrthonization
      ( bas, 0, bas.ne, tri, weirdToo, dot, 1.0e-7, verbose);
  }

void SPSpline_SemiOrthizeAgainstPrev
  ( Basis cls,
    Basis bas,
    Triangulation *tri,
    Metric dot,
    bool_t verbose
  )
  { SPFunction *c0 = cls.e[0];
    SPSpline *cpw0 = SPSpline_Cast(c0);
    Basis sub = SPFunctionRef_vec_new(50);
        
    affirm(cpw0 != NULL, "not PW function");
    affirm(! SPSpline_WeirdSupport(c0, tri), "weird class");
    
    /* Collect the relevant subset {sub} of {bas}, consisting of those
      elements against which the class {cls} is to be made orthogonal
      --- i.e. all elements {bas[j]} whose support is (properly)
      contained in the support of {cls[0]}.
      
      However single-face elements can be processed right away, and
      omitted from {sub} --- since, by hypothesis, they are orthogonal
      to all other relevant elements. */
    { int nTot = 0; /* Tot number of elems {bas[j]} contained in {cls[0]}. */
      int nSub = 0; /* Ditto, excluding single-face elements. */
      int j; 
      for (j = 0; j < bas.ne; j++)
        { SPFunction *bj = bas.e[j];
          affirm(! SPSpline_WeirdSupport(bj, tri), "weird basis");
          if (SPSpline_StrictSubSupport(bj, c0, tri))
            { SPSpline *bpwj = SPSpline_Cast(bj);
              affirm(bpwj != NULL, "unsorted basis");
              if (verbose) { fprintf(stderr, "%s%d", (nTot > 0 ? "," : " +["), j); }
              if (bpwj->d->pd.ne == 1)
                { /* Single-face element: we can dispose of it right away. */
                  int i;
                  if (verbose) { fprintf(stderr, "*"); }
                  for (i = 0; i < cls.ne; i++)
                    { SPFunction *ci = cls.e[i];
                      affirm(SPSpline_StrictSubSupport(bj, ci, tri), "duh?");
                      SPFunction_GenOrthize(ci, bj, dot, dot(bj,bj));
                    }
                }
              else
                { /* Multi-face element: save in {idx} to orthonize later: */
                  Basis_expand(&sub, nSub);
                  sub.e[nSub] = bj; nSub++;
                }
              nTot++;
            }
        }
      SPFunctionRef_vec_trim(&sub, nSub);
      if (verbose && (nTot > 0)) { fprintf(stderr, "](%d)", nTot); }
    }
    
    /* Normalize elements of {cls}, should be good for numerics: */
    { int i;
      for (i = 0; i < cls.ne; i++) { SPFunction_GenNormalize(cls.e[i], dot); }
    }
    
    /* Now make each class element {cls[i]} orthogonal to all {sub[k]}: */
    SPFunction_GenOrthizeBasisAgainstBasis(cls, sub, dot, verbose, NULL, NULL);
    free(sub.e);
  }
  
void SPSpline_CheckSemiOrthonization
  ( Basis bas, 
    int ini, int lim, 
    Triangulation *tri,
    bool_t weirdToo,
    Metric dot,
    double maxErr,
    bool_t verbose
  )
  { int i, j;
    double diiWorst = 1.0, dijWorst = 0.0;
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
                    if (err > fabs(diiWorst - 1.0)) { diiWorst = dij; }
                  }
                else
                  { /* Check orthogonality: */
                    err = fabs(dij);
                    if (fabs(dij) > fabs(dijWorst)) { dijWorst = dij; }
                  }
                if (err > maxErr)
                  { fprintf
                      ( stderr, 
                        " ** dot error <%04d|%04d> = %14.8e\n", 
                        i, j, dij
                      );
                    affirm(FALSE, "aborted"); 
                  }
                nChk++;
              }
          }
        if (verbose) { fprintf(stderr, "](%d)\n", nChk); }
      }
    fprintf(stderr, "worst <b[i]|b[i]> = %14.8e\n", diiWorst);
    fprintf(stderr, "worst <b[i]|b[j]> = %14.8e\n", dijWorst);
  }

/* SPProcFunction.c */

    SPHarmonic_Trig(SPPF_Hr75_m, &x, &y);
    double fz = SPHarmonic_Legendre(SPPF_Hr75_d, SPPF_Hr75_m, z);
    double alpha = SPHarmonic_LegendreNormFactor(SPPF_Hr75_d, SPPF_Hr75_m);
    return s*alpha*x*fz;

/* SPapprox.h */

void SPApprox_ReadBasisCoeffs(FILE *rd, SPVector coeff);
  /* Reads, from the given file, a list of coefficients
    of some spherical function in terms of some basis. */

#define BasisCoeffs_FileFormat "2005-08-27"

void SPApprox_ReadBasisCoeffs(FILE *rd, SPVector coeff)
  { 
    int N, i;
    filefmt_read_header(rd, "BasisCoeffs", BasisCoeffs_FileFormat);
    fget_skip_formatting_chars(rd);
    N = nget_int(rd, "coeffs"); fget_eol(rd);
    affirm(N == coeff.ne, "wrong frame size");
    for (i = 0; i < N; i++)
      { coeff.e[i] = fget_double(rd);
        fget_eol(rd);
      }
    filefmt_read_footer(rd, "BasisCoeffs");
  }

void SPApprox_WriteBasisCoeffs(FILE *wr, SPVector coeff);
  /* Writes to the given file a list of coefficients
    of some spherical function in terms of some basis,
    in a format accepted by {SPApprox_ReadBasisCoeffs}. */

void SPApprox_WriteBasisCoeffs(FILE *wr, SPVector coeff)
  { int i;
    filefmt_write_header(wr, "BasisCoeffs", BasisCoeffs_FileFormat);
    fprintf(wr, "coeffs = %d\n", coeff.ne);
    for (i = 0; i < coeff.ne; i++)
      { fprintf(wr, "%22.16e\n", coeff.e[i]); }
    filefmt_write_footer(wr, "BasisCoeffs");
    fflush(wr);
  }



/* SPVector.h */

SPVector SPVector_Null(nat_t rows, nat_t cols);  
  /* Returns the null vector with given dimensions. */
  
SPVector SPVector_Scale(double alpha, SPVector A);
  /* Returns the product of {alpha} by {A} */
     
SPVector SPVector_Mix(double alpha, SPVector A, double beta, SPVector B);
  /* Returns the vector  {alpha*A + beta*B}. The two vectors must have
    the same index bounds. */
   
SPVector SPVector_RelDiff(SPVector A, SPVector B);  
  /* Returns the vector {M[i] = rdiff(A[i],B[i])}, where
    {rdiff(x,y) = (x - y)/sqrt(x^2 + y^2)/2}. This number is {0} when
    {x == y}, is {0.5} when only one of them is zero, and has extremal
    values {-1} and {-1} when {x} and {y} are equal and opposite). */
  
SPVector SPVector_BinaryOp
  ( SPVector A, SPVector B,
    int dPos,
    double func(int pos, double va, double vb),
    bool_t skip00
  );
  /* A generic elementwise binary operation that combines {A} with {B}
    shifted by {dPos}. Namely, returns the vector {C}, with the
    same size as {A}, such that {C[i] = func(i,A[i],B[i'])},
    where {i'= i-dPos}. Elements {B[i']} that
    don't exist are assumed to be zero, and those outside the bounds
    of {A} are ignored. If {skip00 == TRUE}, assumes that {C[i]} is
    zero when {A[i]} and {B[i']} are both zero. */

double SPVector_GetValue(SPVector A, nat_t pos);
  /* Returns the value of the element {A[pos]}. In particular,
     returns 0 if there is no {MatEntry} with those indices. */



/* SPMultiSolve.c */

    bool_t cholesky;     /* TRUE uses direct method based on Cholesky factzn. */
    bool_t SVD;          /* TRUE uses direct method based on SVD factzn. */
    bool_t gaussSeidel;  /* TRUE uses Gauss-Seidel iteration. */
    bool_t conjGrad;     /* TRUE uses the conjugate gradient method. */
    double omega;        /* Relaxation factor for Gauss-Seidel. */
    /* Options for solution non-linear system: */
    int maxIter;         /* Maximum iterations at this scale. */
    double absTol;       /* Max error compared to true sol, at this scale. */
    double relTol;       /* Max change in {u} coeffs between iters, at this scale. */


            if (sor.cholesky)
              { SPApprox_CholeskySolve(mr.HL, b, mr.y, aNew, TRUE); }
            else if (sor.SVD)
              { SPApprox_SVDSolve(mr.HL, mr.HD, mr.HR, b, mr.y, aNew, TRUE); }
            else if (sor.conjGrad)
              { SPApprox_ConjugateGradientSolve(mr.H, b, a, aNew, TRUE); }
            else if (sor.gaussSeidel)
              { SPApprox_GaussSeidelSolve(mr.H, b, sor.omega, 1, 0.0, 0.0, a, aNew, TRUE); }
            else
              { affirm(FALSE , "no solution method?"); }

/* SPApprox.c */

    if (cholesky)
      { SPMatrix HL;
        fprintf(stderr, "computing the Cholesky factorization...\n");
        HL = SPMatrix_Cholesky(H, 0.0); 
        free(H.ents.e);
        return HL;
      }
    else
      { return H; }

/* TIME BASIS EVALUATION */

/* In the following procedures, {TD(f,diff)} means the {diff}th time derivative
   of function {f}; {tg} is the maximum degree of a time element (3); and 
   {pg} is the maximum degree of a gauge element (1). */

double EvalTimeBasis(int u, double t, double tj, double tStep, int diff);
  /* Evaluates {TD(tau[j,u],diff)(t)}, i.e. {TD(tau[0,u],diff)(t - tj)}. */

void TimeBasisCoeffs(int k, int u, double tStep, int diff, double *tc);
  /* Returns in {tc[0..tg-diff]} the coefficients of {TD(tau[j-k,u],diff)(t)}
    for {t} in the interval {[tj-tStep _ tj]}, that is, 
    of {TD(tau[0,u],diff)(t - tj + k*tStep)}
    
    The coefficients are such that 
      {TD(tau[j-k,u],diff)(t) = SUM { tc[i]·z^i : i = 0..tg-diff}} 
    where {z = (t-tj)/tStep+1} ranges in {[0 _ 1]}. Note that the array {tc} 
    must have {tg+1} elements, even though only {tg-diff+1} are needed 
    for the result. */

double EvalTimeGauge(int v, double t, double tj, double tStep);
  /* Evaluates {pi[0,v](t)}. */

void TimeGaugeCoeffs(int v, double *pc);
  /* Returns in {pc[0..pg]} the coefficients of the time gauge element
    {pi[0,v]} in the time interval {[-1 _ 0]}.  The coefficients 
    are are such that {pi[0,v](t) = SUM { pc[i]·z^i : i = 0..1}}
    where {z = t/tStep+1} ranges in {[0 _ 1]}. */
    
double TimeBasisProduct(int k, int u, int v, double tStep, int diff, bool_t verbose);
  /* These procedures compute the time-wise scalar products 
    { <TD(tau[0,u],diff) | pi[k,v]> }, where {tau[0,u]} is the standard
    time basis element, and {pi} is the gauge basis element. Note that 
    the result is independent of the epoch {j}. */
