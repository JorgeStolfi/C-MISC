#define PROG_NAME "compute_interpolators"
#define PROG_DESC "computes coeffs of interpolating splines"
#define PROG_VERS "1.0"

/* Last edited on 2024-12-21 11:56:46 by stolfi */
/* Created on 2007-07-11 by J. Stolfi, UNICAMP */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include <affirm.h>
#include <bool.h>
#include <jsmath.h>
#include <rn.h>
#include <rmxn.h>
#include <gausol_solve.h>

typedef struct spline_t
  { int m;        /* Total support width. */
    int g;        /* Degree of each piece. */
    double **a;   /* Coeffcients of each piece. */
  } spline_t;
  /* A spline consisting of {m} pieces, each a polynomial of 
    degree {g} whose domain is an unit interval with integer length.
    
    The element {a[k][i]} is the exponent of {z^i} in piece {k}, where {z}
    is the relative position inside interval {k}, between 0 and 1. */

spline_t *spline_alloc(int m, int g);
  /* Allocates a spline with {m >= 1} pieces of degree {g >= 0}. */

void spline_free(spline_t *S);
  /* Reclaims all storage of {S}, including the descriptor {*S} itself. */

void poly_print(FILE *wr, int g, int r, int t, double h, double P[], char *v);
  /* Writes to {wr} the derivative of order {r} of a 
    polynomial of degree {g} on the variable {v}.  If {t=0}, the coefficients 
    of {v^i} is {P[i]/h^i}.  If {t>0}, the coefficient of {v^i}
    is {P[i]/h^i*falpow(i+t,t)}.*/

void spline_print(FILE *wr, char *tag, int c, double h, spline_t *S);
  /* Writes to {wr} the formulas for spline {S}, saying "hopefully of
    order {c}", assuming each cell has width {h}. Also writes the
    formulas for its derivatives from order 1 up to order {c+1}. The
    formulas are compatible with {gnuplot}; they define functions called
    "{tag}{c}(x)" for the spline and "D{r}{tag}{c}(x)" for the
    derivative of order {r}, where {g = S->g}.  When {c} is {-1},
    the letter "n" is used instead of {c} in the function names.
    Also prints the interpolation weights as comments. */

void spline_print_pieces(FILE *wr, char *dtag, int r, spline_t *S);
  /* Writes to {wr} the formulas for the pieces of the derivative of order {r}
    of the spline {S}. The formulas are compatible with {gnuplot}; they define
    functions "{dtag}x{k}(z)" for {k} in {0..S->m-1}. The argument {z}
    is the relative position of the argument in the cell,
    in {[0_1]}. */

void spline_print_whole(FILE *wr, char *dtag, int m, int r, double h);
  /* Writes to {wr} the formula for the derivative of order {r} of the
    spline {S}, in terms of its pieces, assuming each cell has width
    {h}. The formula is compatible with {gnuplot}; it defines a function
    "{dtag}(x)", and assumes that the pieces are called "{dtag}x{k}" for
    {k} in {0..m-1}. */

void spline_print_weights(FILE *wr, char *wtag, spline_t *S, double h);
  /* Writes to {wr} the weights for interpolation using the 
    spline {S}, assuming each cell has width {h}. The weights
    are labeled "{wtag}[{k}]".  If {h==1}, then {k} ranges {0..S->m-1}.
    If {h==0.5}, then {k} ranges in {0..S->m/2-1}, and separate
    definitions are given for {z < 0.5} and {z >= 0.5}. */

void system_print(FILE *wr, int ne, int nv, int g, double M[], double b[]);
  /* Prints the system {M*z = b} to {wr}, assumed to habe
    {ne} equations and {nv} unknowns.  Assumes {M} mas {ne*nv} elements,
    linearized by rows, and {b} has {ne} elements.
    
    Expects {nv} to be a multiple of {g+1}.  Prints each unknown as
    as "z{k}{t}" standing for {z[k*(g+1)+t]}. */

spline_t *compute_B_spline_pulse(int c);
  /* A /B-spline pulse of order {c}/ is a non-negative spline {F} of degree
    {g=c+1}, that is continuous to order {c} and whose copies translated
    by integer displacements can be combined to reproduce any polynomial of degree
    {g}.
    
    This procedure computes and prints all B-splines pulses up to a given order {c},
    and returns the B-spline pulse of order {c}. */

spline_t *compute_I_spline_pulse(int c);
  /* An /I-spline pulse/ or order {c} is a spline that is continuous to
    order {c}, whose copies translated by *even* integer displacements
    can can be combined to reproduce any polynomial of degree {c+1}, and
    which is interpolating at even integers -- that is, {F(0) = 1}, and
    {F(x) = 0} for even non-zero integer {x}.
    
    Computes and prints an I-spline pulse with order {c}, and minimal
    support width. Returns NULL if it cannot not find such a spline. */ 
      
spline_t *compute_hspline_pulse(int c);
  /* An /I-spline pulse/ or order {c} is a spline that is continuous to
    order {c}, whose copies translated by *even* integer displacements
    can can be combined to reproduce any polynomial of degree {c+1}, and
    which is interpolating at even integers -- that is, {F(0) = 1}, and
    {F(x) = 0} for even non-zero integer {x}.
    
    Computes and prints an I-spline pulse with order {c}, and minimal
    support width. Returns NULL if it cannot not find such a spline. */ 
      
spline_t *compute_O_spline_pulse(int c);
  /* An /O-spline pulse/ or order {c} is a
    spline that is continuous to order {c}, whose
    copies translated by integer displacements can can be
    combined to reproduce any polynomial of degree {c+1}, and which is
    interpolating at integers -- that is, {F(0) = 1}, and {F(x) =
    0} for non-zero integer {x}.
    
    Computes and prints an O-spline pulse with order {c},
    and minimal support width. */ 
      
spline_t *compute_I_spline_pulse_specific(int c, int g, int w);
  /* Computes (does not print) an I-spline pulse of order {c} and degree {g},
     with {m=2*w} pieces.  Returns NULL if it cannot not find
    such a spline. */

spline_t *compute_O_spline_pulse_specific(int c, int g, int w);
  /* Computes (does not print) an O-spline pulse of order {c} and degree {g},
     with {m=2*w} pieces.  Returns NULL if it cannot not find
    such a spline. */

void poly_mirror(int g, double p[], double q[]);
  /* Stores into {q[0..g]} the polynomial {p[0..g]}, with its argument {z} replaced
     by {1-z}. */

void test_poly_mirror(void);
int main(int argc, char **argv);

int main(int argc, char **argv)
  {
    /* Compute an print all B-spline pulses up to {c=3}: */
    spline_t *S = compute_B_spline_pulse(3);
    spline_free(S);
    
    /* Compute and print all I-spline pulses up to {c = 3}: */
    int cmin = -1;
    int cmax = 2;  /* !!! For now. */
    int c;
    for (c = cmin; c <= cmax; c++)
      { spline_t *S = NULL;
        /* I-spline (with knots at half-integers): */
        S = compute_I_spline_pulse(c);
        spline_free(S);
        if (c >= 0)
          { /* O-spline (with knots at integers): */
            S = compute_O_spline_pulse(c);
            spline_free(S);
          }
      }
    return 0;
  }

spline_t *spline_alloc(int m, int g)
  { spline_t *S = notnull(malloc(sizeof(spline_t)), "no mem");
    S->m = m;
    S->g = g;
    S->a = notnull(malloc(m*sizeof(double *)), "no mem");
    int k;
    for (k = 0; k < m; k++) { S->a[k] = notnull(malloc((g+1)*sizeof(double)), "no mem"); }
    return S;
  }

void spline_free(spline_t *S)
  { int k;
    for (k = 0; k < S->m; k++) { free(S->a[k]); }
    free(S->a);
    free(S);
  }

void poly_print(FILE *wr, int g, int r, int t, double h, double P[], char *v)
  { 
    assert(g >= 0);
    assert(r >= 0);
    
    if (g < r)
      { fprintf(wr, "0"); }
    else if (r > 0)
      { poly_print(wr, g-1, r-1, t+1, h, &(P[1]), v); }
    else
      { 
        auto void doprint(int k, double A);
          /* Prints the polynomial from term {k} on, mutiplied by {A} and divided by {z^k}. */
        
        auto void doprint(int k, double A)
          { if (k < g)
              { fprintf(wr, "(");
                doprint(k+1, A/h);
                fprintf(wr, ")*%s", v);
              }
            fprintf(wr, "%+.6f", A*P[k]*falpow(k+t,t));
          }
          
        doprint(0, 1.0);
      }
  }

void spline_print(FILE *wr, char *tag, int c, double h, spline_t *S)
  {
    fprintf(wr, "# %s-spline with %d pieces of degree %d, hopefully of order %d\n", tag, S->m, S->g, c);
    int r;
    char *ctag = NULL;
    char *ctag = jsprintf((c < 0 ? "%c" : "%d"), (c < 0 ? 'n' : c));
    char *ptag = NULL;
    char *ptag = jsprintf("%s%s", tag, ctag);
    for (r = 0; r <= c+1; r++)
      { 
        char *dtag = NULL;
        char *dtag = jsprintf("D%d%s", r, ptag);
        /* Print the pieces: */ 
        spline_print_pieces(wr, dtag, r, S);
        spline_print_whole(wr, dtag, S->m, r, h);
        if (r == 0) { fprintf(wr, "  %s(x) = %s(x)\n", ptag, dtag); }
        free(dtag);
      }
    free(ctag);
    free(ptag);
    char *wtag = NULL;
    char *wtag = jsprintf((c < 0 ? "w%s%c" : "w%s%d"), tag, (c < 0 ? 'n' : c));
    spline_print_weights(wr, wtag, S, h);
    free(wtag);
  }

void spline_print_pieces(FILE *wr, char *dtag, int r, spline_t *S)
  { int k;
    for (k = 0; k < S->m; k++)
      { fprintf(wr, "  %sx%d(z) = ", dtag, k);
        poly_print(wr, S->g, r, 0, 1.0, S->a[k], "z");
        fprintf(wr, "\n");
      }
  }

void spline_print_whole(FILE *wr, char *dtag, int m, int r, double h)
  {
    /* Write an auxiliary function without {h}-scaling: */
    double xlo = -0.5*((double)m);
    double xhi = +0.5*((double)m);
    fprintf(wr, "  %sR(x) = ", dtag);
    fprintf(wr, "((x < %+.1f)||(x > %+.1f) ? 0 : ", xlo, xhi);

    auto void print_tree(int i, int j);
      /* Prints the spline formula for {x} in cells {i..j}. */

    void print_tree(int i, int j)
      { 
        if (i == j)
          { fprintf(wr, "%sx%d(x%+.1f)", dtag, i, -(xlo+i)); }
        else
          { int t = (i+j)/2;
            fprintf(wr, "(x < %+.1f ? ", xlo+t+1);
            print_tree(i,t);
            fprintf(wr, " : ");
            print_tree(t+1,j);
            fprintf(wr, ")");
          }
      }

    print_tree(0,m-1);

    fprintf(wr, ")");
    fprintf(wr, "\n");
    /* Write the function with {h}-scaling: */
    if (h == 1)
      { fprintf(wr, "  %s(x) = %sR(x)\n", dtag, dtag); }
    else if (r == 0)
      { fprintf(wr, "  %s(x) = %sR(x/%.1f)\n", dtag, dtag, h); }
    else
      { fprintf(wr, "  %s(x) = %f*%sR(x/%.1f)\n", dtag, 1/pow(h,r), dtag, h); }
  }        

void spline_print_weights(FILE *wr, char *wtag, spline_t *S, double h)
  {
    int i;
    if (h == 1.0)
      { 
        int nw = S->m;
        for (i = 0; i < nw; i++)
          { fprintf(wr, "# %s[%2d] = ", wtag, i);
            poly_print(wr, S->g, 0, 0, h, S->a[S->m - i - 1], "fz");
            fprintf(wr, "\n");
          }
      }
    else if (h == 0.5)
      { 
        assert((S->m % 2) == 0);
        int nw = S->m/2;
        fprintf(wr, "# if (z < 0.5) {\n");
        for (i = 0; i < nw; i++)
          { fprintf(wr, "#   %s[%2d] = ", wtag, i);
            poly_print(wr, S->g, 0, 0, h, S->a[S->m - 2*i - 2], "fz");
            fprintf(wr, "\n");
          }
        fprintf(wr, "# } else {\n");
        /* Assumes the spline is symmetric: */
        for (i = 0; i < nw; i++)
          { fprintf(wr, "#   %s[%2d] = ", wtag, i);
            poly_print(wr, S->g, 0, 0, h, S->a[2*i], "gz");
            fprintf(wr, "\n");
          }
        fprintf(wr, "# }\n");
      }
    else
      { 
        affirm(FALSE, "invalid {h}");
      }
  }

void system_print(FILE *wr, int ne, int nv, int g, double M[], double b[])
  {
    assert(nv % (g+1) == 0);
    int m = nv/(g+1);
    fprintf(wr, "system  ne = %d  nv = %d  g = %d  m = %d\n", ne, nv, g, m);
    int k, i, t;
    for (i = 0; i < ne; i++)
      { double *Mi = &(M[i*nv]);
        double *bi = &(b[i]);
        for (k = 0; k < m; k++)
          { fprintf(wr, " ");
            for (t = 0; t <= g; t++)
              { int j = k*(g+1) + t;
                if (Mi[j] != 0)
                  { fprintf(wr, " %+6.2f*z%d%d", Mi[j], k, t); }
                else
                  { fprintf(wr, " %10s", ""); }
              }
          }
        fprintf(wr, " = %+6.2f", *bi);
        fprintf(wr, "\n");
      }
     fprintf(wr, "\n");
  }
  
spline_t *compute_B_spline_pulse(int c)
  { 
    /* When {g} is zero the spline is the rectangular pulse, 
      with {a[0,0]=1}. The spline for any other degree {g>0} is
      obtained by constructing a spline {F} of degree {g-1},
      subtracting from {F} a copy of {F} shifted by {+1},
      and integrating the result */

    int g = c+1;
    int m = g+1;

    spline_t *S = spline_alloc(m,g);

    int i, k;
    if (g == 0)
      { S->a[0][0] = 1; }
    else
      { spline_t *F = compute_B_spline_pulse(c-1);
        assert(F->m == m-1);
        assert(F->g == g-1);
        /* Subtract from {F} a shifted copy of {F}: */
        for (k = 0; k < m; k++)
          { for (i = 0; i < g; i++) 
              { S->a[k][i] = (k==(m-1) ? 0 : F->a[k][i]) - (k==0 ? 0 : F->a[k-1][i]); }
            S->a[k][g] = 0;
          }
        /* Now integrate {S} from the start of the first interval: */
        double C = 0;  /* The initial value of the integral in the current interval. */
        for (k = 0; k < m; k++)
          { /* Integrate the monomials from {z^0} to {z^{g-1}}, evaluate at 1: */
            double T = 0;
            for (i = g; i > 0; i--) { S->a[k][i] = S->a[k][i-1]/i;  T += S->a[k][i]; }
            /* Add the integration constant: */
            S->a[k][0] = C;
            /* Set {C} to the value at the end of the interval: */
            C += T;
          }
        spline_free(F);
      }
    /* Print: */
    spline_print(stdout, "B", c, 1.0, S);
    return S;
  }

spline_t *compute_I_spline_pulse(int c)
  { 
    /* The interpolating spline is supported on some even number {m = 2*w}
      of intervals, to be determined.  */

    fprintf(stderr, "--- compute_I_spline_pulse(%d) ---------------------------------\n", c); 
    spline_t *S = NULL;
    int w;
    for (w = c+2; w <= 2*c+4; w++)
      { int gmin = c+1; /* !!! For now. */
        int gmax = c+1; /* !!! For now. */
        int g;
        for (g = gmin; g <= gmax; g++)
          { S = compute_I_spline_pulse_specific(c,g,w);
            if (S != NULL) { break; }
          }
        if (S != NULL) { break; }
      }
    if (S == NULL) { fprintf(stderr, "** FAILED\n"); assert(FALSE); }
    spline_print(stdout, "I", c, 0.5, S);
    fprintf(stderr, "------------------------------------------------------------------\n"); 
    return S;
  }

spline_t *compute_O_spline_pulse(int c)
  { 
    /* The interpolating spline is supported on some even number {m = 2*w}
      of intervals, to be determined.  */
    demand(c >= 0, "invalid order");

    fprintf(stderr, "--- compute_O_spline_pulse(%d) ---------------------------------\n", c); 
    spline_t *S = NULL;
    int w;
    for (w = c+2; w <= c+5; w++)
      { int gmin = c+1;   /* !!! For now. */
        int gmax = 2*c+2; /* !!! For now. */
        int g;
        for (g = gmin; g <= gmax; g++)
          { S = compute_O_spline_pulse_specific(c,g,w);
            if (S != NULL) { break; }
          }
        if (S != NULL) { break; }
      }
    if (S == NULL) { fprintf(stderr, "** FAILED\n"); assert(FALSE); }
    spline_print(stdout, "O", c, 1.0, S);
    fprintf(stderr, "------------------------------------------------------------------\n"); 
    return S;
  }

/* CONSTRAINT EQUATIONS

  Used internally by {compute_I_spline_pulse_specific} 
  and {compute_O_spline_pulse_specific} to build the linear system.
  They all assume that {Mi[0..nv-1]} and {*b} have been set to zero. */

void iO_spline_constraint_cont_at_knot(int g, int w, int x, int k, double Mi[], double *bi);
  /* Stores into {Mi[0..nv-1]} and {*b} the equation for the constraint that the 
    {k}th derivative of the piece with domain {[x_x+1]},
    evaluated at {z=0}, must equal the same derivative 
    of the piece with domain {[x-1_x]} evaluated at {z=1}. 
    Special cases for {x=0} or {x=w}. */

void iO_spline_constraint_interp_at_knot(int g, int w, int x, double Mi[], double *bi);
  /* Stores into {Mi[0..nv-1]} and {*b} the equation for the constraint
    that the spline be interpolating at {x}, that is equal to 1 if {x=0}
    and 0 otherwise. Requires {x < w} since the condition for {x = w}
    follows from continuity. */
    
void I_spline_constraint_reconstruction(int g, int w, int k, double Mi[], double *bi);  
  /* Stores into {Mi[0..nv-1]} and {*b} the equation for the constraint
    that the {k}th derivative of the interpolation of a polynomial of
    degree {k}, evaluated at {x=0}, is the constant {2^k*k!}. Requires 
    {x < w} since the condition for {x = w} follows from continuity. */
     
void O_spline_constraint_reconstruction(int g, int w, int k, double Mi[], double *bi);  
  /* Stores into {Mi[0..nv-1]} and {*b} the equation for the constraint
    that the {k}th derivative of the interpolation of a polynomial of
    degree {k}, evaluated at {x=0}, is the constant {k!}. Requires 
    {x < w} since the condition for {x = w} follows from continuity. */
     
spline_t *iO_spline_from_system(int ne, int nv, int w, int g, double M[], double b[]);
  /* Solves the constraints system {M*z = b} and packages the solution
    {z} as a spline. Assumes that {M} has {ne*nv} coeffs, linearized by
    rows, that {b} has {ne} elements. Assumes that the solution
    {z[0..nv-1]} consists of the coefficients of {w} polynomial pieces
    of the positive half of the spline, each of degree {g}. If the
    system can't be solved, returns NULL. Also prints the system to
    {stderr}. */

spline_t *compute_I_spline_pulse_specific(int c, int g, int w)
  { 
    fprintf(stderr, "--- compute_I_spline_pulse_specific(%d,%d,%d) -----------------------\n", c, g, w); 
    demand(c >= -1, "invalid continuity order");
    
    if (g < c+1) 
      { fprintf(stderr, "degree is too low\n"); 
        return NULL;
      }
    
    int nv = w*(g+1);  /* Number of independent coeffs in symmetric pulse. */

    int ne_z = (c + 1)/2;  /* Continuity constraints at {x=0} (odd derivs only). */
    int ne_c = w*(c + 1);  /* Continuity constraints at positive knots. */
    int ne_i = (w + 1)/2;  /* Interpolation constraints (even int knots except {x==w}). */
    int ne_r = c + 1;      /* Reconstruction constraints for degree {c+1} polys. */
    
    /* Not sure whether the reconstruction constraints are independent: */
    if (ne_z + ne_c + ne_i > nv)
      { fprintf(stderr, "more constraints than unknowns");
        fprintf(stderr, "  ne_z = %d  ne_c = %d  ne_i = %d  ne_r = %d  nv = %d\n", ne_z, ne_c, ne_i, ne_r, nv);
        return NULL;
      }

    int ne = ne_z + ne_c + ne_i + ne_r; /* Total number of constraints: */
    if (ne < nv) 
      { fprintf(stderr, "system is underdetermined ne = %d nv = %d\n", ne, nv); 
        return NULL;
      }
    
    /* Assemble the linear system of constraints: */
    double *M = rmxn_alloc(ne,nv); rn_zero(ne*nv,M);
    double *b = rn_alloc(ne); rn_zero(ne,b);
    int ie = 0;  /* Number of equations already added to {M}. */
    
    /* Continuity contraints at {x=0}: */
    int x, kc;
    for (kc = 1; kc <= c; kc += 2)
      { iO_spline_constraint_cont_at_knot(g, w, 0, kc, &(M[ie*nv]), &(b[ie]));
        ie++;
      }

    /* Continuity contraints at positive knots: */
    for (x = 1; x <= w; x++)
      { for (kc = 0; kc <= c; kc++)
          { iO_spline_constraint_cont_at_knot(g, w, x, kc, &(M[ie*nv]), &(b[ie]));
            ie++;
          }
      }

    /* Interpolation constraints at even int knots except {x==w}: */
    for (x = 0; x < w; x += 2)
      { iO_spline_constraint_interp_at_knot(g, w, x, &(M[ie*nv]), &(b[ie]));
        ie++;
      }

    /* Reconstruction constraints for polys of degree {c+1}: */
    for (kc = 1; kc <= c+1; kc++)
      { I_spline_constraint_reconstruction(g, w, kc, &(M[ie*nv]), &(b[ie]));
        ie++;
      }
      
    assert(ie == ne);
    
    spline_t *S = iO_spline_from_system(ne, nv, w, g, M, b);
    
    fprintf(stderr, "------------------------------------------------------------------\n"); 
    return S;
  }

spline_t *compute_O_spline_pulse_specific(int c, int g, int w)
  { 
    fprintf(stderr, "--- compute_O_spline_pulse_specific(%d,%d,%d) -----------------------\n", c, g, w); 
    demand(c >= -1, "invalid continuity order");
    
    if (g < c+1) 
      { fprintf(stderr, "degree is too low\n"); 
        return NULL;
      }
    
    int nv = w*(g+1);  /* Number of independent coeffs in symmetric pulse. */

    int ne_z = (c + 1)/2;  /* Continuity constraints at {x=0} (odd derivs only). */
    int ne_c = w*(c + 1);  /* Continuity constraints at positive knots. */
    int ne_i = w;          /* Interpolation constraints (int knots except {x==w}). */
    int ne_r = g;          /* Reconstruction constraints for degree {g} polys. */
    
    /* Not sure whether the reconstruction constraints are independent: */
    if (ne_z + ne_c + ne_i > nv) 
      { fprintf(stderr, "more constraints than unknowns");
        fprintf(stderr, "  ne_z = %d  ne_c = %d  ne_i = %d  ne_r = %d  nv = %d\n", ne_z, ne_c, ne_i, ne_r, nv);
        return NULL;
      }

    int ne = ne_z + ne_c + ne_i + ne_r; /* Total number of constraints: */
    if (ne < nv) 
      { fprintf(stderr, "system is underdetermined ne = %d nv = %d\n", ne, nv); 
        return NULL;
      }
    
    /* Assemble the linear system of constraints: */
    double *M = rmxn_alloc(ne,nv); rn_zero(ne*nv,M);
    double *b = rn_alloc(ne); rn_zero(ne,b);
    int ie = 0;  /* Number of equations already added to {M}. */
    
    /* Continuity contraints at {x=0}: */
    int x, kc;
    for (kc = 1; kc <= c; kc += 2)
      { iO_spline_constraint_cont_at_knot(g, w, 0, kc, &(M[ie*nv]), &(b[ie]));
        ie++;
      }

    /* Continuity contraints at positive knots: */
    for (x = 1; x <= w; x++)
      { for (kc = 0; kc <= c; kc++)
          { iO_spline_constraint_cont_at_knot(g, w, x, kc, &(M[ie*nv]), &(b[ie]));
            ie++;
          }
      }

    /* Interpolation constraints at int knots except {x==w}: */
    for (x = 0; x < w; x++)
      { iO_spline_constraint_interp_at_knot(g, w, x, &(M[ie*nv]), &(b[ie]));
        ie++;
      }

    /* Reconstruction constraints for polys of degree {g}: */
    for (kc = 1; kc <= g; kc++)
      { O_spline_constraint_reconstruction(g, w, kc, &(M[ie*nv]), &(b[ie]));
        ie++;
      }
      
    assert(ie == ne);

    spline_t *S = iO_spline_from_system(ne, nv, w, g, M, b);
    
    fprintf(stderr, "------------------------------------------------------------------\n"); 
    return S;
  }

spline_t *iO_spline_from_system(int ne, int nv, int w, int g, double M[], double b[])
  {
    assert(nv == w*(g+1));
    
    /* Output the system: */
    system_print(stderr, ne, nv, g, M, b);

    /* Try to solve the system: */
    double *z = rn_alloc(nv);  /* Coeff vector: */
    if (ne > nv) 
      { fprintf(stderr, "system is overdetermined ne = %d nv = %d, trying LSQ\n", ne, nv); 
        /* Build the system {M^t*M*z = M^t*b} namely {MtM*z = Mtb}: */
        double *MtM = rmxn_alloc(nv,nv);
        double *Mtb = rn_alloc(nv);
        rmxn_tr_mul(ne, nv, nv, M, M, MtM);
        rmxn_map_row(ne, nv, b, M, Mtb);
        uint32_t r;
        gausol_solve(nv, nv, MtM, 1, Mtb, z, TRUE,TRUE, 0.0, NULL, &r);
        free(MtM); free(Mtb);
        if (r < nv)
          { fprintf(stderr, "LSQ system is not definite, r = %d\n", r); 
            free(M); free(b); free(z);
            return NULL;
          }
      }
    else
      { assert(ne == nv);
        fprintf(stderr, "solving the system with ne = nv = %d\n", ne); 
        uint32_t r;
        gausol_solve(ne, nv, M, 1, b, z, TRUE,TRUE, 0.0, NULL, &r);
        if (r < nv)
          { fprintf(stderr, "system is not definite, r = %d\n", r); 
            free(M); free(b); free(z);
            return NULL;
          }
      }
    /* Check is all eqs were satisfied: */
    fprintf(stderr, "computing the residual\n"); 
    double *e = rn_alloc(ne);
    gausol_residual(ne, nv, M, 1, b, z, e);
    int bad = 0;
    double tol = 1.0e-8;
    int ie;
    for (ie = 0; ie < ne; ie++)
      { if (fabs(e[ie]) > tol) 
          { fprintf(stderr, "  equation %d not satisfied (e = %f)\n", ie, e[ie]); bad++; }
      }
    free(e);
    
    free(M); free(b);
    
    if (bad > 0)
      { fprintf(stderr, "system solution failed, bad = %d\n", bad); 
         free(z);
         return NULL;
      }
            
    /* Create the spline: */
    fprintf(stderr, "packaging as a spline\n"); 
    int m = 2*w;
    spline_t *S = spline_alloc(m, g);
    { int k;
      for (k = 0; k < m; k++)
        { double *ak = S->a[k]; /* Coeffs of piece {k} of spline. */
          if (k < w)
            { /* Negative half of spline, take positive half and flip it: */
              double *zk = &(z[(w-1-k)*(g+1)]);
              poly_mirror(g, zk, ak);
            }
          else
            { double *zk = &(z[(k-w)*(g+1)]);
              rn_copy(g+1, zk, ak);
            }
        }
    }
    free(z);
    return S;
  }

void iO_spline_constraint_cont_at_knot(int g, int w, int x, int k, double Mi[], double *bi)
  { 
    assert((k >= 0) && (k <= g));
    assert((x >= 0) && (x <= w));
    if (x < w) 
      { Mi[x*(g+1) + k] = -falpow(k,k); }
    if (x > 0) 
      { int j; 
        for (j = k; j <= g; j++)
          { Mi[(x-1)*(g+1) + j] = falpow(j,k); }
      }
  }


void iO_spline_constraint_interp_at_knot(int g, int w, int x, double Mi[], double *bi)
  { 
    assert((x >= 0) && (x < w));
    Mi[x*(g+1)] = 1;
    if (x == 0) { (*bi) = 1; }
  }

void I_spline_constraint_reconstruction(int g, int w, int k, double Mi[], double *bi)
  {
    assert((k > 0) && (k <= g));
    int i, x;
    for (x = 1; x < w; x++)
      { if ((x%2) == 0)
          { Mi[x*(g+1) + k] = falpow(k,k)*pow(-x,k); }
        else
          { for (i = k; i <= g; i++)
              { Mi[x*(g+1)+i] = falpow(i,k)*pow(-(x+1),k); }
          }
      }
    (*bi) = falpow(k,k);
  }

void O_spline_constraint_reconstruction(int g, int w, int k, double Mi[], double *bi)
  {
    assert((k > 0) && (k <= g));
    int i, x;
    for (x = 0; x < w; x++)
      { for (i = k; i <= g; i++)
          { Mi[x*(g+1)+i] = falpow(i,k)*pow(-(x+1),k); }
        if (x > 0)
          { Mi[x*(g+1) + k] += falpow(k,k)*pow(-x,k); }
      }
    (*bi) = falpow(k,k);
  }

void poly_mirror(int g, double p[], double q[])
  { 
    q[0] = p[0];
    if (g == 0) { return; }
    poly_mirror(g-1, &(p[1]), &(q[1]));
    int i;
    for (i = 1; i<= g; i++) { q[i-1] += q[i]; q[i] = -q[i]; }
  }

void test_poly_mirror(void)
  { 
    /* Test poly_mirror: */
    fprintf(stderr, "----------------------------------------------------------------------\n");
    fprintf(stderr, "testing {poly_mirror} ...\n");
    int g = 5;
    double p[g+1], q[g+1];
    int t;
    for (t = 0; t <= g; t++)
      { rn_zero(g+1, p);
        p[t] = 1;
        fprintf(stderr, "P(z)   = ");
        poly_print(stderr, g, 0, 0, 1.0, p, "z");
        fprintf(stderr, "\n");
        poly_mirror(g, p, q);
        fprintf(stderr, "P(1-z) = ");
        poly_print(stderr, g, 0, 0, 1.0, q, "z");
        fprintf(stderr, "\n\n");
      }
    fprintf(stderr, "----------------------------------------------------------------------\n");
  }

